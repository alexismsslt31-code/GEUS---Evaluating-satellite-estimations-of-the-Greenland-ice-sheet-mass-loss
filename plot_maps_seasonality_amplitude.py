"""
Amplitude moyenne de la saisonnalité par localité, en points de grille,
avec le même habillage visuel que `subplot_anomaly_uncertainty` (mini-carte
de localisation Groenland, colorbar individuelle par case), mais sans fond
ArcticDEM.

Historique des corrections (par rapport à la V1) :
  1) Effet de moiré/grillage : la V1 reprojetait chaque produit satellite
     en EPSG:3413 "natif", puis alignait les 5 grilles entre elles avec
     interp_like(method="nearest"). Deux rééchantillonnages "au plus proche"
     enchaînés sur des grilles à résolutions/origines différentes = motif
     de moiré. Ici, chaque produit est reprojeté UNE SEULE FOIS, directement
     sur la grille cible commune (rio.reproject_match, resampling bilinéaire).
  2) Grille coincée dans un coin : la fenêtre est centrée sur le centroïde
     des stations de la localité, avec une demi-largeur fixe (pas le bbox
     min/max des stations).
  3) Mini-carte + colorbar par case : réutilisation directe de
     `add_location_inset`, tel qu'utilisé dans `plot_regional_variation`.
     Si elle vit déjà dans un module partagé de ton côté, remplace la
     définition ci-dessous par un import direct depuis ce module plutôt
     que de la dupliquer.

Le rendu utilise un pcolormesh sur la grille commune (cellules carrées
jointives, un pixel = un point de grille) : comme la grille est désormais
régulière et alignée pour les 5 produits, les cellules se touchent
proprement, sans effet de moiré.

Historique V3 (par rapport à la V2) :
  4) Généralisation de la métrique spatiale : les fonctions de calcul
     (compute_*_on_target_grid, compute_mean_*_for_locality) sont
     factorisées autour d'un `fit_fn` interchangeable, pour pouvoir tracer
     soit l'amplitude saisonnière (fit_sinus_amplitude_fast), soit la MAE
     des résidus tendance+saisonnalité retirés (fit_mae_residuals_fast),
     avec le même habillage visuel (mini-carte, colorbar, terminus PROMICE,
     etc.). Les trois fonctions de figure de haut niveau
     (plot_amplitude_map_by_locality_grid, plot_amplitude_map_greenland,
     plot_amplitude_map_at_location) acceptent désormais un paramètre
     `compute_mean_fn` (+ `metric_label`/`metric_prefix` pour les libellés
     et noms de fichier par défaut) plutôt que d'appeler en dur le calcul
     d'amplitude. Des wrappers `plot_mae_map_*` sont fournis pour retrouver
     exactement le même confort d'appel que les fonctions d'amplitude.

Historique V4 (par rapport à la V3) :
  5) Titre trop long / carte minuscule sur les grilles étroites (n_cols=1) :
     un group_label listant plusieurs satellites débordait du cadre d'une
     figure à 1-2 colonnes, et bbox_inches='tight' agrandissait alors le
     PNG final bien au-delà de la carte elle-même. Le titre est maintenant
     wrappé selon la largeur réelle de la figure (wrap_title), chaque
     sous-graphique s'agrandit quand n_cols est petit (subplot_size_inches),
     et le bandeau titre+légende est dimensionné en pouces absolus à partir
     du nombre réel de lignes du titre (ancré depuis le haut : titre
     d'abord, légende juste en dessous).
  6) Chemins de sortie Windows trop longs (MAX_PATH) : safe_output_path
     tronque + hash le nom de fichier si le chemin complet dépasse ~240
     caractères.
  7) Performance : _crop_native_before_reproject recadre la grille native
     géographique (lat/lon) à une bbox englobant la fenêtre EPSG:3413 visée
     AVANT tout appel à rio.reproject (coûteux), au lieu de reprojeter
     l'intégralité du Groenland pour ne garder ensuite qu'une fenêtre de
     quelques dizaines de km. Des print() de progression par satellite
     (ouverture/recadrage, reprojection, fit) évitent aussi l'impression que
     le script est bloqué sur les calculs longs (carte Groenland entier).

Historique V5 (par rapport à la V4) :
  8) Restriction temporelle Copernicus : le produit Copernicus Climate Data
     Store est jugé non fiable avant 2015 (cf. anomalies SEC identifiées
     précédemment sur la station SDL). _restrict_time_range coupe la série
     temporelle à 2015-01-01 pour ce produit uniquement, AVANT le recadrage
     spatial et le fit. Comme toutes les fonctions de tracé passent par
     compute_metric_on_target_grid, ce filtre s'applique uniformément à
     l'amplitude et à la MAE, sur la grille par localité, le Groenland
     entier et les points libres, sans rien dupliquer.
"""
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import rioxarray  # noqa: F401  (active l'accessor .rio)
from rasterio.enums import Resampling
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.patch import geos_to_path
from matplotlib.transforms import Affine2D
from pyproj import Transformer

import AWS_data
import interpolation_altimetry_AWS
from plot_trend_seasonality_residues import FIGURES_DIR_ROOT
from paths import ICE_MASK_DATA_DIR

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Station/locality helpers ─────────────────────────────────────────────
"""
def group_stations_by_locality(stations_dict) -- groups station names by locality (the part of the name before the first "_", e.g. "KAN_L"/"KAN_U" -> "KAN").

def get_station_coordinates(stations_dict, lat_key, lon_key) -- returns each station's (lon, lat), either directly from stations_dict or averaged from its hourly_data CSV.

def locality_centroid_xy(station_names, coords) -- EPSG:3413 centroid (x, y) of a locality, plus its stations' individually projected coordinates.
"""

# ── Map background elements ──────────────────────────────────────────────
"""
def add_location_inset(fig, ax, x_min, x_max, y_min, y_max, ...) -- small Greenland-wide location inset in the lower-left corner of ax, with the plotted extent outlined in red.

def get_land_ocean_paths(scale) -- precomputes Greenland's land/coastline outlines projected to EPSG:3413, as matplotlib Paths (the main axes are plain matplotlib Axes, not cartopy GeoAxes, so ax.add_feature can't be used directly) -- computed once per figure, reused on every subplot.

def add_land_ocean_mask(ax, land_paths, coast_paths, x_center, y_center, res, ...) -- draws the land (grey) / ocean (light blue) background and coastline under the data, translating/scaling the precomputed absolute-EPSG:3413 paths into the subplot's centre-relative (km) coordinates.

def get_ice_terminus(gpkg_path) -- loads (and caches) the PROMICE ice terminus line, reprojecting to EPSG:3413 only if the source file isn't already in that CRS.

def add_ice_terminus(ax, gdf_terminus, x_center, y_center, half_window, res, buffer_m, ...) -- overlays the PROMICE ice terminus line on the amplitude/MAE grid, in the same centre-relative coordinates as the rest of the subplot, cropped early to the window (+ a buffer) for performance. Returns legend handles (one per terminus type encountered), meant to be aggregated across localities into one shared legend.
"""

# ── Reprojection / cropping (metric grid) ────────────────────────────────
"""
def _restrict_time_range(ds, sat_name, min_date) -- restricts the time dimension for one named satellite (Copernicus Climate Data Store, judged unreliable before 2015 -- see the SDL station SEC anomalies flagged earlier in this project). No-op for any other satellite, or if the dataset has no "time" dimension. Called before any cropping/fit, so it applies uniformly to the amplitude and MAE metrics, on the per-locality grid, the Greenland-wide grid and the free-point queries.

def _to_epsg3413(ds, var_name) -- reprojects a satellite product's native grid to EPSG:3413.

def _latlon_bbox_from_epsg3413_window(x_center, y_center, buf, pad_deg) -- lat/lon bounding box covering an EPSG:3413 window, with a degree padding, used to crop the native lat/lon grid before reprojecting.

def _crop_native_before_reproject(ds, var_name, x_center, y_center, buf, pad_deg) -- crops the native geographic (lat/lon) grid to a bbox enclosing the target EPSG:3413 window BEFORE calling rio.reproject (expensive), instead of reprojecting the whole of Greenland just to keep a window of a few tens of km.

def crop_to_extent(da, x_min, x_max, y_min, y_max) -- crops a reprojected DataArray to an exact EPSG:3413 extent.

def build_target_grid(x_center, y_center, half_window, res) -- builds the common target grid (EPSG:3413, regular spacing) every satellite product is reprojected onto, so the five products end up aligned pixel-for-pixel (avoids the moire pattern produced by chaining two independent nearest-neighbour resamplings -- see the V1 fix in the module docstring above).
"""

# ── Per-pixel metric fits ────────────────────────────────────────────────
"""
def fit_sinus_amplitude_fast(da) -- fits a single sinusoid per pixel and returns the seasonal amplitude.

def fit_mae_residuals_fast(da, poly_degree, n_harmonics) -- fits a polynomial trend + harmonic seasonality per pixel (like decompose_trend_seasonal) and returns the mean absolute error of the residuals.
"""

# ── Metric computation on the target grid ────────────────────────────────
"""
def compute_metric_on_target_grid(satellite_names, x_center, y_center, half_window, target_res, fit_fn, ...) -- generic per-satellite pipeline: restrict time range, crop native grid, reproject onto the common target grid, apply fit_fn per pixel. Shared by the amplitude and MAE computations below via an interchangeable fit_fn.

def compute_amplitude_on_target_grid(...) -- compute_metric_on_target_grid with fit_fn=fit_sinus_amplitude_fast, for one satellite.

def compute_mae_on_target_grid(...) -- compute_metric_on_target_grid with fit_fn=fit_mae_residuals_fast, for one satellite.

def compute_mean_metric_for_locality(satellite_names, x_center, y_center, half_window, target_res, fit_fn, ...) -- averages compute_metric_on_target_grid's result across every requested satellite, for one locality/window.

def compute_mean_amplitude_for_locality(...) -- compute_mean_metric_for_locality with the amplitude fit.

def compute_mean_mae_for_locality(...) -- compute_mean_metric_for_locality with the MAE fit.
"""

# ── Small formatting / validation helpers ────────────────────────────────
"""
def order_localities(locality_names, order) -- sorts locality names according to LOCALITY_ORDER (falling back to alphabetical for anything not listed).

def build_group_label(satellite_names, max_listed) -- builds a short, human-readable label listing the requested satellite products (elided beyond max_listed).

def sanitize_for_filename(label) -- strips/replaces characters unsafe for a filename.

def safe_output_path(out_dir, fig_name, max_path_len) -- truncates + hashes the filename if the full output path would exceed max_path_len characters (Windows MAX_PATH).

def wrap_title(title, fig_width_inches, chars_per_inch, min_chars) -- wraps a figure title to the actual figure width, so a long group_label doesn't overflow a narrow (n_cols=1) figure.

def validate_satellite_names(satellite_names) -- raises if any requested satellite name isn't a known product.
"""

# ── Figures: grid of subplots by locality ────────────────────────────────
"""
def plot_amplitude_map_by_locality_grid(...) -- grid of subplots, one per locality (grouped by group_stations_by_locality), each a metric map with location inset, colorbar, and ice terminus overlay.

def plot_mae_map_by_locality_grid(...) -- wrapper around plot_amplitude_map_by_locality_grid, passing the MAE fit/labels instead of the amplitude ones.

def render_amplitude_panel(ax, mean_metric, x_center, y_center, half_window, land_paths, coast_paths, gdf_terminus, cmap, vmin, vmax, var_alpha, ...) -- renders one metric panel (pcolormesh + land/ocean mask + ice terminus + optional station markers) onto a given axis; the shared drawing routine behind every figure function in this module.
"""

# ── Figures: Greenland-wide ───────────────────────────────────────────────
"""
def plot_amplitude_map_greenland(satellite_names, ..., compute_mean_fn, metric_label, metric_prefix) -- metric map over the whole of Greenland (no per-locality split), on a coarser default grid to stay tractable.

def plot_mae_map_greenland(satellite_names, poly_degree, n_harmonics, **kwargs) -- shortcut for plot_amplitude_map_greenland using the MAE-of-residuals fit instead of the seasonal amplitude.
"""

# ── Figures: arbitrary point ──────────────────────────────────────────────
"""
def plot_amplitude_map_at_location(satellite_names, lat, lon, radius_km, ..., compute_mean_fn, metric_label, metric_prefix) -- metric map centred on a freely chosen (lat, lon) point (WGS84) rather than a predefined stations_dict locality.

def plot_mae_map_at_location(satellite_names, lat, lon, poly_degree, n_harmonics, **kwargs) -- shortcut for plot_amplitude_map_at_location using the MAE-of-residuals fit instead of the seasonal amplitude.
"""

# ── Paths and constants ───────────────────────────────────────────────────
"""
ICE_TERMINUS_GPKG_PATH, ICE_TERMINUS_COLORS, ICE_TERMINUS_LINEWIDTH -- PROMICE ice terminus line file (imported via paths.py's ICE_MASK_DATA_DIR) and its display style.

LAT_KEY, LON_KEY, EPSG3413_CRS, TRANSFORMER_4326_TO_3413, TRANSFORMER_3413_TO_4326, HALF_WINDOW_M, TARGET_RES_M, RESAMPLING_METHOD, CROP_BUFFER_FACTOR, VAR_ALPHA, COPERNICUS_SAT_NAME, COPERNICUS_MIN_DATE, LOCALITY_ORDER, GREENLAND_CENTER_XY, GREENLAND_HALF_WINDOW_M, GREENLAND_TARGET_RES_M -- shared column-name keys, coordinate transformers/CRS, and the adjustable parameters controlling window size, grid resolution, resampling method and the Copernicus pre-2015 time restriction (see _restrict_time_range above).
"""


LAT_KEY = "lat"
LON_KEY = "lon"

EPSG3413_CRS = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)
TRANSFORMER_4326_TO_3413 = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
TRANSFORMER_3413_TO_4326 = Transformer.from_crs("EPSG:3413", "EPSG:4326", always_xy=True)

# --- Paramètres ajustables -------------------------------------------------
HALF_WINDOW_M = 60_000          # demi-largeur de la fenêtre autour de chaque localité (m)
TARGET_RES_M = 5_000            # résolution de la grille commune d'amplitude (m)
RESAMPLING_METHOD = Resampling.bilinear  # bilinear (upsampling) ou average (downsampling)
CROP_BUFFER_FACTOR = 1.4        # marge de recadrage précoce avant le fit (perf)

VAR_ALPHA = 1.0                 # opacité de l'overlay d'amplitude

# Filtrage temporel spécifique à Copernicus : produit jugé non fiable avant
# 2015 (cf. anomalies SEC identifiées précédemment) -- exclu des fits en
# amont, pour rester valable sur tous les plots (amplitude, MAE, ...) sans
# dupliquer le filtre dans chaque fit_fn.
COPERNICUS_SAT_NAME = "Copernicus_Climate_Data_Store"
COPERNICUS_MIN_DATE = "2015-01-01"

# Ordre d'affichage voulu pour les localités (au lieu de l'ordre alphabétique)
LOCALITY_ORDER = ["KPC", "TAS", "QAS", "NUK", "KAN", "SWC", "JAR", "UPE", "THU"]

# Trait de terminus PROMICE (déjà en EPSG:3413) -- à adapter au chemin local
ICE_TERMINUS_GPKG_PATH = ICE_MASK_DATA_DIR / "01-PROMICE-2022-IceMask-line-v3.gpkg"
ICE_TERMINUS_COLORS = {"marine": "#1f5fa8", "land": "#a0522d"}
ICE_TERMINUS_LINEWIDTH = 1.6

# Emprise Groenland entier en EPSG:3413 (dérivée des bounds réelles du
# fichier terminus PROMICE : x [-620000, 837000], y [-3308000, -790000]).
# Fenêtre carrée centrée, avec marge, pour couvrir toute la calotte.
GREENLAND_CENTER_XY = (108_500, -2_049_000)
GREENLAND_HALF_WINDOW_M = 1_300_000
GREENLAND_TARGET_RES_M = 5_000  # grille plus grossière que par localité (perf)


# --- Utilitaires stations ---------------------------------------------------

def group_stations_by_locality(stations_dict):
    localities = {}
    for station_name in stations_dict:
        locality = station_name.split("_")[0]
        localities.setdefault(locality, []).append(station_name)
    return localities


def get_station_coordinates(stations_dict, lat_key=LAT_KEY, lon_key=LON_KEY):
    coords = {}
    for station_name, info in stations_dict.items():
        lat, lon = None, None

        if lat_key in info and lon_key in info:
            lat, lon = info[lat_key], info[lon_key]
        elif "hourly_data" in info:
            try:
                df = pd.read_csv(info["hourly_data"], usecols=[lat_key, lon_key])
                lat = df[lat_key].dropna().mean()
                lon = df[lon_key].dropna().mean()
            except Exception:
                pass

        if lat is None or lon is None:
            continue

        coords[station_name] = (lon, lat)

    return coords


def locality_centroid_xy(station_names, coords):
    """Centre (x, y) EPSG:3413 d'une localité + coordonnées projetées de ses stations."""
    lons, lats = [], []
    for s in station_names:
        if s in coords:
            lon, lat = coords[s]
            lons.append(lon)
            lats.append(lat)

    x_st, y_st = TRANSFORMER_4326_TO_3413.transform(np.array(lons), np.array(lats))
    return x_st.mean(), y_st.mean(), x_st, y_st


def add_location_inset(
    fig, ax, x_min, x_max, y_min, y_max,
    inset_width_frac=0.22, inset_height_frac=0.22, margin_frac=0.03,
):
    """Mini-carte Groenland dans le coin inférieur gauche de `ax`, avec l'emprise en rouge."""
    bbox = ax.get_position()

    inset_w = bbox.width * inset_width_frac
    inset_h = bbox.height * inset_height_frac
    margin_x = bbox.width * margin_frac
    margin_y = bbox.height * margin_frac

    left = bbox.x0 + margin_x
    bottom = bbox.y0 + margin_y

    ax_inset = fig.add_axes([left, bottom, inset_w, inset_h], projection=EPSG3413_CRS)
    ax_inset.set_extent([-60, -30, 59, 84], crs=ccrs.PlateCarree())

    ax_inset.add_feature(cfeature.LAND.with_scale("50m"), facecolor="0.85", edgecolor="none", zorder=1)
    ax_inset.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="#cfe8f3", zorder=0)
    ax_inset.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.4, zorder=2)

    rect = mpatches.Rectangle(
        (x_min, y_min), x_max - x_min, y_max - y_min,
        linewidth=1.1, edgecolor="red", facecolor="none",
        transform=EPSG3413_CRS, zorder=5,
    )
    ax_inset.add_patch(rect)

    ax_inset.set_xticks([])
    ax_inset.set_yticks([])
    for spine in ax_inset.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(0.6)

    return ax_inset


def get_land_ocean_paths(scale="50m"):
    """
    Précalcule les tracés terre/côte du Groenland projetés en EPSG:3413.
    Les axes principaux sont des Axes matplotlib "plats" (pas de GeoAxes
    cartopy, cf. plot_regional_variation), donc on ne peut pas utiliser
    ax.add_feature(...) directement : on reprojette les géométries
    Natural Earth (lat/lon) vers EPSG:3413 avec project_geometry, puis on
    les convertit en Path matplotlib avec geos_to_path pour pouvoir les
    ajouter comme patches ordinaires. Calculé une seule fois pour toute la
    figure, réutilisé sur chaque sous-graphe.
    """
    land_geoms = list(cfeature.LAND.with_scale(scale).geometries())
    coast_geoms = list(cfeature.COASTLINE.with_scale(scale).geometries())

    land_paths = []
    for geom in land_geoms:
        proj_geom = EPSG3413_CRS.project_geometry(geom, ccrs.PlateCarree())
        land_paths.extend(geos_to_path(proj_geom))

    coast_paths = []
    for geom in coast_geoms:
        proj_geom = EPSG3413_CRS.project_geometry(geom, ccrs.PlateCarree())
        coast_paths.extend(geos_to_path(proj_geom))

    return land_paths, coast_paths


def add_land_ocean_mask(
    ax, land_paths, coast_paths, x_center, y_center, res=1000.0,
    ocean_color="#d8eef7", land_color="0.88",
):
    """
    Fond terre (gris) / mer (bleu clair) + trait de côte, sous les données.

    x_center, y_center, res : les tracés terre/côte sont précalculés une
    seule fois en coordonnées EPSG:3413 absolues (mètres). Comme les axes
    principaux sont maintenant en coordonnées relatives au centre de la
    localité (km, pour permettre sharex/sharey), on translate + met à
    l'échelle chaque tracé au moment de l'ajouter à cet axe précis :
    (x - x_center) / res, (y - y_center) / res.
    """
    transform = Affine2D().translate(-x_center, -y_center).scale(1.0 / res)

    ax.set_facecolor(ocean_color)
    for path in land_paths:
        ax.add_patch(mpatches.PathPatch(path.transformed(transform), facecolor=land_color, edgecolor="none", zorder=0))
    for path in coast_paths:
        ax.add_patch(mpatches.PathPatch(path.transformed(transform), facecolor="none", edgecolor="black", linewidth=0.5, zorder=1))


# --- Trait de terminus PROMICE ----------------------------------------------

_ICE_TERMINUS_CACHE = {}


def get_ice_terminus(gpkg_path=ICE_TERMINUS_GPKG_PATH):
    """
    Charge (et met en cache) le trait de terminus PROMICE. Le fichier est
    déjà en EPSG:3413 -- reprojeté seulement si ce n'est pas le cas (ex. si
    tu changes de fichier source pour une version en lat/lon).
    """
    if gpkg_path not in _ICE_TERMINUS_CACHE:
        gdf = gpd.read_file(gpkg_path)
        if gdf.crs is None or gdf.crs.to_epsg() != 3413:
            gdf = gdf.to_crs("EPSG:3413")
        _ICE_TERMINUS_CACHE[gpkg_path] = gdf
    return _ICE_TERMINUS_CACHE[gpkg_path]


def add_ice_terminus(
    ax, gdf_terminus, x_center, y_center, half_window, res=1000.0,
    buffer_m=5_000, colors=ICE_TERMINUS_COLORS, linewidth=ICE_TERMINUS_LINEWIDTH, zorder=3,
):
    """
    Superpose le trait de terminus (PROMICE) sur la grille d'amplitude, dans
    les mêmes coordonnées relatives (km depuis le centre) que le reste du
    subplot -- permet de voir si les pixels de la grille sont en amont
    (intérieur de la calotte, pas de trait visible dans la fenêtre) ou à
    cheval sur la limite (trait traversant des cellules colorées).

    Recadrage précoce via `.cx[...]` (bbox géopandas) pour ne convertir en
    Path matplotlib que les segments réellement dans la fenêtre + une marge
    `buffer_m`, plutôt que les 1677 features du fichier complet.

    Retourne les handles de légende (un par type de terminus rencontré dans
    cette fenêtre), à agréger entre localités pour une légende commune.
    """
    x_min = x_center - half_window - buffer_m
    x_max = x_center + half_window + buffer_m
    y_min = y_center - half_window - buffer_m
    y_max = y_center + half_window + buffer_m

    cropped = gdf_terminus.cx[x_min:x_max, y_min:y_max]
    if cropped.empty:
        return {}

    transform = Affine2D().translate(-x_center, -y_center).scale(1.0 / res)

    handles = {}
    for _, row in cropped.iterrows():
        termini_type = row.get("Termini", "unknown")
        color = colors.get(termini_type, "black")

        for path in geos_to_path(row.geometry):
            ax.add_patch(
                mpatches.PathPatch(
                    path.transformed(transform), facecolor="none",
                    edgecolor=color, linewidth=linewidth, zorder=zorder,
                )
            )

        if termini_type not in handles:
            handles[termini_type] = mlines.Line2D(
                [], [], color=color, linewidth=linewidth, label=f"{termini_type} terminus"
            )

    return handles


# --- Reprojection / recadrage (grille d'amplitude) --------------------------

def _restrict_time_range(ds, sat_name, min_date=COPERNICUS_MIN_DATE):
    """Restreint la dimension temporelle pour un satellite donné (ex.
    Copernicus à partir de 2015). No-op pour tout autre satellite, ou si le
    dataset n'a pas de dimension 'time'. Appelé avant tout recadrage/fit,
    donc s'applique uniformément à l'amplitude et à la MAE, sur la grille
    par localité, le Groenland entier et les points libres."""
    if sat_name != COPERNICUS_SAT_NAME or "time" not in ds.dims:
        return ds
    return ds.sel(time=slice(min_date, None))


def _to_epsg3413(ds, var_name):
    """
    Reprojette le DataArray en EPSG:3413 (grille native du produit).
    Gère les grilles géographiques (lat/lon, ex. Copernicus) et les grilles
    déjà projetées x/y (ex. Zhang).
    """
    da = ds[var_name]

    if "time" in da.dims and da.dims[0] != "time":
        da = da.transpose("time", ...)

    lat_dim = next((d for d in ["lat", "latitude", "y"] if d in da.dims), None)
    lon_dim = next((d for d in ["lon", "longitude", "x"] if d in da.dims), None)

    if lat_dim is None or lon_dim is None:
        raise ValueError(
            f"Impossible de trouver des dimensions spatiales (lat/lon ou x/y) pour "
            f"'{var_name}'. Dimensions disponibles : {tuple(da.dims)}"
        )

    is_geographic = lat_dim in ("lat", "latitude") and lon_dim in ("lon", "longitude")

    da = da.rio.write_crs("EPSG:4326" if is_geographic else "EPSG:3413", inplace=True)
    da = da.rio.set_spatial_dims(x_dim=lon_dim, y_dim=lat_dim, inplace=True)

    if is_geographic:
        da = da.rio.reproject("EPSG:3413")

    return da


def _latlon_bbox_from_epsg3413_window(x_center, y_center, buf, pad_deg=2.0):
    """
    Bbox lat/lon (avec marge pad_deg) englobant une fenêtre carrée
    EPSG:3413 [x_center±buf, y_center±buf]. Transforme les 4 coins + les 4
    milieux de côtés (pas seulement les coins) vers WGS84 avant de prendre
    min/max, car une fenêtre carrée en projection polaire stéréographique
    ne correspond pas à un rectangle propre en lat/lon -- les coins seuls
    peuvent sous-estimer l'étendue réelle. La marge pad_deg absorbe le
    reste de la distorsion ; le recadrage final précis est de toute façon
    refait après reprojection par crop_to_extent.
    """
    xs = [
        x_center - buf, x_center - buf, x_center + buf, x_center + buf,
        x_center, x_center, x_center - buf, x_center + buf,
    ]
    ys = [
        y_center - buf, y_center + buf, y_center - buf, y_center + buf,
        y_center - buf, y_center + buf, y_center, y_center,
    ]
    lons, lats = TRANSFORMER_3413_TO_4326.transform(xs, ys)
    lon_min, lon_max = min(lons), max(lons)
    lat_min, lat_max = min(lats), max(lats)
    return lon_min - pad_deg, lon_max + pad_deg, lat_min - pad_deg, lat_max + pad_deg


def _crop_native_before_reproject(ds, var_name, x_center, y_center, buf, pad_deg=2.0):
    """
    Recadre le dataset natif à une bbox lat/lon englobant la fenêtre
    EPSG:3413 visée, AVANT toute reprojection -- évite de reprojeter
    l'intégralité du Groenland (rio.reproject, coûteux) pour ne garder
    ensuite qu'une fenêtre de quelques dizaines de km.

    Ne s'applique que si la grille native est géographique (lat/lon, ex.
    Copernicus) sur une grille régulière 1D. Pour une grille déjà projetée
    (ex. Zhang, déjà en EPSG:3413), _to_epsg3413 ne reprojette pas -- rien à
    gagner, on retourne ds tel quel. Idem pour une grille curviligne
    (lat/lon 2D dépendant de x/y) : pas de recadrage précoce simple
    possible, on laisse la reprojection complète gérer le cas (comportement
    inchangé, juste pas d'optimisation).
    """
    da = ds[var_name]
    lat_dim = next((d for d in ["lat", "latitude", "y"] if d in da.dims), None)
    lon_dim = next((d for d in ["lon", "longitude", "x"] if d in da.dims), None)

    if lat_dim is None or lon_dim is None:
        return ds
    is_geographic = lat_dim in ("lat", "latitude") and lon_dim in ("lon", "longitude")
    if not is_geographic:
        return ds

    lon_coord = ds[lon_dim]
    lat_coord = ds[lat_dim]
    if lon_coord.ndim != 1 or lat_coord.ndim != 1:
        return ds  # grille curviligne -- pas de recadrage précoce simple

    lon_min, lon_max, lat_min, lat_max = _latlon_bbox_from_epsg3413_window(x_center, y_center, buf, pad_deg)

    # gère le cas où les longitudes natives sont en 0-360 plutôt que -180/180
    lon_vals = lon_coord.values
    if lon_vals.max() > 180.0 and lon_min < 0:
        lon_min, lon_max = lon_min % 360, lon_max % 360

    lon_mask = (lon_vals >= lon_min) & (lon_vals <= lon_max)
    lat_mask = (lat_coord.values >= lat_min) & (lat_coord.values <= lat_max)

    lon_idx = np.where(lon_mask)[0]
    lat_idx = np.where(lat_mask)[0]
    if lon_idx.size == 0 or lat_idx.size == 0:
        return ds  # rien dans la fenêtre -- laisse la suite lever l'erreur normalement

    try:
        return ds.isel({
            lon_dim: slice(lon_idx.min(), lon_idx.max() + 1),
            lat_dim: slice(lat_idx.min(), lat_idx.max() + 1),
        })
    except Exception:
        return ds  # repli sécurisé : comportement inchangé si l'isel échoue


def crop_to_extent(da, x_min, x_max, y_min, y_max):
    """Découpe sécurisée d'un DataArray projeté (avec repli si clip_box échoue)."""
    try:
        return da.rio.clip_box(minx=x_min, miny=y_min, maxx=x_max, maxy=y_max)
    except Exception:
        x_slice = slice(x_min, x_max) if da.x[0] < da.x[-1] else slice(x_max, x_min)
        y_slice = slice(y_min, y_max) if da.y[0] < da.y[-1] else slice(y_max, y_min)
        return da.sel(x=x_slice, y=y_slice)


def build_target_grid(x_center, y_center, half_window, res):
    """Grille régulière carrée centrée sur (x_center, y_center)."""
    xs = np.arange(x_center - half_window, x_center + half_window, res)
    ys = np.arange(y_center - half_window, y_center + half_window, res)

    target = xr.DataArray(
        np.full((len(ys), len(xs)), np.nan), coords={"y": ys, "x": xs}, dims=["y", "x"]
    )
    target = target.rio.write_crs("EPSG:3413", inplace=True)
    target = target.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
    return target


# --- Fits pixel par pixel (vectorisés) --------------------------------------

def fit_sinus_amplitude_fast(da):
    """Fit sinusoïdal pixel par pixel, résistant aux NaNs (moindres carrés via lstsq).

    Modèle : y = A*sin(ωt) + B*cos(ωt) + offset, ω fixé au cycle annuel
    (365.25 j). Pas de retrait de tendance polynomiale ici (cohérent avec le
    comportement historique de cette figure) : cf. fit_mae_residuals_fast
    pour une version avec tendance retirée, utilisée pour la MAE des
    résidus.
    """
    time = pd.to_datetime(da["time"].values)
    t_days = (time - time[0]).days.astype(float)
    T = len(t_days)

    omega = 2 * np.pi / 365.25
    S = np.sin(omega * t_days)
    C = np.cos(omega * t_days)
    O = np.ones_like(t_days)

    X = np.vstack([S, C, O]).T

    vals = da.values
    shape_spatial = vals.shape[1:]
    y_flat = vals.reshape(T, -1)

    coeffs, _, _, _ = np.linalg.lstsq(X, np.nan_to_num(y_flat, nan=0.0), rcond=None)

    A_sin = coeffs[0, :]
    B_cos = coeffs[1, :]
    amplitude = np.sqrt(A_sin**2 + B_cos**2)

    mask_nan = np.isnan(y_flat).all(axis=0)
    amplitude[mask_nan] = np.nan

    amp_grid = amplitude.reshape(shape_spatial)

    spatial_coords = {k: da[k] for k in da.dims if k != "time"}
    spatial_dims = [d for d in da.dims if d != "time"]

    amp_da = xr.DataArray(amp_grid, coords=spatial_coords, dims=spatial_dims, name="amplitude")

    src_crs = da.rio.crs if hasattr(da, "rio") else None
    if src_crs is not None:
        amp_da = amp_da.rio.write_crs(src_crs, inplace=True)
        x_dim = "x" if "x" in amp_da.dims else da.rio.x_dim
        y_dim = "y" if "y" in amp_da.dims else da.rio.y_dim
        amp_da = amp_da.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim, inplace=True)

    return amp_da


def fit_mae_residuals_fast(da, poly_degree=2, n_harmonics=1):
    """Fit tendance polynomiale + saisonnalité harmonique, pixel par pixel
    (moindres carrés via lstsq, vectorisé), puis MAE des résidus
    (signal - tendance - saisonnalité) sur toute la période disponible.

    Même décomposition que `decompose_trend_seasonal` dans plot_trend_seasonality_residues.py
    (tendance polynomiale de degré `poly_degree` + `n_harmonics` harmoniques
    sin/cos), mais appliquée simultanément à tous les pixels d'une grille
    plutôt qu'à une seule série station/satellite.

    Résistant aux NaNs : les pixels sans aucune observation valide
    ressortent en NaN ; les NaNs ponctuels dans le temps sont ignorés dans
    le calcul de la MAE (moyenne des |résidu| sur les seules dates valides).
    """
    time = pd.to_datetime(da["time"].values)
    t_years = (time - time[0]).days.astype(float) / 365.25
    T = len(t_years)

    cols = [t_years**d for d in range(poly_degree + 1)]
    for k in range(1, n_harmonics + 1):
        cols.append(np.sin(2 * np.pi * k * t_years))
        cols.append(np.cos(2 * np.pi * k * t_years))
    X = np.column_stack(cols)  # (T, n_params)

    vals = da.values
    shape_spatial = vals.shape[1:]
    y_flat = vals.reshape(T, -1)  # (T, n_pixels), NaNs conservés

    valid = np.isfinite(y_flat)
    y_filled = np.nan_to_num(y_flat, nan=0.0)

    coeffs, *_ = np.linalg.lstsq(X, y_filled, rcond=None)  # (n_params, n_pixels)
    fitted = X @ coeffs  # (T, n_pixels)

    residual = np.where(valid, y_flat - fitted, np.nan)
    abs_residual = np.abs(residual)

    with np.errstate(invalid="ignore"):
        mae = np.nanmean(abs_residual, axis=0)  # (n_pixels,), NaN si aucune donnée

    mask_all_nan = ~valid.any(axis=0)
    mae[mask_all_nan] = np.nan

    mae_grid = mae.reshape(shape_spatial)

    spatial_coords = {k: da[k] for k in da.dims if k != "time"}
    spatial_dims = [d for d in da.dims if d != "time"]

    mae_da = xr.DataArray(mae_grid, coords=spatial_coords, dims=spatial_dims, name="mae_residuals")

    src_crs = da.rio.crs if hasattr(da, "rio") else None
    if src_crs is not None:
        mae_da = mae_da.rio.write_crs(src_crs, inplace=True)
        x_dim = "x" if "x" in mae_da.dims else da.rio.x_dim
        y_dim = "y" if "y" in mae_da.dims else da.rio.y_dim
        mae_da = mae_da.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim, inplace=True)

    return mae_da


# --- Calcul de la métrique moyenne, par localité, sur grille commune -------

def compute_metric_on_target_grid(
    sat_name, x_center, y_center, half_window, res, resampling,
    fit_fn, fit_kwargs=None, buffer_factor=CROP_BUFFER_FACTOR,
):
    """
    Métrique spatiale (amplitude, MAE des résidus, ...) d'UN satellite,
    recadrée tôt autour de la localité (perf) puis reprojetée en une seule
    passe sur la grille cible commune.

    fit_fn : fonction (da_crop, **fit_kwargs) -> DataArray 2D (une valeur
        par pixel), ex. fit_sinus_amplitude_fast ou fit_mae_residuals_fast.

    Perf : le recadrage précoce se fait maintenant en DEUX temps --
    _crop_native_before_reproject recadre en lat/lon (large, approximatif)
    AVANT tout appel à rio.reproject (l'étape coûteuse pour les grilles
    géographiques), puis crop_to_extent fait le recadrage précis une fois
    en EPSG:3413. Avant ce changement, rio.reproject tournait sur
    l'intégralité de la grille native avant tout recadrage.

    Filtrage temporel : _restrict_time_range est appliqué juste après
    l'ouverture du dataset (avant tout recadrage spatial/reprojection/fit),
    pour couper la série Copernicus à partir de 2015 -- no-op pour tout
    autre satellite.
    """
    fit_kwargs = fit_kwargs or {}
    info = interpolation_altimetry_AWS.SATELLITE[sat_name]
    buf = half_window * buffer_factor

    print(f"    - {sat_name}: opening + cropping native grid...", flush=True)
    ds = xr.open_dataset(info["file"])
    ds = _restrict_time_range(ds, sat_name)
    ds = _crop_native_before_reproject(ds, info["var"], x_center, y_center, buf)

    print(f"    - {sat_name}: reprojecting...", flush=True)
    da_native = _to_epsg3413(ds, info["var"])

    da_crop = crop_to_extent(da_native, x_center - buf, x_center + buf, y_center - buf, y_center + buf)

    if da_crop.sizes.get("x", 0) == 0 or da_crop.sizes.get("y", 0) == 0:
        raise ValueError("no data in the window for this product")

    print(f"    - {sat_name}: fitting ({da_crop.sizes.get('x','?')}x{da_crop.sizes.get('y','?')} px)...", flush=True)
    metric_native = fit_fn(da_crop, **fit_kwargs)

    target = build_target_grid(x_center, y_center, half_window, res)
    metric_on_target = metric_native.rio.reproject_match(target, resampling=resampling)
    print(f"    - {sat_name}: done.", flush=True)

    return metric_on_target


def compute_amplitude_on_target_grid(
    sat_name, x_center, y_center, half_window, res, resampling, buffer_factor=CROP_BUFFER_FACTOR,
):
    """Amplitude saisonnière d'UN satellite sur la grille cible commune.
    Cf. compute_metric_on_target_grid pour le détail (recadrage + reprojection)."""
    return compute_metric_on_target_grid(
        sat_name, x_center, y_center, half_window, res, resampling,
        fit_fn=fit_sinus_amplitude_fast, buffer_factor=buffer_factor,
    )


def compute_mae_on_target_grid(
    sat_name, x_center, y_center, half_window, res, resampling,
    poly_degree=2, n_harmonics=1, buffer_factor=CROP_BUFFER_FACTOR,
):
    """MAE des résidus (tendance+saisonnalité retirées), sur toute la
    période disponible, d'UN satellite sur la grille cible commune."""
    return compute_metric_on_target_grid(
        sat_name, x_center, y_center, half_window, res, resampling,
        fit_fn=fit_mae_residuals_fast,
        fit_kwargs={"poly_degree": poly_degree, "n_harmonics": n_harmonics},
        buffer_factor=buffer_factor,
    )


def compute_mean_metric_for_locality(
    satellite_names, x_center, y_center,
    half_window, res, resampling, fit_fn, fit_kwargs=None,
):
    """Moyenne pixel à pixel d'une métrique (amplitude, MAE, ...) sur
    plusieurs produits satellite, pour une localité donnée."""
    metrics = []
    for name in satellite_names:
        try:
            metric = compute_metric_on_target_grid(
                name, x_center, y_center, half_window, res, resampling,
                fit_fn=fit_fn, fit_kwargs=fit_kwargs,
            )
            metrics.append(metric)
        except Exception as e:
            print(f"  [!] {name} skipped for this location: {e}")

    if not metrics:
        return None

    stack = xr.concat(metrics, dim="product")
    mean_metric = stack.mean("product", skipna=True)
    mean_metric.name = "mean_metric"
    return mean_metric


def compute_mean_amplitude_for_locality(
    satellite_names, x_center, y_center,
    half_window=HALF_WINDOW_M, res=TARGET_RES_M, resampling=RESAMPLING_METHOD,
):
    mean_amp = compute_mean_metric_for_locality(
        satellite_names, x_center, y_center, half_window, res, resampling,
        fit_fn=fit_sinus_amplitude_fast,
    )
    if mean_amp is not None:
        mean_amp.name = "mean_amplitude"
    return mean_amp


def compute_mean_mae_for_locality(
    satellite_names, x_center, y_center,
    half_window=HALF_WINDOW_M, res=TARGET_RES_M, resampling=RESAMPLING_METHOD,
    poly_degree=2, n_harmonics=1,
):
    """Comme compute_mean_amplitude_for_locality, mais pour la MAE des
    résidus (tendance+saisonnalité retirées), moyennée pixel à pixel sur
    plusieurs produits satellite."""
    mean_mae = compute_mean_metric_for_locality(
        satellite_names, x_center, y_center, half_window, res, resampling,
        fit_fn=fit_mae_residuals_fast,
        fit_kwargs={"poly_degree": poly_degree, "n_harmonics": n_harmonics},
    )
    if mean_mae is not None:
        mean_mae.name = "mean_mae_residuals"
    return mean_mae


def order_localities(locality_names, order=LOCALITY_ORDER):
    """
    Trie les localités selon `order` (ex. LOCALITY_ORDER). Toute localité
    présente dans les données mais absente de `order` est ajoutée à la fin,
    par ordre alphabétique, plutôt que d'être silencieusement supprimée.
    """
    present = set(locality_names)
    ordered = [loc for loc in order if loc in present]
    remaining = sorted(present - set(ordered))
    return ordered + remaining


def build_group_label(satellite_names, max_listed=3):
    """
    Construit un libellé lisible pour un groupe de satellites, utilisé dans
    le titre, la colorbar et le nom de fichier — pour que deux appels avec
    des combinaisons différentes (1 satellite, 2, 3, tous...) ne s'écrasent
    pas et restent identifiables sans ambiguïté.
    """
    if len(satellite_names) == 1:
        return satellite_names[0]
    if len(satellite_names) <= max_listed:
        return " + ".join(satellite_names)
    return f"{len(satellite_names)} products ({', '.join(satellite_names[:max_listed])}, ...)"


def sanitize_for_filename(label):
    """Nettoie un libellé pour l'utiliser dans un nom de fichier (pas d'espaces/plus/virgules)."""
    return (
        label.replace(" + ", "-").replace(", ", "-").replace(" ", "_")
        .replace("(", "").replace(")", "").replace(",", "")
    )


def safe_output_path(out_dir, fig_name, max_path_len=240):
    """
    Évite le FileNotFoundError Windows quand le chemin complet dépasse
    MAX_PATH (~260 caractères) -- typiquement quand build_group_label liste
    plusieurs noms de satellites complets dans un fig_name auto-généré, sur
    un figures_dir déjà profond (Documents/ENM - Toulouse/.../figures/...).

    Si le chemin complet (out_dir / fig_name) tient dans max_path_len, il
    est retourné tel quel. Sinon, le nom de fichier est tronqué et un court
    hash (8 caractères, dérivé du nom complet) est ajouté pour éviter les
    collisions entre deux combinaisons de satellites différentes qui
    tronqueraient au même préfixe.
    """
    out_path = out_dir / fig_name
    if len(str(out_path)) <= max_path_len:
        return out_path

    import hashlib
    stem = Path(fig_name).stem
    suffix = Path(fig_name).suffix
    dir_len = len(str(out_dir)) + 1  # +1 pour le séparateur
    h = hashlib.md5(fig_name.encode()).hexdigest()[:8]
    budget = max_path_len - dir_len - len(suffix) - len(h) - 1  # -1 pour le "_" avant le hash
    budget = max(budget, 10)
    short_stem = stem[:budget]
    short_name = f"{short_stem}_{h}{suffix}"
    print(f"  [i] Nom de fichier tronqué (chemin trop long) : {fig_name} -> {short_name}")
    return out_dir / short_name


def wrap_title(title, fig_width_inches, chars_per_inch=6.0, min_chars=20):
    """
    Retour à la ligne automatique d'un titre (potentiellement multi-lignes
    via '\\n') selon la largeur réelle de la figure -- évite qu'un
    group_label très long (plusieurs satellites listés) ne déborde du cadre
    de la figure, ce qui forçait bbox_inches='tight' à agrandir le PNG
    final bien au-delà de la carte elle-même (carte minuscule au milieu
    d'un immense espace vide).

    Retourne (titre_wrappé, nombre_de_lignes) -- le nombre de lignes sert à
    dimensionner le bandeau d'en-tête (cf. plot_amplitude_map_by_locality_grid).
    """
    max_chars = max(int(fig_width_inches * chars_per_inch), min_chars)
    wrapped_lines = []
    for line in title.split("\n"):
        wrapped_lines.extend(textwrap.wrap(line, width=max_chars) or [""])
    return "\n".join(wrapped_lines), len(wrapped_lines)


def validate_satellite_names(satellite_names):
    """Vérifie que chaque nom demandé existe bien dans interpolation_altimetry_AWS.SATELLITE."""
    available = set(interpolation_altimetry_AWS.SATELLITE.keys())
    unknown = [name for name in satellite_names if name not in available]
    if unknown:
        raise ValueError(
            f"Unknown satellite product(s): {unknown}. "
            f"Available: {sorted(available)}"
        )


# --- Figure ------------------------------------------------------------------

def plot_amplitude_map_by_locality_grid(
    stations_dict,
    satellite_names,
    n_cols=3,
    cmap="viridis",
    half_window=HALF_WINDOW_M,
    target_res=TARGET_RES_M,
    var_alpha=VAR_ALPHA,
    vmin=0.0,
    vmax=3.0,
    cbar_label=None,
    group_label=None,
    title=None,
    show_ice_terminus=True,
    ice_terminus_path=ICE_TERMINUS_GPKG_PATH,
    figures_dir=None,
    fig_name=None,
    save=True,
    compute_mean_fn=compute_mean_amplitude_for_locality,
    metric_label="Seasonal mean amplitude",
    metric_prefix="Amplitude",
):
    """
    satellite_names : liste d'1 à N clés de interpolation_altimetry_AWS.SATELLITE.
        1 élément -> pas de moyenne (le produit est affiché tel quel).
        2+ éléments -> moyenne pixel à pixel des produits listés.
        Permet de comparer facilement un seul produit, un sous-groupe
        (ex. 2 des 5 satellites), ou l'ensemble, sans dupliquer le code.
    group_label : libellé du groupe pour la colorbar/le nom de fichier par
        défaut (si title/fig_name ne sont pas fournis). Si None, généré
        automatiquement à partir de satellite_names (ex. "Andersen et al., 2025",
        "Andersen + Nilsson", "5 produits (...)").
    title : texte du suptitle. Si None, généré automatiquement à partir de
        group_label/half_window/target_res (comportement précédent).
    show_ice_terminus : superpose le trait de terminus PROMICE (marine/land)
        sur chaque case, pour voir si les pixels sont en amont de la
        délimitation ou à cheval dessus.
    fig_name : nom du fichier de sortie. Si None, généré automatiquement à
        partir de group_label, pour éviter d'écraser la figure d'une
        combinaison précédente.
    compute_mean_fn : fonction (satellite_names, x_center, y_center,
        half_window, target_res) -> DataArray moyenne (ou None), utilisée
        pour calculer la métrique tracée. Par défaut, amplitude saisonnière
        moyenne (compute_mean_amplitude_for_locality). Passer
        compute_mean_mae_for_locality (éventuellement figé sur
        poly_degree/n_harmonics via functools.partial) pour tracer la MAE
        des résidus à la place — cf. plot_mae_map_by_locality_grid pour un
        raccourci prêt à l'emploi.
    metric_label, metric_prefix : utilisés pour générer les valeurs par
        défaut de cbar_label/title (metric_label) et fig_name (metric_prefix)
        quand la métrique n'est pas l'amplitude.
    vmin, vmax : bornes FIXES de la colorbar (m), identiques sur toutes les
        cases de la grille et d'un appel à l'autre -- comparaison directe
        entre localités/exécutions. Défaut 0-3m (amplitude). Les wrappers
        plot_mae_map_* passent automatiquement 0-4m.
    """
    validate_satellite_names(satellite_names)
    if group_label is None:
        group_label = build_group_label(satellite_names)
    if cbar_label is None:
        cbar_label = f"{metric_label} (m)"  # group_label déjà dans le titre ; évite un label colorbar trop long
    if fig_name is None:
        fig_name = f"{metric_prefix}_map_by_locality_grid_EPSG3413_{sanitize_for_filename(group_label)}.png"
    if title is None:
        title = (
            f"{metric_label} by location, grid points\n"
            f"({group_label}, EPSG:3413, "
            f"window {2 * half_window / 1000:.0f} km, grid resolution {target_res / 1000:.1f} km)"
        )

    localities = group_stations_by_locality(stations_dict)
    coords = get_station_coordinates(stations_dict)

    locality_names = order_localities(list(localities.keys()))
    n_localities = len(locality_names)
    n_rows = int(np.ceil(n_localities / n_cols))

    # --- Dimensionnement de la figure --------------------------------------
    # Pour une grille étroite (peu de colonnes), une largeur de 5.5" par
    # colonne ne suffit pas à contenir un titre long (plusieurs satellites
    # listés) : le texte débordait du cadre de la figure, et
    # bbox_inches='tight' agrandissait alors le PNG final bien au-delà de la
    # carte elle-même (carte minuscule au milieu d'un grand espace vide).
    # On agrandit donc CHAQUE sous-graphique (pas seulement la largeur totale)
    # quand n_cols est petit, pour que l'espace gagné serve la carte plutôt
    # que du vide -- sans effet sur les grilles déjà larges (3+ colonnes).
    BASE_SUBPLOT_INCHES = 5.5
    MIN_FIG_WIDTH_INCHES = 9.0
    subplot_size_inches = max(BASE_SUBPLOT_INCHES, MIN_FIG_WIDTH_INCHES / n_cols)
    fig_width_inches = subplot_size_inches * n_cols
    fig_height_inches = subplot_size_inches * n_rows

    title, n_title_lines = wrap_title(title, fig_width_inches)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(fig_width_inches, fig_height_inches),
        squeeze=False,
        sharex=True, sharey=True,
        constrained_layout=True,
    )
    axes = axes.ravel()

    # --- Bandeau réservé en haut (titre + légende terminus) --------------
    # Dimensionné en pouces à partir du nombre RÉEL de lignes du titre après
    # retour à la ligne (n_title_lines), donc s'adapte aussi bien à un titre
    # court (peu de satellites) qu'à un titre qui a dû être wrappé sur
    # plusieurs lignes. Converti en fraction de la hauteur de figure -- cette
    # fraction diminue automatiquement pour les grandes grilles (9 cartes),
    # où le bandeau ne représente qu'une petite portion d'une figure haute.
    LINE_HEIGHT_INCHES = 0.26
    TITLE_TOP_PAD_INCHES = 0.10
    LEGEND_GAP_INCHES = 0.05
    LEGEND_HEIGHT_INCHES = 0.24
    BOTTOM_PAD_INCHES = 0.08

    title_height_inches = n_title_lines * LINE_HEIGHT_INCHES
    header_inches = (
        TITLE_TOP_PAD_INCHES + title_height_inches + LEGEND_GAP_INCHES
        + LEGEND_HEIGHT_INCHES + BOTTOM_PAD_INCHES
    )
    top_frac = 1.0 - header_inches / fig_height_inches
    top_frac = min(max(top_frac, 0.5), 0.93)
    fig.get_layout_engine().set(rect=(0, 0, 1, top_frac))

    land_paths, coast_paths = get_land_ocean_paths()
    half_window_km = half_window / 1000.0

    gdf_terminus = get_ice_terminus(ice_terminus_path) if show_ice_terminus else None
    terminus_handles = {}

    # --- Passe 1 : calcul de toutes les grilles -----------------------------
    # (vmin/vmax sont désormais fixes, passés en paramètre -- comparables
    # entre localités et entre exécutions -- donc plus besoin du percentile
    # 98 global calculé ici auparavant.)
    locality_data = {}

    for locality in locality_names:
        station_names = localities[locality]
        if not any(s in coords for s in station_names):
            continue

        x_center, y_center, x_st, y_st = locality_centroid_xy(station_names, coords)
        print(f"Computing {metric_label.lower()} grid for {locality}...")
        mean_metric = compute_mean_fn(
            satellite_names, x_center, y_center, half_window, target_res
        )
        if mean_metric is None:
            print(f"  [!] {locality} skipped (no data)")
            continue

        locality_data[locality] = {
            "mean_metric": mean_metric, "x_center": x_center, "y_center": y_center,
            "x_st": x_st, "y_st": y_st, "station_names": station_names,
        }

    # --- Passe 2 : tracé, échelle commune, une seule colorbar pour la figure ---
    plotted_axes = []
    im_ref = None
    visible_axes = []

    for idx, locality in enumerate(locality_names):
        ax = axes[idx]
        data = locality_data.get(locality)
        if data is None:
            ax.set_visible(False)
            continue

        visible_axes.append(ax)
        mean_metric = data["mean_metric"]
        x_center, y_center = data["x_center"], data["y_center"]
        x_st, y_st = data["x_st"], data["y_st"]
        station_names = data["station_names"]
        x_min, x_max = x_center - half_window, x_center + half_window
        y_min, y_max = y_center - half_window, y_center + half_window

        # --- fond terre / mer (translaté + mis à l'échelle en km relatifs) ---
        add_land_ocean_mask(ax, land_paths, coast_paths, x_center, y_center, res=1000.0)

        # --- cellules de grille jointives, en km depuis le centre de la localité ---
        vals = mean_metric.values
        x_rel = (mean_metric["x"].values - x_center) / 1000.0
        y_rel = (mean_metric["y"].values - y_center) / 1000.0

        if np.isfinite(vals).any():
            im = ax.pcolormesh(
                x_rel, y_rel, vals,
                cmap=cmap, vmin=vmin, vmax=vmax,
                shading="nearest",  # une cellule par point de grille, bords jointifs
                alpha=var_alpha, zorder=2,
            )
            if im_ref is None:
                im_ref = im

        # --- trait de terminus PROMICE, par-dessus les cellules ---
        if gdf_terminus is not None:
            new_handles = add_ice_terminus(ax, gdf_terminus, x_center, y_center, half_window, res=1000.0)
            terminus_handles.update(new_handles)

        x_st_rel = (x_st - x_center) / 1000.0
        y_st_rel = (y_st - y_center) / 1000.0

        ax.scatter(x_st_rel, y_st_rel, marker="+", s=200, linewidths=2.2, color="black", zorder=5)
        for xs, ys, label in zip(x_st_rel, y_st_rel, station_names):
            ax.annotate(
                label, xy=(xs, ys), xytext=(xs + 2, ys + 2),
                fontsize=9, fontweight="bold", color="black", zorder=6,
            )

        ax.set_xlim(-half_window_km, half_window_km)
        ax.set_ylim(-half_window_km, half_window_km)
        ax.set_aspect("equal")
        ax.set_xlabel("x (km from the center)")
        ax.set_ylabel("y (km from the center)")
        ax.set_title(locality, fontsize=13, fontweight="bold")

        plotted_axes.append((ax, x_min, x_max, y_min, y_max))

    if im_ref is not None and visible_axes:
        cbar = fig.colorbar(im_ref, ax=visible_axes, shrink=0.8, pad=0.02)
        cbar.set_label(cbar_label, fontsize=11)

    if terminus_handles:
        legend_y = 1.0 - (TITLE_TOP_PAD_INCHES + title_height_inches + LEGEND_GAP_INCHES) / fig_height_inches
        fig.legend(
            handles=list(terminus_handles.values()),
            loc="upper center", bbox_to_anchor=(0.5, legend_y),
            ncol=len(terminus_handles), frameon=False, fontsize=11,
        )

    for idx in range(n_localities, len(axes)):
        axes[idx].set_visible(False)

    title_y = 1.0 - TITLE_TOP_PAD_INCHES / fig_height_inches
    fig.suptitle(title, fontsize=15, fontweight="bold", y=title_y, va="top")

    # les insets sont ajoutés APRES constrained_layout / draw, une fois les
    # positions finales des axes connues (même logique que subplot_anomaly_uncertainty)
    fig.canvas.draw()
    for ax, x_min, x_max, y_min, y_max in plotted_axes:
        add_location_inset(fig, ax, x_min, x_max, y_min, y_max)

    if save:
        out_dir = figures_dir or FIGURES_DIR_ROOT
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = safe_output_path(out_dir, fig_name)
        fig.savefig(out_path, dpi=180, bbox_inches="tight")
        print(f"Figure saved: {out_path}")

    plt.show()
    return fig, axes


def plot_mae_map_by_locality_grid(
    stations_dict,
    satellite_names,
    poly_degree=2,
    n_harmonics=1,
    **kwargs,
):
    """Raccourci pour plot_amplitude_map_by_locality_grid, mais pour la MAE
    des résidus (tendance polynomiale + saisonnalité retirées) plutôt que
    pour l'amplitude saisonnière. Même habillage visuel, mêmes paramètres
    additionnels acceptés via **kwargs (n_cols, cmap, half_window, ...).

    poly_degree, n_harmonics : passés à decompose_trend_seasonal-équivalent
    (fit_mae_residuals_fast) pour le retrait de tendance/saisonnalité avant
    calcul de la MAE, pixel par pixel.
    """
    def _compute_mean_fn(sat_names, x_center, y_center, half_window, target_res):
        return compute_mean_mae_for_locality(
            sat_names, x_center, y_center, half_window, target_res,
            poly_degree=poly_degree, n_harmonics=n_harmonics,
        )

    kwargs.setdefault("metric_label", "Residual MAE (trend + seasonal removed)")
    kwargs.setdefault("metric_prefix", "MAE_residuals")
    kwargs.setdefault("vmin", 0.0)
    kwargs.setdefault("vmax", 4.0)

    return plot_amplitude_map_by_locality_grid(
        stations_dict, satellite_names,
        compute_mean_fn=_compute_mean_fn,
        **kwargs,
    )


# --- Panneau unique réutilisable (carte Groenland entier / point libre) ----

def render_amplitude_panel(
    ax, mean_amp, x_center, y_center, half_window,
    land_paths, coast_paths, gdf_terminus,
    cmap, vmin, vmax, var_alpha,
    marker_xy_list=None, marker_labels=None,
    ice_terminus_buffer_m=5_000,
):
    """
    Dessine un panneau de métrique unique (amplitude, MAE, ...) : fond
    terre/mer, cellules de grille jointives, trait de terminus, marqueur(s)
    optionnels. Factorisé pour être réutilisé par plot_amplitude_map_greenland
    et plot_amplitude_map_at_location (la grille par localité garde sa propre
    boucle, plus proche de subplot_anomaly_uncertainty).

    Retourne (im, terminus_handles) pour construire colorbar/légende côté appelant.
    """
    add_land_ocean_mask(ax, land_paths, coast_paths, x_center, y_center, res=1000.0)

    vals = mean_amp.values
    x_rel = (mean_amp["x"].values - x_center) / 1000.0
    y_rel = (mean_amp["y"].values - y_center) / 1000.0

    im = None
    if np.isfinite(vals).any():
        im = ax.pcolormesh(
            x_rel, y_rel, vals,
            cmap=cmap, vmin=vmin, vmax=vmax,
            shading="nearest", alpha=var_alpha, zorder=2,
        )

    terminus_handles = {}
    if gdf_terminus is not None:
        terminus_handles = add_ice_terminus(
            ax, gdf_terminus, x_center, y_center, half_window,
            res=1000.0, buffer_m=ice_terminus_buffer_m,
        )

    if marker_xy_list:
        xs_rel = [(x - x_center) / 1000.0 for x, y in marker_xy_list]
        ys_rel = [(y - y_center) / 1000.0 for x, y in marker_xy_list]
        ax.scatter(xs_rel, ys_rel, marker="+", s=200, linewidths=2.2, color="black", zorder=5)
        if marker_labels:
            for xr_, yr_, label in zip(xs_rel, ys_rel, marker_labels):
                ax.annotate(
                    label, xy=(xr_, yr_), xytext=(xr_ + 2, yr_ + 2),
                    fontsize=9, fontweight="bold", color="black", zorder=6,
                )

    half_window_km = half_window / 1000.0
    ax.set_xlim(-half_window_km, half_window_km)
    ax.set_ylim(-half_window_km, half_window_km)
    ax.set_aspect("equal")
    ax.set_xlabel("x (km from the center)")
    ax.set_ylabel("y (km from the center)")

    return im, terminus_handles


def plot_amplitude_map_greenland(
    satellite_names,
    cmap="viridis",
    center_xy=GREENLAND_CENTER_XY,
    half_window=GREENLAND_HALF_WINDOW_M,
    target_res=GREENLAND_TARGET_RES_M,
    var_alpha=VAR_ALPHA,
    vmin=0.0,
    vmax=3.0,
    cbar_label=None,
    group_label=None,
    title=None,
    show_ice_terminus=True,
    ice_terminus_path=ICE_TERMINUS_GPKG_PATH,
    stations_dict=None,
    figures_dir=None,
    fig_name=None,
    save=True,
    compute_mean_fn=compute_mean_amplitude_for_locality,
    metric_label="Seasonal mean amplitude",
    metric_prefix="Amplitude",
):
    """
    Carte de métrique sur l'ensemble du Groenland (pas de découpage par
    localité). Grille plus grossière par défaut (15 km) pour rester gérable
    en mémoire/temps de calcul sur une emprise aussi large.

    stations_dict : si fourni (ex. AWS_data.STATIONS_ablation), superpose
        toutes les stations comme repères (sans libellé, pour éviter la
        surcharge visuelle sur une carte aussi dense).
    title : texte du titre. Si None, généré automatiquement.
    fig_name : nom du fichier de sortie. Si None, généré automatiquement à
        partir de group_label.
    compute_mean_fn, metric_label, metric_prefix : cf.
        plot_amplitude_map_by_locality_grid. Utiliser
        plot_mae_map_greenland pour tracer la MAE des résidus directement.
    vmin, vmax : bornes FIXES de la colorbar (m). Défaut 0-3m (amplitude) ;
        plot_mae_map_greenland passe automatiquement 0-4m.
    """
    validate_satellite_names(satellite_names)
    if group_label is None:
        group_label = build_group_label(satellite_names)
    if cbar_label is None:
        cbar_label = f"{metric_label} (m)"  # group_label déjà dans le titre ; évite un label colorbar trop long
    if fig_name is None:
        fig_name = f"{metric_prefix}_map_Greenland_{sanitize_for_filename(group_label)}.png"

    x_center, y_center = center_xy

    print(f"Computing Greenland-wide {metric_label.lower()} grid...")
    mean_metric = compute_mean_fn(
        satellite_names, x_center, y_center, half_window, target_res
    )
    if mean_metric is None:
        raise ValueError("No data could be computed for the requested satellites/extent.")

    land_paths, coast_paths = get_land_ocean_paths()
    gdf_terminus = get_ice_terminus(ice_terminus_path) if show_ice_terminus else None

    marker_xy_list = None
    if stations_dict is not None:
        coords = get_station_coordinates(stations_dict)
        marker_xy_list = [
            TRANSFORMER_4326_TO_3413.transform(lon, lat) for lon, lat in coords.values()
        ]

    fig, ax = plt.subplots(figsize=(9, 11), constrained_layout=True)

    im, terminus_handles = render_amplitude_panel(
        ax, mean_metric, x_center, y_center, half_window,
        land_paths, coast_paths, gdf_terminus,
        cmap, vmin, vmax, var_alpha,
        marker_xy_list=marker_xy_list, marker_labels=None,
        ice_terminus_buffer_m=0,  # emprise déjà maximale, pas besoin de marge
    )

    if im is not None:
        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label(cbar_label, fontsize=11)

    if terminus_handles:
        ax.legend(handles=list(terminus_handles.values()), loc="lower left", fontsize=9, frameon=True)

    if title is None:
        title = (
            f"{metric_label} — Greenland\n"
            f"({group_label}, EPSG:3413, grid resolution {target_res / 1000:.1f} km)"
        )
    ax.set_title(title, fontsize=14, fontweight="bold")

    if save:
        out_dir = figures_dir or FIGURES_DIR_ROOT
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = safe_output_path(out_dir, fig_name)
        fig.savefig(out_path, dpi=180, bbox_inches="tight")
        print(f"Figure saved: {out_path}")

    plt.show()
    return fig, ax


def plot_mae_map_greenland(
    satellite_names,
    poly_degree=2,
    n_harmonics=1,
    **kwargs,
):
    """Raccourci pour plot_amplitude_map_greenland, mais pour la MAE des
    résidus plutôt que pour l'amplitude saisonnière."""
    def _compute_mean_fn(sat_names, x_center, y_center, half_window, target_res):
        return compute_mean_mae_for_locality(
            sat_names, x_center, y_center, half_window, target_res,
            poly_degree=poly_degree, n_harmonics=n_harmonics,
        )

    kwargs.setdefault("metric_label", "Residual MAE (trend + seasonal removed)")
    kwargs.setdefault("metric_prefix", "MAE_residuals")
    kwargs.setdefault("vmin", 0.0)
    kwargs.setdefault("vmax", 4.0)

    return plot_amplitude_map_greenland(
        satellite_names,
        compute_mean_fn=_compute_mean_fn,
        **kwargs,
    )


def plot_amplitude_map_at_location(
    satellite_names,
    lat,
    lon,
    radius_km=40,
    location_label=None,
    cmap="viridis",
    target_res=TARGET_RES_M,
    var_alpha=VAR_ALPHA,
    vmin=0.0,
    vmax=3.0,
    cbar_label=None,
    group_label=None,
    title=None,
    show_ice_terminus=True,
    ice_terminus_path=ICE_TERMINUS_GPKG_PATH,
    figures_dir=None,
    fig_name=None,
    save=True,
    compute_mean_fn=compute_mean_amplitude_for_locality,
    metric_label="Seasonal mean amplitude",
    metric_prefix="Amplitude",
):
    """
    Carte de métrique centrée sur un point choisi librement (lat/lon WGS84),
    plutôt que sur une localité prédéfinie de stations_dict. Utile pour
    explorer un point précis sans avoir besoin d'une station AWS à proximité.

    lat, lon : coordonnées WGS84 (degrés décimaux) du centre voulu.
    radius_km : demi-largeur de la fenêtre carrée autour de ce point.
    title : texte du titre. Si None, généré automatiquement.
    fig_name : nom du fichier de sortie. Si None, généré automatiquement à
        partir de location_label et group_label.
    compute_mean_fn, metric_label, metric_prefix : cf.
        plot_amplitude_map_by_locality_grid. Utiliser
        plot_mae_map_at_location pour tracer la MAE des résidus directement.
    vmin, vmax : bornes FIXES de la colorbar (m). Défaut 0-3m (amplitude) ;
        plot_mae_map_at_location passe automatiquement 0-4m.
    """
    validate_satellite_names(satellite_names)
    if group_label is None:
        group_label = build_group_label(satellite_names)
    if location_label is None:
        location_label = f"{lat:.3f}N, {lon:.3f}E"
    if cbar_label is None:
        cbar_label = f"{metric_label} (m)"  # group_label déjà dans le titre ; évite un label colorbar trop long
    if fig_name is None:
        fig_name = (
            f"{metric_prefix}_map_{sanitize_for_filename(location_label)}_"
            f"{sanitize_for_filename(group_label)}.png"
        )

    x_center, y_center = TRANSFORMER_4326_TO_3413.transform(lon, lat)
    half_window = radius_km * 1000.0

    print(f"Computing {metric_label.lower()} grid at ({lat:.4f}, {lon:.4f})...")
    mean_metric = compute_mean_fn(
        satellite_names, x_center, y_center, half_window, target_res
    )
    if mean_metric is None:
        raise ValueError("No data available for this location/extent.")

    land_paths, coast_paths = get_land_ocean_paths()
    gdf_terminus = get_ice_terminus(ice_terminus_path) if show_ice_terminus else None

    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)

    im, terminus_handles = render_amplitude_panel(
        ax, mean_metric, x_center, y_center, half_window,
        land_paths, coast_paths, gdf_terminus,
        cmap, vmin, vmax, var_alpha,
        marker_xy_list=[(x_center, y_center)], marker_labels=[location_label],
    )

    if im is not None:
        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label(cbar_label, fontsize=11)

    if terminus_handles:
        ax.legend(handles=list(terminus_handles.values()), loc="lower right", fontsize=9, frameon=True)

    if title is None:
        title = (
            f"{metric_label} — {location_label}\n"
            f"({group_label}, EPSG:3413, window {2 * radius_km:.0f} km, "
            f"grid resolution {target_res / 1000:.1f} km)"
        )
    ax.set_title(title, fontsize=13, fontweight="bold")

    fig.canvas.draw()
    add_location_inset(
        fig, ax,
        x_center - half_window, x_center + half_window,
        y_center - half_window, y_center + half_window,
    )

    if save:
        out_dir = figures_dir or FIGURES_DIR_ROOT
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = safe_output_path(out_dir, fig_name)
        fig.savefig(out_path, dpi=180, bbox_inches="tight")
        print(f"Figure saved: {out_path}")

    plt.show()
    return fig, ax


def plot_mae_map_at_location(
    satellite_names,
    lat,
    lon,
    poly_degree=2,
    n_harmonics=1,
    **kwargs,
):
    """Raccourci pour plot_amplitude_map_at_location, mais pour la MAE des
    résidus plutôt que pour l'amplitude saisonnière."""
    def _compute_mean_fn(sat_names, x_center, y_center, half_window, target_res):
        return compute_mean_mae_for_locality(
            sat_names, x_center, y_center, half_window, target_res,
            poly_degree=poly_degree, n_harmonics=n_harmonics,
        )

    kwargs.setdefault("metric_label", "Residual MAE (trend + seasonal removed)")
    kwargs.setdefault("metric_prefix", "MAE_residuals")
    kwargs.setdefault("vmin", 0.0)
    kwargs.setdefault("vmax", 4.0)

    return plot_amplitude_map_at_location(
        satellite_names, lat, lon,
        compute_mean_fn=_compute_mean_fn,
        **kwargs,
    )


if __name__ == "__main__":
    # plot_amplitude_map_by_locality_grid(
    #     stations_dict=AWS_data.STATION_THU_L,
    #     satellite_names=['Nilsson and Gardner, 2026'], # ,, 'Khan et al., 2025','Copernicus_Climate_Data_Store',  'Andersen et al., 2025', 'Zhang et al., 2022'
    #     fig_name = 'Amplitude Greenland, Nilsson and Garnder, 2026',
    #     n_cols=1,
    # )

    # Carte Groenland entier — amplitude saisonnière (décommenter pour tester
    # sur l'ensemble du Groenland -- calcul long, cf. print() de progression)
    plot_amplitude_map_greenland(
        satellite_names=['Khan et al., 2025'], # 'Khan et al., 2025', 'Nilsson and Gardner, 2026'
        stations_dict=AWS_data.STATIONS_ablation,  # optionnel, mettre None pour l'omettre
        fig_name="Amplitude Greenland, Khan et al., 2025",
        vmin=None,
        vmax=None
    )

    # Carte Groenland entier — MAE des résidus (même habillage, même appel)
    plot_mae_map_greenland(
        satellite_names=['Khan et al., 2025'], #, 'Khan et al., 2025', 'Nilsson and Gardner, 2026'
        stations_dict=AWS_data.STATIONS_ablation,
        fig_name="MAE Greenland, Khan et al., 2025",
        poly_degree=2,
        n_harmonics=1,
        vmin=None,
        vmax=None
    )

    # Grille par localité — MAE des résidus (test actif : une seule station)
    # plot_mae_map_by_locality_grid(
    #     stations_dict=AWS_data.STATION_THU_L,
    #     satellite_names=['Khan et al., 2025', 'Nilsson and Gardner, 2026'], # ,,'Copernicus_Climate_Data_Store',  'Andersen et al., 2025', 'Zhang et al., 2022'
    #     fig_name = 'MAE THU_L',
    #     n_cols=1,
    # )

    # Carte centrée sur un point choisi librement (lat/lon WGS84)
    # plot_amplitude_map_at_location(
    #     satellite_names=list(interpolation_altimetry_AWS.SATELLITE.keys()),
    #     lat=69.5,
    #     lon=-49.9,
    #     radius_km=40,
    #     location_label="Custom point",  # optionnel
    # )