from pathlib import Path
import xarray as xr
import AWS_data
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import rioxarray
import interpolation_altimetry_AWS
from matplotlib.colors import LightSource
from pyproj import Transformer
from pystac_client import Client
from rasterio.io import MemoryFile
from rasterio.merge import merge

from paths import FIGURES_DIR as _BASE_FIGURES_DIR

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── ArcticDEM mosaic fetching ────────────────────────────────────────────
"""
def fetch_arcticdem(bounds_3413, resolution, version, target_resolution_m) -- fetches and merges the ArcticDEM mosaic tiles covering the given EPSG:3413 bounding box (static multi-epoch mosaic, not individually-dated strips -- see transect.py for the strip-based, time-resolved alternative).

def _auto_target_resolution(radius_km, native_resolution, max_pixels_per_side) -- picks a target resolution (m) so the fetched DEM never exceeds max_pixels_per_side pixels per side, whatever the requested radius.
"""

# ── Region geometry helpers ───────────────────────────────────────────────
"""
def _compute_center_and_bounds(aws_path, dataset, dataset_variable, start_date, end_date, radius_km) -- computes the EPSG:3413 center and bounding box (x_min, y_min, x_max, y_max) of the region to plot around a station.

def add_location_inset(fig, ax, x_min, x_max, y_min, y_max, ...) -- adds a small Greenland-wide location inset showing where the plotted region sits.
"""

# ── Diagnostics ───────────────────────────────────────────────────────────
"""
def print_extrema_with_uncertainty(aws_name, radius, label, anomaly_da, uncertainty_da) -- prints the min/max of anomaly_da to the terminal, alongside the uncertainty at the same location if uncertainty_da is given.
"""

# ── Regional map plotting ─────────────────────────────────────────────────
"""
def plot_regional_variation(aws, dataset, dataset_variable, start_date, end_date, radius_km, aws_name, ...) -- plots the spatial variation of a satellite variable around an AWS station, overlaid on an ArcticDEM hillshade. Accepts a precomputed center_xy/dem (e.g. from a caller that fetches the DEM once per radius) to avoid re-fetching the same DEM for every dataset in a column sharing that radius.

def subplot_anomaly_uncertainty(aws_path, radius, aws_name, ...) -- grid of subplots showing plot_regional_variation's anomaly map next to its uncertainty map, for one or more datasets/variables at one station.
"""

# ── Paths ─────────────────────────────────────────────────────────────────
"""
FIGURES_DIR, STAC_URL, EPSG3413_CRS, VALID_ARCTICDEM_RESOLUTIONS -- output folder (subfolder of paths.py's shared FIGURES_DIR), PGC STAC endpoint used by fetch_arcticdem, the cartopy CRS matching EPSG:3413, and the ArcticDEM mosaic's valid native resolutions (v4.1).
"""


FIGURES_DIR = _BASE_FIGURES_DIR
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

STAC_URL = "https://stac.pgc.umn.edu/api/v1/"

# Projection équivalente à EPSG:3413 pour cartopy (Polaire Stéréo Nord)
EPSG3413_CRS = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)

# Résolutions natives disponibles pour ArcticDEM mosaics v4.1
VALID_ARCTICDEM_RESOLUTIONS = (2, 10, 32)


def fetch_arcticdem(
    bounds_3413, resolution=32, version="v4.1", target_resolution_m=None
):
    """
    Récupère et fusionne les tuiles ArcticDEM couvrant l'emprise donnée.

    Parameters
    ----------
    bounds_3413 : tuple
        (x_min, y_min, x_max, y_max) en EPSG:3413.
    resolution : int
        Résolution native à interroger sur le STAC (doit être 2, 10 ou 32 pour v4.1).
    target_resolution_m : float, optional
        Si fourni et > resolution, sous-échantillonne après récupération pour
        atteindre approximativement cette résolution finale (ex: 128 pour un
        gros rayon).
    """
    if resolution not in VALID_ARCTICDEM_RESOLUTIONS:
        raise ValueError(
            f"resolution={resolution} invalide. "
            f"Valeurs autorisées pour ArcticDEM {version} : {VALID_ARCTICDEM_RESOLUTIONS}."
        )

    x_min, y_min, x_max, y_max = bounds_3413

    transformer = Transformer.from_crs("EPSG:3413", "EPSG:4326", always_xy=True)
    lon_min, lat_min = transformer.transform(x_min, y_min)
    lon_max, lat_max = transformer.transform(x_max, y_max)
    bbox = [
        min(lon_min, lon_max),
        min(lat_min, lat_max),
        max(lon_min, lon_max),
        max(lat_min, lat_max),
    ]

    catalog = Client.open(STAC_URL)
    collection_id = f"arcticdem-mosaics-{version}-{resolution}m"
    search = catalog.search(collections=[collection_id], bbox=bbox)
    items = list(search.items())

    if not items:
        raise ValueError("Aucune tuile ArcticDEM trouvée pour cette emprise.")

    datasets = []
    for item in items:
        asset = item.assets.get("dem") or item.assets.get("data")
        if asset is None:
            continue
        datasets.append(rasterio.open(asset.href))

    if not datasets:
        raise ValueError(
            "Aucun asset DEM exploitable trouvé dans les tuiles renvoyées."
        )

    mosaic_arr, mosaic_transform = merge(datasets, bounds=(x_min, y_min, x_max, y_max))

    profile = datasets[0].profile.copy()
    profile.update(
        {
            "height": mosaic_arr.shape[1],
            "width": mosaic_arr.shape[2],
            "transform": mosaic_transform,
        }
    )

    for src in datasets:
        src.close()
    del datasets  # libère les handles rasterio dès que possible

    with MemoryFile() as memfile:
        with memfile.open(**profile) as dst:
            dst.write(mosaic_arr)
        dem_da = rioxarray.open_rasterio(memfile, masked=True).squeeze(
            "band", drop=True
        )
        dem_da = dem_da.load()

    del mosaic_arr  # libère le tableau brut, plus besoin une fois dem_da chargé

    # Sous-échantillonnage optionnel si on demande une résolution finale plus grossière
    if target_resolution_m is not None and target_resolution_m > resolution:
        factor = max(1, round(target_resolution_m / resolution))
        dem_da = dem_da.coarsen(x=factor, y=factor, boundary="trim").mean()

    return dem_da


def _auto_target_resolution(radius_km, native_resolution, max_pixels_per_side=2000):
    """
    Calcule une résolution cible (en m) telle que le DEM final ne dépasse
    jamais `max_pixels_per_side` pixels de côté, quel que soit le rayon
    demandé.

    Sans ce plafond, un grand rayon (ex. 200 km) combiné à la résolution
    native la plus fine (32 m) donne un DEM de l'ordre de 12500x12500
    pixels (~150 millions de points, plusieurs centaines de Mo en mémoire
    par figure) -- multiplié par le nombre de sous-figures d'une grille,
    c'est ce qui fait planter la machine.
    """
    diameter_m = 2 * radius_km * 1000
    min_resolution_needed = diameter_m / max_pixels_per_side
    target = max(native_resolution, min_resolution_needed)
    # arrondi à un multiple entier du natif, pour un coarsen() propre
    factor = max(1, round(target / native_resolution))
    return native_resolution * factor


def _compute_center_and_bounds(
    aws_path, dataset, dataset_variable, start_date, end_date, radius_km
):
    """
    Calcule le centre (en EPSG:3413) et l'emprise (x_min, y_min, x_max, y_max)
    autour d'une station AWS, à partir de la trajectoire satellite la plus
    proche de start_date/end_date.

    Factorisé hors de `plot_regional_variation` pour pouvoir être appelé
    UNE SEULE FOIS par rayon (et non une fois par combinaison dataset x
    rayon) dans `subplot_datasets_radius`.
    """
    df = interpolation_altimetry_AWS.satellite_on_aws(aws_path, dataset, dataset_variable)
    s_date = pd.to_datetime(start_date)
    e_date = pd.to_datetime(end_date)
    start_index = (df["time"] - s_date).abs().idxmin()
    end_index = (df["time"] - e_date).abs().idxmin()
    lat_start, lon_start = df["lat"][start_index], df["lon"][start_index]
    lat_end, lon_end = df["lat"][end_index], df["lon"][end_index]

    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    center_x, center_y = transformer.transform(
        (lon_start + lon_end) / 2, (lat_start + lat_end) / 2
    )

    radius_m = radius_km * 1000
    bounds = (
        center_x - radius_m,
        center_y - radius_m,
        center_x + radius_m,
        center_y + radius_m,
    )
    return (center_x, center_y), bounds


def add_location_inset(
    fig,
    ax,
    x_min,
    x_max,
    y_min,
    y_max,
    inset_width_frac=0.22,
    inset_height_frac=0.22,
    margin_frac=0.03,
):
    """
    Ajoute une petite carte de localisation du Groenland dans le coin
    inférieur gauche de l'axe principal `ax`, avec un rectangle indiquant
    l'emprise (x_min, x_max, y_min, y_max) en EPSG:3413.

    inset_width_frac / inset_height_frac : taille de l'inset, en fraction
        de la largeur/hauteur de l'axe principal.
    margin_frac : marge intérieure par rapport au coin de l'axe principal.
    """
    bbox = ax.get_position()  # position de l'axe principal en coordonnées figure (0-1)

    inset_w = bbox.width * inset_width_frac
    inset_h = bbox.height * inset_height_frac
    margin_x = bbox.width * margin_frac
    margin_y = bbox.height * margin_frac

    left = bbox.x0 + margin_x
    bottom = bbox.y0 + margin_y

    ax_inset = fig.add_axes([left, bottom, inset_w, inset_h], projection=EPSG3413_CRS)
    ax_inset.set_extent([-60, -30, 59, 84], crs=ccrs.PlateCarree())

    ax_inset.add_feature(
        cfeature.LAND.with_scale("50m"), facecolor="0.85", edgecolor="none", zorder=1
    )
    ax_inset.add_feature(
        cfeature.OCEAN.with_scale("50m"), facecolor="#cfe8f3", zorder=0
    )
    ax_inset.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.4, zorder=2)

    rect = mpatches.Rectangle(
        (x_min, y_min),
        x_max - x_min,
        y_max - y_min,
        linewidth=1.1,
        edgecolor="red",
        facecolor="none",
        transform=EPSG3413_CRS,
        zorder=5,
    )
    ax_inset.add_patch(rect)

    ax_inset.set_xticks([])
    ax_inset.set_yticks([])
    for spine in ax_inset.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(0.6)

    return ax_inset

def print_extrema_with_uncertainty(aws_name, radius, label, anomaly_da, uncertainty_da=None):
    """
    Affiche dans le terminal le minimum et le maximum de anomaly_da,
    accompagnés de l'incertitude au même endroit (si uncertainty_da fourni).
    """
    data = anomaly_da.values
    if np.all(np.isnan(data)):
        print(f"  {label} : aucune donnée valide dans l'emprise.")
        return

    imin = np.unravel_index(np.nanargmin(data), data.shape)
    imax = np.unravel_index(np.nanargmax(data), data.shape)
    vmin, vmax = data[imin], data[imax]

    if uncertainty_da is not None:
        u = uncertainty_da.values
        umin = u[imin] if u.shape == data.shape else np.nan
        umax = u[imax] if u.shape == data.shape else np.nan
        print(f"  {label} : min = {vmin:.3f} ± {umin:.3f} m | max = {vmax:.3f} ± {umax:.3f} m")
    else:
        print(f"  {label} : min = {vmin:.3f} m | max = {vmax:.3f} m (incertitude indisponible)")

def plot_regional_variation(
    aws,
    dataset,
    dataset_variable,
    start_date,
    end_date,
    radius_km,
    aws_name,
    label="Station",
    cross_color="black",
    cross_size=200,
    label_offset_km=(2, 2),
    label_fontsize=11,
    title=None,
    cbar_label=None,
    dem_resolution=32,
    dem_target_resolution=None,
    max_dem_pixels_per_side=2000,
    dem_alpha=0,
    var_alpha=0.75,
    vmin=-3,
    vmax=3,
    add_colorbar=True,
    show_inset=True,
    save=True,
    fig=None,
    ax=None,
    center_xy=None,
    dem=None,
):
    """
    Trace la variation spatiale d'une variable satellite autour d'une station AWS,
    superposée à un hillshade ArcticDEM.

    center_xy, dem : si fournis (ex. depuis `subplot_datasets_radius`, qui les
        précalcule une seule fois par rayon), ils sont réutilisés tels quels
        au lieu d'être recalculés/retéléchargés -- évite de refaire le même
        fetch DEM pour chaque dataset d'une même colonne (même rayon).
    max_dem_pixels_per_side : plafond de taille du DEM chargé (voir
        `_auto_target_resolution`), utilisé seulement si `dem` n'est pas
        fourni et que `dem_target_resolution` n'est pas fourni non plus.
    """
    if dem_resolution not in VALID_ARCTICDEM_RESOLUTIONS:
        raise ValueError(
            f"dem_resolution={dem_resolution} invalide. "
            f"Valeurs autorisées : {VALID_ARCTICDEM_RESOLUTIONS}."
        )

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(7, 7))

    df = interpolation_altimetry_AWS.satellite_on_aws(aws, dataset, dataset_variable)
    ds = interpolation_altimetry_AWS.satellite_opening(
        interpolation_altimetry_AWS.SATELLITE[dataset]["file"]
    )

    s_date = pd.to_datetime(start_date)
    e_date = pd.to_datetime(end_date)

    if center_xy is None:
        (center_x, center_y), bounds = _compute_center_and_bounds(
            aws, dataset, dataset_variable, start_date, end_date, radius_km
        )
    else:
        center_x, center_y = center_xy
        radius_m = radius_km * 1000
        bounds = (
            center_x - radius_m,
            center_y - radius_m,
            center_x + radius_m,
            center_y + radius_m,
        )

    x_min, y_min, x_max, y_max = bounds

    if dem is None:
        if dem_target_resolution is None:
            dem_target_resolution = _auto_target_resolution(
                radius_km, dem_resolution, max_dem_pixels_per_side
            )
        dem = fetch_arcticdem(
            bounds,
            resolution=dem_resolution,
            target_resolution_m=dem_target_resolution,
        )

    actual_start = pd.Timestamp(
        ds["time"].sel(time=s_date, method="nearest").values
    ).date()
    actual_end = pd.Timestamp(
        ds["time"].sel(time=e_date, method="nearest").values
    ).date()

    var = ds[dataset_variable].sel(time=e_date, method="nearest") - ds[
        dataset_variable
    ].sel(time=s_date, method="nearest")
    y_descending = ds["y"].values[0] > ds["y"].values[-1]
    var_zoom = var.sel(
        x=slice(x_min, x_max),
        y=slice(y_max, y_min) if y_descending else slice(y_min, y_max),
    )

    dem_vals = dem.values
    ls = LightSource(azdeg=315, altdeg=45)
    hillshade = ls.hillshade(
        np.nan_to_num(dem_vals, nan=np.nanmin(dem_vals)), vert_exag=1.5
    )

    ax.pcolormesh(
        dem["x"],
        dem["y"],
        hillshade,
        cmap="gray",
        shading="auto",
        alpha=dem_alpha,
        zorder=1,
    )

    im = var_zoom.plot.pcolormesh(
        x="x",
        y="y",
        ax=ax,
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        add_colorbar=False,
        alpha=var_alpha,
        zorder=2,
    )

    if add_colorbar:
            cbar = fig.colorbar(im, ax=ax, shrink=0.70)
            cbar.set_label(cbar_label if cbar_label is not None else f"Δ {dataset_variable}")

    ax.scatter(
        center_x,
        center_y,
        marker="+",
        s=cross_size,
        linewidths=2.2,
        color=cross_color,
        zorder=5,
    )
    ax.annotate(
        label,
        xy=(center_x, center_y),
        xytext=(
            center_x + label_offset_km[0] * 1000,
            center_y + label_offset_km[1] * 1000,
        ),
        fontsize=label_fontsize,
        fontweight="bold",
        color=cross_color,
        zorder=6,
    )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m, EPSG:3413)")
    ax.set_ylabel("y (m, EPSG:3413)")
    ax.set_title(
        title
        if title is not None
        else (
            f"Spatialized Surface Elevation Change around {aws_name} station \n"
            f"between {actual_start} and {actual_end} with a {radius_km} km radius"
        ),
        fontsize=13,
        fontweight="bold",
    )

    if standalone:
        plt.tight_layout()

    if show_inset:
        add_location_inset(fig, ax, x_min, x_max, y_min, y_max)

    if standalone:
        plt.tight_layout()
        if save:
            out_path = FIGURES_DIR / (
                f"Spatialized Surface Elevation Change around {aws_name} station "
                f"between {actual_start} and {actual_end}.png"
            )
            fig.savefig(out_path, dpi=300, bbox_inches="tight")
            print(f"Figure sauvegardée : {out_path}")
        plt.show()

    return fig, ax, im, var_zoom


def subplot_anomaly_uncertainty(
    aws_path,
    radius,
    aws_name,
    #dataset variable
    copernicus_anomaly_var="dh",
    andersen_anomaly_var="ZZ",
    nilsson_anomaly_var="dh",
    khan_anomaly_var="dh_vol",
    zhang_anomaly_var="elev_interp",
    #dataset uncertainty variable
    copernicus_uncertainty_var="dh_uncert",
    andersen_uncertainty_var="ZZer",
    nilsson_uncertainty_var="rms",
    khan_uncertainty_var=None,
    zhang_uncertainty_var="elev_uncer_interp",
    #dataset label
    copernicus_label="Copernicus_Climate_Data_Store",
    andersen_label="Andersen et al., 2025",
    nilsson_label="Nilsson and Gardner, 2026",
    khan_label="Khan et al., 2025",
    zhang_label="Zhang et al., 2022",
    start_date=None,
    end_date=None,
    dem_resolution=32,
    max_dem_pixels_per_side=2000,
    # -- échelle de la colorbar (une paire vmin/vmax par ligne) --
    anomaly_vmin=-3, anomaly_vmax=3,
    uncertainty_vmin=-3, uncertainty_vmax=3,
    auto_scale=False,
    # auto_scale=False (défaut) -> les bornes ci-dessus sont utilisées telles
    #   quelles pour toute la ligne (ex. anomaly_vmin=-10, anomaly_vmax=10).
    # auto_scale=True -> ignore anomaly_vmin/vmax et uncertainty_vmin/vmax,
    #   et recalcule pour chaque ligne des bornes symétriques [-m, m] où m
    #   est le max absolu observé sur tous les subplots de la ligne
    #   (comportement de la version précédente du script).
    save=True,
    output_path=None,
):
    anomaly_row = [
        (copernicus_label, copernicus_anomaly_var, copernicus_label),
        (andersen_label, andersen_anomaly_var, andersen_label),
        (nilsson_label, nilsson_anomaly_var, nilsson_label),
        (khan_label, khan_anomaly_var, khan_label),
        (zhang_label, zhang_anomaly_var, zhang_label),
    ]
    uncertainty_row = [
        (copernicus_label, copernicus_uncertainty_var, copernicus_label),
        (andersen_label, andersen_uncertainty_var, andersen_label),
        (nilsson_label, nilsson_uncertainty_var, nilsson_label),
        (khan_label, khan_uncertainty_var, khan_label),
        (zhang_label, zhang_uncertainty_var, zhang_label),
    ]

    n_rows, n_cols = 2, 5
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows), squeeze=False,
        constrained_layout=True,
    )

    ref_dataset, ref_variable, _ = anomaly_row[0]
    center_xy, bounds = _compute_center_and_bounds(
        aws_path, ref_dataset, ref_variable, start_date, end_date, radius
    )
    target_res = _auto_target_resolution(radius, dem_resolution, max_dem_pixels_per_side)
    dem = fetch_arcticdem(bounds, resolution=dem_resolution, target_resolution_m=target_res)

    def _plot_row(row_axes, row_items, row_kind, row_vmin, row_vmax, row_auto_scale):
        images, plotted_axes, var_zooms = [], [], {}
        for ax, (dataset, dataset_variable, subtitle) in zip(row_axes, row_items):
            if dataset_variable is None:
                ax.axis("off")
                ax.text(
                    0.5, 0.5,
                    f"No {row_kind} data\navailable for\n{subtitle}",
                    ha="center", va="center", fontsize=11, style="italic",
                    transform=ax.transAxes,
                )
                continue

            _, _, im, var_zoom = plot_regional_variation(
                aws_path,
                dataset,
                dataset_variable,
                start_date,
                end_date,
                radius,
                aws_name,
                label=aws_name,
                var_alpha=0.75,
                vmin=row_vmin,
                vmax=row_vmax,
                add_colorbar=False,
                show_inset=False,
                title=f"{subtitle}\n({row_kind})",
                fig=fig,
                ax=ax,
                center_xy=center_xy,
                dem=dem,
                save=False,
            )
            images.append(im)
            plotted_axes.append(ax)
            var_zooms[subtitle] = var_zoom

        if images:
            if row_auto_scale:
                # bornes symétriques recalculées à partir du max absolu de la ligne
                row_max_abs = 0.0
                for im in images:
                    data = np.ma.filled(im.get_array().astype(float), np.nan)
                    row_max_abs = max(row_max_abs, np.nanmax(np.abs(data)))
                for im in images:
                    im.set_clim(-row_max_abs, row_max_abs)
                cbar_vmin, cbar_vmax = -row_max_abs, row_max_abs
            else:
                # bornes manuelles, appliquées telles quelles à tous les subplots de la ligne
                for im in images:
                    im.set_clim(row_vmin, row_vmax)
                cbar_vmin, cbar_vmax = row_vmin, row_vmax

            cbar = fig.colorbar(images[0], ax=row_axes, shrink=0.75, location="right")
            cbar.set_label(f"Δ SEC {row_kind}")

        return plotted_axes, var_zooms

    plotted_axes = []
    axes_anomaly, zooms_anomaly = _plot_row(
        axes[0], anomaly_row, "anomaly", anomaly_vmin, anomaly_vmax, auto_scale
    )
    axes_uncert, zooms_uncert = _plot_row(
        axes[1], uncertainty_row, "uncertainty", uncertainty_vmin, uncertainty_vmax, auto_scale
    )
    plotted_axes = axes_anomaly + axes_uncert

    fig.suptitle(
        f"{aws_name} — SEC anomaly and uncertainties ({end_date} - {start_date}), radius {radius} km",
        fontsize=15,
        fontweight="bold",
        y=1.02,
    )

    fig.canvas.draw()

    x_min, y_min, x_max, y_max = (
        center_xy[0] - radius * 1000,
        center_xy[1] - radius * 1000,
        center_xy[0] + radius * 1000,
        center_xy[1] + radius * 1000,
    )
    for ax in plotted_axes:
        add_location_inset(fig, ax, x_min, x_max, y_min, y_max)

    if save:
        if output_path is None:
            output_path = (
                FIGURES_DIR
                / f"anomaly_uncertainty_grid_{aws_name}_{start_date}_{end_date}_{radius}km.png"
            )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Figure sauvegardée : {output_path}")

    plt.show()

    print(f"\n--- Extrema — {aws_name}, rayon {radius} km ({start_date} → {end_date}) ---")
    for label, anomaly_da in zooms_anomaly.items():
        uncertainty_da = zooms_uncert.get(label)
        print_extrema_with_uncertainty(aws_name, radius, label, anomaly_da, uncertainty_da)

    return fig, axes




for station, info in AWS_data.STATION_UPE_L.items():
    subplot_anomaly_uncertainty(
        info["file"],
        radius=75,
        aws_name=station,
        start_date="2020-04-01",
        end_date="2021-09-01",
        auto_scale=True,
        # anomaly_vmin=-3, anomaly_vmax=3,
        # uncertainty_vmin=-3, uncertainty_vmax=3,
    )