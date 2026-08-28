import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import ruptures as rpt
from scipy import stats
from pathlib import Path

import AWS_data
import interpolation_altimetry_AWS
import plot_time_series_altimetry_AWS

from paths import FIGURES_DIR as _BASE_FIGURES_DIR

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Small numeric/date helpers ───────────────────────────────────────────
"""
def _to_scalar(v) -- coerces a value to float, returning NaN if it can't be converted.

def _as_list(x) -- forces x into a list if it isn't already one (so satellite_name / variable / uncertainty_variable can be passed as either a single value or a list).

def _decimal_years(time_index) -- converts a DatetimeIndex to decimal years relative to its first point, so a linear-regression slope comes out directly in units/year.

def _years_to_dates(x_years, time_index) -- converts decimal years (same origin as _decimal_years) back to dates, to plot segments on the time axis.

def _dates_to_indices(dates, time_index) -- converts a list of dates to the nearest matching indices (positions) in time_index.
"""

# ── Trend-segment / breakpoint detection ──────────────────────────────────
"""
def _detect_trend_segments(x_years, y, max_breakpoints, min_segment_size, penalty, verbose) -- detects slope breakpoints in (x_years, y) via piecewise-linear regression (ruptures' Pelt/Dynp), auto-selecting the number of segments up to max_breakpoints.

def _segment_slope(x_years, y, i0, i1, cm_factor) -- linear regression on segment [i0:i1); returns (slope in mm/year, x_seg in years, y_fit), or None if the segment is too short.

def _fixed_trend_segments(n, candidate_idx, min_segment_size) -- segments bounded EXACTLY by candidate_idx (no selection, no penalty): every given date becomes a segment boundary.

def _segment_ssr(x_years, y, i0, i1) -- sum of squared residuals of a linear regression on [i0:i1).

def _detect_trend_segments_from_candidates(x_years, y, candidate_idx, min_segment_size, penalty) -- selects, from candidate_idx only, the optimal subset of breakpoints (dynamic programming, cost = sum of per-segment SSR + penalty per added breakpoint).
"""

# ── Main figure ───────────────────────────────────────────────────────────
"""
def plot_variable_trend_breakpoints(STATIONS, fig_name, source, satellite_name, variable, uncertainty_variable, aws_variable, ..., detect_breakpoints, breakpoint_dates, breakpoint_dates_mode, max_breakpoints, min_segment_years, breakpoint_penalty, ...) -- for each station, plots one or more series (satellite(s) or AWS) with their uncertainty band, detects slope breakpoints (or uses fixed/candidate breakpoint_dates) and overlays per-segment trend lines annotated with their slope. source selects "satellite" (satellite_name/variable/uncertainty_variable, each scalar or list to overlay several products) or "aws" (aws_variable, from monthly_std_hourly).
"""

# ── Paths ─────────────────────────────────────────────────────────────────
"""
FIGURES_DIR, DEFAULT_SATELLITE_COLORS -- output folder (subfolder of paths.py's shared FIGURES_DIR) and the default per-satellite colour mapping used across this module's figures.
"""


FIGURES_DIR = _BASE_FIGURES_DIR / "TUN minimum Andersen et al"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Couleurs par défaut par satellite, à aligner avec celles utilisées dans
# tes autres scripts (plots_variables, plots_taylor, ...). Modifie ce dict
# si tu veux changer une couleur globalement, ou passe `satellite_colors=`
# à plot_variable_trend_breakpoints pour ne changer que pour un appel donné.
DEFAULT_SATELLITE_COLORS = {
    "Copernicus_Climate_Data_Store": "#1f77b4",
    "Nilsson and Gardner, 2026": "#ff7f0e",
    "Andersen et al., 2025": "#2ca02c",
    "Khan et al., 2025": "#d62728",
    "Zhang et al., 2022": "#9467bd",
}

def _to_scalar(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _as_list(x):
    """Force x en liste si ce n'est pas déjà une liste/tuple (pour permettre
    de passer soit une valeur scalaire, soit une liste de valeurs, pour
    satellite_name / variable / uncertainty_variable)."""
    return list(x) if isinstance(x, (list, tuple)) else [x]


def _decimal_years(time_index):
    """Convertit un DatetimeIndex en années décimales (relatif au 1er point),
    pour que la pente d'une régression linéaire soit directement en
    unité/an."""
    t0 = pd.Timestamp(time_index[0])
    seconds = np.array([(pd.Timestamp(t) - t0).total_seconds() for t in time_index])
    return seconds / (365.25 * 24 * 3600)


def _years_to_dates(x_years, time_index):
    """Reconvertit des années décimales (même origine que _decimal_years)
    en dates, pour pouvoir tracer les segments sur l'axe temporel."""
    t0 = pd.Timestamp(time_index[0])
    return pd.DatetimeIndex([t0 + pd.Timedelta(days=365.25 * x) for x in x_years])


def _detect_trend_segments(x_years, y, max_breakpoints=3, min_segment_size=10,
                            penalty=None, verbose=False):
    """
    Détecte les ruptures de pente dans (x_years, y) via une régression
    linéaire par morceaux (ruptures.Pelt / Dynp, model='linear').

    Renvoie une liste de tuples (i_start, i_end) (i_end exclu) délimitant
    chaque segment. Si `ruptures` est absent ou la série trop courte,
    renvoie un seul segment couvrant toute la série.
    """
    n = len(y)
    if n < 2 * min_segment_size:
        return [(0, n)]

    signal = np.column_stack([x_years, y])
    # ruptures propose plusieurs façon de mesurer si un segment est homogèn :
    #       "l2" déctecte des ruptures de moyenne
    #       "linear" détecte des ruptures dans une relation de régression entre une ou plusieurs colonnes "régresseurs" et une colonne "cible".
    # Pour 'signal', ruptures impose une convention : 1ère colonne = variable expliquée et 2ème colonne = régresseur (ici le temps).
    # Ruptures fait alors une régression sur le temps et cherche où la pente change (coeff de régression temps->valeur)

    if penalty is None:
        # pénalité type BIC par défaut, recommandée par la doc de ruptures pour Pelt
        # BIC : Bayesian Information Criterion : formule : BIC = -2·log(vraisemblance du modèle) + p·log(n)
        # Pelt : Pruned Exact Linear Time : algo de recherche (exacte) qui trouve le découpage optimal pour la pénalité donnée tout en restant rapide.
        sigma2 = np.var(y)
        penalty = np.log(n) * sigma2 if sigma2 > 0 else 1.0

    algo = rpt.Pelt(model="linear", min_size=min_segment_size, jump=1).fit(signal)
        # model : cf précédemment
        # min_size : aucun segment ne peut faire moins que min_segment_size
        # jump=1 : l'aglo teste toutes les positions possibles (pas de sous-échantillonnage)
    bkps = algo.predict(pen=penalty)
        # renvoie la liste des indices de ruptures : par convention le dernier vaut toujours n (fin de la série, rupture en soi). Ce sont les bornes.

    n_bkps_pelt = len(bkps) - 1
    if verbose:
        print(f"[_detect_trend_segments] Pelt (pen={penalty:.4g}) -> "
              f"{n_bkps_pelt} rupture(s) trouvée(s) "
              f"(plafond max_breakpoints={max_breakpoints})")

    if n_bkps_pelt > max_breakpoints:
        # Si trop de ruptures détectées avec la pénalité auto : on fait tourner avec un autre algorithme Dynp : Dynamic Programming
        # On impose cette fois_ci un nombre fixe de segments (n_bkps=max_breakpoints)
        if verbose:
            print(f"[_detect_trend_segments] -> repli sur Dynp avec "
                  f"n_bkps={max_breakpoints} (penalty ignorée à partir d'ici)")

        algo2 = rpt.Dynp(model="linear", min_size=min_segment_size, jump=1).fit(signal)
        bkps = algo2.predict(n_bkps=max_breakpoints)

    bounds = [0] + list(bkps)  # Ajoute 0 à la liste des bornes de ruptures
    segments = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]
    # on élimine d'éventuels segments vides/trop courts issus de l'algo
    segments = [(a, b) for (a, b) in segments if b - a >= 2]
    return segments if segments else [(0, n)]


def _segment_slope(x_years, y, i0, i1, cm_factor=1000.0):
    """Régression linéaire sur le segment [i0:i1). Renvoie
    (slope_mm_par_an, x_seg (années), y_fit) ou None si segment trop court."""
    xs = x_years[i0:i1]
    ys = y[i0:i1]
    if len(xs) < 2:
        return None
    slope, intercept, r, p, se = stats.linregress(xs, ys)
    y_fit = slope * xs + intercept
    return slope * cm_factor, xs, y_fit


def _dates_to_indices(dates, time_index):
    """Convertit une liste de dates en indices (positions) les plus proches
    dans time_index (DatetimeIndex de la série déjà filtrée/temporelle)."""
    idx = pd.DatetimeIndex([pd.to_datetime(d) for d in dates])
    positions = time_index.get_indexer(idx, method="nearest")
    return sorted(set(int(p) for p in positions if p != -1))


def _fixed_trend_segments(n, candidate_idx, min_segment_size):
    """Segments délimités EXACTEMENT par candidate_idx (pas de sélection,
    pas de pénalité) : chaque date fournie devient une frontière de segment."""
    bounds = sorted(set([0] + [c for c in candidate_idx if 0 < c < n] + [n]))
    segments = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]
    segments = [(a, b) for a, b in segments if b - a >= 2]
    return segments if segments else [(0, n)]


def _segment_ssr(x_years, y, i0, i1):
    """Somme des carrés des résidus d'une régression linéaire sur [i0:i1)."""
    xs, ys = x_years[i0:i1], y[i0:i1]
    if len(xs) < 2:
        return 0.0
    slope, intercept, r, p, se = stats.linregress(xs, ys)
    y_fit = slope * xs + intercept
    return float(np.sum((ys - y_fit) ** 2))


def _detect_trend_segments_from_candidates(x_years, y, candidate_idx, min_segment_size, penalty):
    """
    Sélectionne, PARMI candidate_idx uniquement, le sous-ensemble optimal de
    ruptures (programmation dynamique, coût = somme des SSR par segment +
    penalty par rupture ajoutée) -- même principe que Pelt, mais l'espace de
    recherche est restreint aux dates fournies au lieu de tous les points.
    """
    n = len(y)
    candidates = sorted(set(
        [0] + [c for c in candidate_idx if min_segment_size <= c <= n - min_segment_size] + [n]
    ))
    m = len(candidates)
    if m <= 2:
        return [(0, n)]

    INF = float("inf")
    cost = [INF] * m
    back = [-1] * m
    cost[0] = 0.0

    for j in range(1, m):
        for i in range(j):
            if candidates[j] - candidates[i] < min_segment_size:
                continue
            c = cost[i] + _segment_ssr(x_years, y, candidates[i], candidates[j])
            if i > 0:  # pénalité seulement pour les ruptures internes, pas les bornes 0/n
                c += penalty
            if c < cost[j]:
                cost[j] = c
                back[j] = i

    bounds, idx = [], m - 1
    while idx >= 0:
        bounds.append(candidates[idx])
        idx = back[idx]
    bounds = sorted(set(bounds))
    segments = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]
    segments = [(a, b) for a, b in segments if b - a >= 2]
    return segments if segments else [(0, n)]

# ─────────────────────────────────────────────────────────────────────────
# Fonction principale
# ─────────────────────────────────────────────────────────────────────────

def plot_variable_trend_breakpoints(
        STATIONS, fig_name,
        source="satellite",                 # "satellite" ou "aws"
        # -- source satellite -- (scalaire OU liste pour superposer plusieurs satellites)
        satellite_name=None, variable=None, uncertainty_variable=None,
        # -- source aws --
        aws_variable=None, aws_hourly_file_key="hourly_data",
        aws_min_valid_fraction=0.5,
        # -- commun --
        color="#1f4e79",
        series_colormap="tab10",            # repli si un satellite n'est pas dans satellite_colors
        satellite_colors=None,              # dict {satellite_name: couleur}, comme dans les anciens scripts
        ref_date=None,
        start_date='2003-01-01', end_date=None,
        band_alpha=0.10,
        show_uncertainty=True,
        value_unit_to_cm_factor=100.0,     # 100 si la série est en m
        # -- détection de ruptures de pente --
        detect_breakpoints=True,
        breakpoint_dates=None,          # liste de dates (str ou Timestamp) -> ignore detect_breakpoints
        breakpoint_dates_mode="fixed",  # "fixed" (dates imposées) ou "candidates" (sélection parmi ces dates)
        max_breakpoints=3,
        min_segment_years=1.0,              # taille min. d'un segment (en années)
        breakpoint_penalty=None,            # None -> pénalité BIC auto
        segment_colormap="tab10",
        segment_linewidth=1.3,
        slope_unit='cm/y',
        slope_fontsize=6,
        annotate_slopes=True,
        ):
    """
    Trace, pour chaque station, une ou plusieurs séries (satellite(s) ou
    AWS) avec leur bandeau d'incertitude, détecte les ruptures de pente et
    superpose les droites de tendance par segment, avec la valeur de la
    pente (en mm/an par défaut) annotée en bout de segment.

    -- Choix de la/des série(s) --
    source : "satellite" ou "aws"
    satellite_name, variable, uncertainty_variable : si source="satellite".
        Chacun peut être une valeur scalaire (1 seule série) OU une liste
        (plusieurs séries superposées sur le même subplot, une couleur par
        série). Les 3 listes doivent alors avoir la même longueur.
        uncertainty_variable peut valoir None pour une série donnée si
        aucune incertitude n'est disponible pour ce satellite.
    aws_variable : "z_surf_combined_anomaly" ou "gps_alt_anomaly", si
        source="aws".
    aws_hourly_file_key : clé dans info (dict STATIONS) pointant vers le
        fichier horaire utilisé par monthly_std_hourly (source="aws").

    -- Commun --
    ref_date : str ou None. Si fourni, chaque série est recalée à 0 à la
        date la plus proche de ref_date. Sinon, elle part de 0 à son
        premier point.
    start_date / end_date : bornes temporelles du plot.
    value_unit_to_cm_factor : facteur de conversion de l'unité de la
        série vers le mm (1000 si la série est en mètres), utilisé
        uniquement pour l'affichage de la pente.
    series_colormap : colormap utilisée pour distinguer les séries
        (satellites) entre elles, quand il y en a plusieurs. Si une seule
        série, `color` est utilisé tel quel.

    -- Détection de ruptures --
    detect_breakpoints : bool, active la détection (sinon un seul segment
        = régression linéaire globale).
    max_breakpoints : nombre max. de ruptures autorisées par station.
    min_segment_years : longueur minimale (en années) d'un segment, pour
        éviter des micro-segments non significatifs.
    breakpoint_penalty : pénalité Pelt (ruptures). None -> calculée
        automatiquement (type BIC) à partir de la variance du signal.
    segment_colormap : colormap matplotlib utilisée pour distinguer les
        segments de RUPTURE (une couleur différente par segment). Si
        plusieurs séries sont superposées, on préfère souvent
        annotate_slopes=False / detect_breakpoints=False pour ne pas
        surcharger le subplot.
    slope_unit : chaîne affichée à côté de la valeur de pente (mm/an par
        défaut ; change juste le texte, pas le calcul).
    """

    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
    N = len(STATIONS)
    NCOLS = 2
    NROWS = (N + NCOLS - 1) // NCOLS

    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()

    cmap = plt.get_cmap(segment_colormap)
    series_cmap = plt.get_cmap(series_colormap)

    # fusion : couleurs passées en argument > couleurs par défaut du module
    resolved_satellite_colors = dict(DEFAULT_SATELLITE_COLORS)
    if satellite_colors:
        if isinstance(satellite_colors, dict):
            resolved_satellite_colors.update(satellite_colors)
        else:
            # liste/tuple de couleurs, alignée avec satellite_name (même ordre)
            sat_names_for_colors = _as_list(satellite_name) if source == "satellite" else []
            colors_list = list(satellite_colors)
            if len(colors_list) != len(sat_names_for_colors):
                raise ValueError(
                    "satellite_colors (liste) doit avoir la même longueur que "
                    "satellite_name pour être associée dans le même ordre."
                )
            resolved_satellite_colors.update(dict(zip(sat_names_for_colors, colors_list)))

    # ── Normalisation satellite_name / variable / uncertainty_variable en listes ──
    if source == "satellite":
        sat_names_all = _as_list(satellite_name)
        variables_all = _as_list(variable)
        uncs_all = _as_list(uncertainty_variable)
        if not (len(sat_names_all) == len(variables_all) == len(uncs_all)):
            raise ValueError(
                "satellite_name, variable et uncertainty_variable doivent "
                "avoir la même longueur (une entrée par série à superposer)."
            )
        series_specs = list(zip(sat_names_all, variables_all, uncs_all))
    else:
        series_specs = [(None, None, None)]  # source == "aws" : une seule série

    multi_series = len(series_specs) > 1

    fig, axes = plt.subplots(
        NROWS, NCOLS,
        figsize=(14, NROWS * 3.5),
        constrained_layout=True,
    )
    axes_flat = np.atleast_1d(axes).ravel()

    for ax, (station, info) in zip(axes_flat, STATIONS.items()):

        any_data_plotted = False

        for series_idx, (sat_name_i, var_i, unc_i) in enumerate(series_specs):

            if multi_series:
                if sat_name_i in resolved_satellite_colors:
                    series_color = resolved_satellite_colors[sat_name_i]
                else:
                    series_color = series_cmap(series_idx % series_cmap.N)
            else:
                series_color = resolved_satellite_colors.get(sat_name_i, color)
            series_label = sat_name_i if sat_name_i is not None else (aws_variable or "AWS")

            # ── Chargement de la série choisie ──────────────────────────
            if source == "satellite":
                Sat = interpolation_altimetry_AWS.satellite_on_aws(info['file'], sat_name_i, var_i)
                series = Sat[var_i].dropna()
                series1 = Sat.iloc[series.index.tolist()]
                time_full = pd.to_datetime(series1['time_sat'])
                values_full = series.values
            elif source == "aws":
                aws_funcs = {
                    "z_surf_combined_anomaly": plot_time_series_altimetry_AWS.z_surf_combined_anomaly,
                    "gps_alt_anomaly": plot_time_series_altimetry_AWS.gps_alt_anomaly,
                }
                aws_df = aws_funcs[aws_variable](info['file'])
                aws_series = aws_df[aws_variable].dropna()
                time_full = pd.DatetimeIndex(aws_df.loc[aws_series.index, "time"])
                values_full = aws_series.values
            else:
                raise ValueError("source doit être 'satellite' ou 'aws'")

            # ── Filtrage temporel ────────────────────────────────────────
            mask = pd.Series(True, index=range(len(time_full)))
            if start_bound is not None:
                mask &= (time_full >= start_bound).values
            if end_bound is not None:
                mask &= (time_full <= end_bound).values
            time_sel = pd.DatetimeIndex(time_full[mask.values])
            values_sel = np.asarray(values_full)[mask.values]

            if len(values_sel) == 0:
                continue  # pas de données pour cette série sur cette station

            any_data_plotted = True

            # ── Alignement sur ref_date ──────────────────────────────────
            if ref_date is not None:
                ref = pd.to_datetime(ref_date)
                closest_idx = (time_sel - ref).abs().argmin()
                baseline = values_sel[closest_idx]
            else:
                baseline = values_sel[0]
            values_aligned = values_sel - baseline

            # ── Bandeau d'incertitude ─────────────────────────────────────
            if show_uncertainty:
                if source == "satellite" and unc_i is not None:
                    Sat_unc = interpolation_altimetry_AWS.satellite_on_aws(info['file'], sat_name_i, unc_i)
                    unc_raw = Sat_unc[unc_i]
                    unc_values = np.array([_to_scalar(v) for v in unc_raw.values], dtype="float64")
                    unc_time = pd.to_datetime(Sat_unc['time_sat'])

                    unc_df = pd.DataFrame({"time": unc_time.values, "unc": unc_values}).dropna()
                    unc_df["time"] = unc_df["time"].astype("datetime64[ns]")
                    main_df = pd.DataFrame({"time": time_sel.values, "value": values_aligned})
                    main_df["time"] = main_df["time"].astype("datetime64[ns]")
                    merged = pd.merge(main_df, unc_df, on="time", how="inner")

                    if not merged.empty:
                        ax.fill_between(
                            merged["time"], merged["value"] - merged["unc"],
                            merged["value"] + merged["unc"],
                            color=series_color, alpha=band_alpha, linewidth=0, zorder=1,
                        )

                elif source == "aws":
                    hourly_path = info.get(aws_hourly_file_key)
                    if hourly_path is not None:
                        base_col = plot_time_series_altimetry_AWS.AWS_BASE_COLUMN.get(aws_variable)
                        std_df = plot_time_series_altimetry_AWS.monthly_std_hourly(
                            hourly_path, columns=(base_col,),
                            min_valid_fraction=aws_min_valid_fraction,
                        )
                        sigma_lookup = std_df[base_col]
                        sigma_lookup.index = sigma_lookup.index.to_period("M")
                        day_periods = time_sel.to_period("M")
                        sigma_values = sigma_lookup.reindex(day_periods).values

                        valid_sigma = ~pd.isna(sigma_values)
                        if valid_sigma.any():
                            t_arr = time_sel[valid_sigma]
                            v_arr = values_aligned[valid_sigma]
                            s_arr = sigma_values[valid_sigma].astype(float)
                            ax.fill_between(
                                t_arr, v_arr - s_arr, v_arr + s_arr,
                                color=series_color, alpha=band_alpha, linewidth=0, zorder=1,
                            )

            # ── Courbe brute ───────────────────────────────────────────────
            ax.plot(time_sel, values_aligned, linewidth=1.2, color=series_color,
                     alpha=0.65, zorder=3, label=series_label)

            # ── Détection des ruptures de pente + régressions par segment ──
            x_years = _decimal_years(time_sel)
            n_pts = len(x_years)
            span_years = x_years[-1] - x_years[0] if n_pts > 1 else 0
            avg_dt_years = span_years / max(n_pts - 1, 1) if n_pts > 1 else 1
            min_segment_size = max(3, int(min_segment_years / avg_dt_years)) if avg_dt_years > 0 else 3

            if breakpoint_dates is not None:
                candidate_idx = _dates_to_indices(breakpoint_dates, time_sel)
                if breakpoint_dates_mode == "fixed":
                    segments = _fixed_trend_segments(n_pts, candidate_idx, min_segment_size)
                elif breakpoint_dates_mode == "candidates":
                    pen = breakpoint_penalty
                    if pen is None:
                        sigma2 = np.var(values_aligned)
                        pen = np.log(n_pts) * sigma2 if sigma2 > 0 else 1.0
                    segments = _detect_trend_segments_from_candidates(
                        x_years, values_aligned, candidate_idx, min_segment_size, pen,
                    )
                else:
                    raise ValueError("breakpoint_dates_mode doit être 'fixed' ou 'candidates'")

            elif detect_breakpoints:
                segments = _detect_trend_segments(
                    x_years, values_aligned,
                    max_breakpoints=max_breakpoints,
                    min_segment_size=min_segment_size,
                    penalty=breakpoint_penalty,
                )
            else:
                segments = [(0, n_pts)]

            for si, (i0, i1) in enumerate(segments):
                result = _segment_slope(x_years, values_aligned, i0, i1,
                                         cm_factor=value_unit_to_cm_factor)
                if result is None:
                    continue
                slope_mm_yr, xs_seg, y_fit = result
                # en mode multi-séries, on garde la couleur de la série pour
                # les segments (sinon ça devient illisible) ; sinon une
                # couleur par segment via segment_colormap comme avant.
                seg_color = series_color if multi_series else cmap(si % cmap.N)
                dates_seg = _years_to_dates(xs_seg, time_sel)

                ax.plot(dates_seg, y_fit, linestyle="--", linewidth=segment_linewidth,
                         color=seg_color, zorder=4)

                if annotate_slopes:
                    ax.annotate(
                        f"{slope_mm_yr:+.1f} {slope_unit}",
                        xy=(dates_seg[-1], y_fit[-1]),
                        xytext=(4, 4), textcoords="offset points",
                        fontsize=slope_fontsize, color=seg_color, fontweight="bold",
                        zorder=5,
                    )

                # marqueur de rupture (sauf au tout dernier point de la série)
                if i1 < n_pts:
                    ax.axvline(dates_seg[-1], color="grey", linestyle=":",
                               linewidth=0.8, alpha=0.6, zorder=0)

        if not any_data_plotted:
            ax.set_title(f"{station} (pas de données)", fontsize=12)
            continue

        ax.set_title(station, fontsize=12, fontweight="bold")
        ax.set_xlabel("Years", fontsize=9)
        ax.set_ylabel('Surface elevation change (m)', fontsize=9)
        ax.tick_params(axis="both", labelsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_major_locator(mdates.YearLocator(1))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
        ax.legend(fontsize=6, loc="best", framealpha=0.6)

    for ax in axes_flat[N:]:
        ax.set_visible(False)

    fig.suptitle(fig_name, fontsize=15, fontweight="bold")
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
    print(f"Figure sauvegardée : {OUTPUT_FILE}")
    plt.show()


# Dates imposées comme ruptures, pour toutes les stations, en superposant
# les 5 jeux de données satellite sur chaque subplot :
""" plot_variable_trend_breakpoints(
    AWS_data.STATION_SWC,
    'Linear trend and slope ruptures - Copernicus Climate Data Store at SWC',
    source="satellite",
    satellite_name=['Copernicus_Climate_Data_Store'], #, 'Nilsson and Gardner, 2026', 'Andersen et al., 2025', 'Khan et al., 2025', 'Zhang et al., 2022'],
    satellite_colors=["#EB2F25"], #, "#17EE17", "#EDE20D", "#D20DF0", "#03FEF1"],
    variable=['dh'], #, 'dh', 'ZZ', 'dh_vol', 'elev_interp'],
    uncertainty_variable=['dh_uncert'], #, 'rms', 'ZZer', None, 'elev_uncer_interp'],
    breakpoint_dates=["2011-01-01", "2013-03-13", "2014-02-01"],
    breakpoint_dates_mode="fixed",
    annotate_slopes=True,   # évite la surcharge visuelle avec 5 séries x segments
) """

""" plot_variable_trend_breakpoints(
    AWS_data.STATION_SWC,
    'Linear trend and slope ruptures - Nilsson and Gardner, 2026 at SWC',
    source="satellite",
    satellite_name=['Nilsson and Gardner, 2026'], #, 'Nilsson and Gardner, 2026', 'Andersen et al., 2025', 'Khan et al., 2025', 'Zhang et al., 2022'],
    satellite_colors=["#EB2F25"], #, "#17EE17", "#EDE20D", "#D20DF0", "#03FEF1"],
    variable=['dh'], #, 'dh', 'ZZ', 'dh_vol', 'elev_interp'],
    uncertainty_variable=['rms'], #, 'rms', 'ZZer', None, 'elev_uncer_interp'],
    breakpoint_dates=["2011-01-01", "2013-03-13", "2014-02-01"],
    breakpoint_dates_mode="fixed",
    annotate_slopes=True,   # évite la surcharge visuelle avec 5 séries x segments
) """

# plot_variable_trend_breakpoints(
#     AWS_data.STATION_SWC,
#     'Linear trend and slope ruptures - Nilsson and Gardner, 2026 at SWC',
#     source="satellite",
#     satellite_name=['Copernicus_Climate_Data_Store', 'Nilsson and Gardner, 2026', 'Andersen et al., 2025', 'Khan et al., 2025', 'Zhang et al., 2022'],
#     satellite_colors=["#EB2F25", "#17EE17", "#EDE20D", "#D20DF0", "#03FEF1"],
#     variable=['dh', 'dh', 'ZZ', 'dh_vol', 'elev_interp'],
#     uncertainty_variable=['dh_uncert' , 'rms', 'ZZer', None, 'elev_uncer_interp'],
#     breakpoint_dates=["2011-01-01", "2013-03-13", "2014-02-01"],
#     breakpoint_dates_mode="fixed",
#     annotate_slopes=True,   # évite la surcharge visuelle avec 5 séries x segments
# )

plot_variable_trend_breakpoints(
    AWS_data.STATION_use,
    'Linear trend and slope ruptures - all datasets at SDL and NSE',
    source="satellite",
    satellite_name=['Copernicus_Climate_Data_Store'], #'Copernicus_Climate_Data_Store', 'Nilsson and Gardner, 2026', 'Andersen et al., 2025', 'Khan et al., 2025', ],
    satellite_colors=["#EB2F25"], #""#03FEF1",#EB2F25", "#17EE17", "#EDE20D", "#D20DF0", ],
    variable=[ 'dh'], #'elev_interp',, 'dh', 'ZZ', 'dh_vol', 'elev_interp'],
    uncertainty_variable=[ 'dh_uncert'], #'elev_uncer_interp','dh_uncert', 'rms', 'ZZer', None, ],
    breakpoint_dates=["2013-01-01", "2014-01-01"],
    breakpoint_dates_mode="fixed",
    annotate_slopes=True,   # évite la surcharge visuelle avec 5 séries x segments
    start_date = '2010-01-01'
)

# Série AWS (z_surf_combined_anomaly), une seule tendance globale :
# plot_variable_trend_breakpoints(
#     AWS_data.STATIONS_accumulation,
#     'Tendance AWS z_surf_combined',
#     source="aws",
#     aws_variable="z_surf_combined_anomaly",
#     color="#000000",
#     detect_breakpoints=False,
# )