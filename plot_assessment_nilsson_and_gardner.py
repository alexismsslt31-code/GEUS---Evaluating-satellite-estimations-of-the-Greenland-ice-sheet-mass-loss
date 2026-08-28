import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess
import interpolation_altimetry_AWS


import AWS_data
import interpolation_altimetry_AWS

from paths import FIGURES_DIR as _BASE_FIGURES_DIR

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Multi-station overview plot ──────────────────────────────────────────
"""
def plots_variables(STATIONS, fig_name, list_satellites, variables, colors, aws_variables, aws_colors, ref_date, start_date, end_date, rms_variables, rms_alpha) -- grid of subplots, one per station, each showing the requested satellite variable(s) and optionally the matching AWS variable(s), all re-anchored to a common ref_date if given, with an optional RMS/uncertainty band per satellite.
"""

# ── Variance change test around a cutoff date ────────────────────────────
"""
def plot_variance_change_test(STATIONS, fig_name, satellite, variable, cutoff_date, lowess_frac, color_before, color_after) -- for each station, tests whether the residual (LOWESS-detrended) variance of a variable changes significantly before/after cutoff_date (Levene's test), with a two-panel figure (raw series + trend on top, coloured before/after residuals with +/-1 std bands on the bottom).
"""

# ── Trend comparison across satellites ───────────────────────────────────
"""
def plot_trend_comparison(STATIONS, fig_name, list_satellites, variables, colors, lowess_frac, output_filename, start_date, end_date) -- for each station, plots the raw series (thin, transparent) and its LOWESS trend (solid) for one or more satellites, to visually compare their respective trends.
"""

# ── Paths ─────────────────────────────────────────────────────────────────
"""
FIGURES_DIR -- output folder for this module's figures, imported from paths.py's shared FIGURES_DIR.
"""


FIGURES_DIR = _BASE_FIGURES_DIR
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def plots_variables(STATIONS, fig_name, list_satellites, variables, colors,
                     aws_variables=None, aws_colors=None, ref_date=None,
                     start_date=None, end_date=None,
                     rms_variables=None, rms_alpha=0.25):
    
    """STATIONS : dict{dict} , fig_name : str , list_satellites : list[str] ,
    variables : list[str] , colors : list[str] , ref_date : str ou None
    aws_variables : list[str] ou None, parmi 'z_surf_combined_anomaly' et/ou 'gps_alt_anomaly'
    aws_colors : list[str] ou None, couleurs associées à aws_variables (même ordre)
    Si ref_date est fourni (ex: '2010-01-15'), chaque série (satellite ET AWS) est recalée
    pour valoir 0 à la date la plus proche de ref_date. Si None, chaque série part de 0
    à son propre premier point.
    start_date / end_date : bornes temporelles du plot (str ou None).
    rms_variables : list[str ou None] ou None, même longueur que list_satellites.
    Pour chaque satellite, nom de la variable RMS/incertitude associée (ex: 'uncertainty'
    pour Yang) à afficher en bandeau (dh ± rms) autour de la courbe. Mettre None pour un
    satellite donné si aucune bande ne doit être tracée pour lui.
    rms_alpha : transparence du bandeau (0 à 1)."""

    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
    N = len(STATIONS)
    NCOLS = 3
    NROWS = (N + NCOLS - 1) // NCOLS

    # ── Bornes temporelles ──────────────────────────────────────────────────
    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()

    fig, axes = plt.subplots(
        NROWS, NCOLS,
        figsize=(14, NROWS * 3.5),
        constrained_layout=True,
    )
    axes_flat = axes.flatten()
    for ax, (station, info) in zip(axes_flat, STATIONS.items()):
        legend_labels = []
        # ── Satellites ──────────────────────────────────────────────────────
        for k in range(len(list_satellites)):
            Satellite = interpolation_altimetry_AWS.satellite_on_aws(info['file'], list_satellites[k], variables[k])
            series = Satellite[variables[k]].dropna()
            series1 = Satellite.iloc[series.index.tolist()]
            time_sat = pd.to_datetime(series1['time_sat'])
            values = series.values

            # filtrage temporel
            mask = pd.Series(True, index=range(len(time_sat)))
            if start_bound is not None:
                mask &= (time_sat >= start_bound).values
            if end_bound is not None:
                mask &= (time_sat <= end_bound).values
            time_sat = time_sat[mask.values]
            values = values[mask.values]

            if len(values) == 0:
                continue

            if ref_date is not None:
                ref = pd.to_datetime(ref_date)
                closest_idx = (time_sat - ref).abs().argmin()
                baseline = values[closest_idx]
            else:
                baseline = values[0]
            values_aligned = values - baseline
            ax.plot(time_sat, values_aligned, linewidth=1.4, color=colors[k])
            legend_labels.append(list_satellites[k])

            # ── Bandeau RMS/incertitude ──────────────────────────────────
            if rms_variables is not None and rms_variables[k] is not None:
                Satellite_rms = interpolation_altimetry_AWS.satellite_on_aws(
                    info['file'], list_satellites[k], rms_variables[k]
                )
                rms_df = Satellite_rms[['time_sat', rms_variables[k]]].dropna()
                rms_df['time_sat'] = pd.to_datetime(rms_df['time_sat'])
                rms_df[rms_variables[k]] = pd.to_numeric(rms_df[rms_variables[k]], errors='coerce')

                # alignement avec les valeurs dh déjà filtrées/recalées
                main_df = pd.DataFrame({'time_sat': time_sat.values, 'value': values_aligned})
                main_df['value'] = pd.to_numeric(main_df['value'], errors='coerce')
                merged = pd.merge(main_df, rms_df, on='time_sat', how='inner')
                merged = merged.dropna(subset=['value', rms_variables[k]])

                if not merged.empty:
                    ax.fill_between(
                        merged['time_sat'],
                        (merged['value'] - merged[rms_variables[k]]).astype(float),
                        (merged['value'] + merged[rms_variables[k]]).astype(float),
                        alpha=rms_alpha, color=colors[k], linewidth=0,
                    )
        # ── AWS (z_surf_combined_anomaly / gps_alt_anomaly) ────────────────
        if aws_variables is not None:
            aws_funcs = {
                "z_surf_combined_anomaly": z_surf_combined_anomaly,
                "gps_alt_anomaly": gps_alt_anomaly,
            }
            for j, aws_var in enumerate(aws_variables):
                aws_df = aws_funcs[aws_var](info['file'])
                aws_series = aws_df[aws_var].dropna()
                time_aws = aws_df.loc[aws_series.index, "time"]
                values_aws = aws_series.values

                # filtrage temporel
                mask_aws = pd.Series(True, index=time_aws.index)
                if start_bound is not None:
                    mask_aws &= (time_aws >= start_bound)
                if end_bound is not None:
                    mask_aws &= (time_aws <= end_bound)
                time_aws = time_aws[mask_aws]
                values_aws = values_aws[mask_aws.values]

                if len(values_aws) == 0:
                    continue

                if ref_date is not None:
                    ref = pd.to_datetime(ref_date)
                    closest_idx = (time_aws - ref).abs().argmin()
                    baseline_aws = values_aws[time_aws.index.get_loc(closest_idx)] \
                        if closest_idx in time_aws.index else values_aws[0]
                else:
                    baseline_aws = values_aws[0]
                values_aws_aligned = values_aws - baseline_aws
                color_aws = aws_colors[j] if aws_colors is not None else None
                ax.plot(time_aws, values_aws_aligned, linewidth=1.4,
                         linestyle="--", color=color_aws)
                legend_labels.append(aws_var)

        ax.legend(legend_labels, fontsize=7)
        ax.set_title(station, fontsize=12, fontweight="bold")
        ax.set_xlabel("Years", fontsize=9)
        ax.set_ylabel('ice surface elevation changes (m)', fontsize=9)
        ax.tick_params(axis="both", labelsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_major_locator(mdates.YearLocator(5))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    for ax in axes_flat[N:]:
        ax.set_visible(False)
    fig.suptitle(fig_name, fontsize=15, fontweight="bold", y=1.01)
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
    print(f"Figure sauvegardée : {OUTPUT_FILE}")
    plt.show()


# plots_variables(
#     AWS_data.STATIONS_accumulation,
#     'Surface Elevation CHange and RMSE associated for Nilsson et al, accumulation',
#     ['Nilsson'],
#     ['dh'],
#     ["#0D09FFFF"],
#     aws_variables=None,
#     aws_colors=None,
#     ref_date=None,
#     rms_variables=['rms'],   # bandeau EC ± uncertainty
# )

#_________________________________________________________________________________________________________________________
 
def plot_variance_change_test(STATIONS, fig_name, satellite, variable,
                               cutoff_date="2010-11-01",
                               lowess_frac=0.15,
                               color_before="#1f77b4", color_after="#d62728"):
 
    """Pour chaque station, teste si la variance du bruit résiduel de
    `variable` change significativement avant/après `cutoff_date`, et
    visualise le résultat en 2 panneaux :
    - Haut : série brute + tendance LOWESS + ligne verticale de coupure
    - Bas  : résidus (brut - tendance) colorés par période, avec bandes
      ±1 écart-type par période, et annotation du test de Levene
      (statistique, p-value, ratio des variances).
 
    STATIONS : dict{dict}, chaque entrée doit contenir 'file'
    fig_name : str, titre de la figure
    satellite : str, nom du satellite à passer à satellite_on_aws
    variable : str, nom de la variable d'élévation (ex: 'dh')
    cutoff_date : str, date de coupure testée (par défaut '2011-11-01')
    lowess_frac : float, fraction de points utilisée pour le lissage LOWESS
        (0.15 = 15% des points dans chaque fenêtre locale ; à ajuster selon
        la densité temporelle de vos données)
    output_filename : str ou None, nom de fichier court pour la sauvegarde
    """
 
    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
 
    cutoff = pd.to_datetime(cutoff_date)
    N = len(STATIONS)
    NCOLS = 2
    NROWS_STATIONS = (N + NCOLS - 1) // NCOLS
 
    # 2 lignes de sous-graphiques (haut/bas) par station -> grille NROWS_STATIONS*2 x NCOLS
    fig, axes = plt.subplots(
        NROWS_STATIONS * 2, NCOLS,
        figsize=(14, NROWS_STATIONS * 5.5),
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2, 1.3] * NROWS_STATIONS},
    )
    fig.get_layout_engine().set(h_pad=0.4, hspace=0.15, w_pad=0.3, wspace=0.15)
 
    results_summary = []
 
    for i, (station, info) in enumerate(STATIONS.items()):
        row_pair = (i // NCOLS) * 2
        col = i % NCOLS
        ax_top = axes[row_pair, col]
        ax_bot = axes[row_pair + 1, col]
 
        Satellite = interpolation_altimetry_AWS.satellite_on_aws(info['file'], satellite, variable)
        series = Satellite[variable].dropna()
        series1 = Satellite.iloc[series.index.tolist()]
        time_sat = pd.to_datetime(series1['time_sat'])
        values = series.values
 
        order = np.argsort(time_sat.values)             #Ressort les indices qui permettrait de trier la liste par ordre croissant
        time_sat = time_sat.iloc[order]                 #Permet d'avoir les dates par ordre croissant
        values = values[order]                          #Les values sont par odre croissant (values est le dataArray de la variable souhaitée)
 
        # ── Tendance LOWESS et résidus ──────────────────────────────────
        t_numeric = (time_sat - time_sat.iloc[0]).dt.days.values.astype(float)
        smoothed = lowess(values, t_numeric, frac=lowess_frac, return_sorted=False)
        residuals = values - smoothed
 
        mask_before = (time_sat < cutoff).values
        mask_after = (time_sat >= cutoff).values
 
        res_before = residuals[mask_before]
        res_after = residuals[mask_after]
 
        # ── Test de Levene (robuste à la non-normalité) ─────────────────
        if len(res_before) > 1 and len(res_after) > 1:
            stat, pval = stats.levene(res_before, res_after)
            std_before = np.std(res_before)
            std_after = np.std(res_after)
            var_ratio = (std_before**2) / (std_after**2) if std_after > 0 else np.nan
        else:
            stat, pval, std_before, std_after, var_ratio = (np.nan,) * 5
 
        results_summary.append({
            "station": station, "levene_stat": stat, "p_value": pval,
            "std_before": std_before, "std_after": std_after,
            "variance_ratio": var_ratio,
            "n_before": len(res_before), "n_after": len(res_after),
        })
 
        # ── Panneau du haut : série brute + tendance + coupure ──────────
        ax_top.plot(time_sat, values, linewidth=0.8, color="grey", alpha=0.6 , label='Nilsson and Gardner, 2026')
        ax_top.legend()
        ax_top.plot(time_sat, smoothed, linewidth=1.6, color="black", label="Tendance LOWESS")
        ax_top.axvline(cutoff, color="black", linestyle="--", linewidth=1)
        ax_top.set_title(station, fontsize=12, fontweight="bold")
        ax_top.set_ylabel(variable, fontsize=8)
        ax_top.tick_params(axis="both", labelsize=7)
        ax_top.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax_top.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)
 
        # ── Panneau du bas : résidus colorés par période + bandes std ───
        ax_bot.scatter(time_sat[mask_before], res_before, s=6, color=color_before,
                         alpha=0.5, label=f"Before {cutoff_date}")
        ax_bot.scatter(time_sat[mask_after], res_after, s=6, color=color_after,
                         alpha=0.5, label=f"After {cutoff_date}")
        ax_bot.axhline(0, color="black", linewidth=0.6)
        ax_bot.axvline(cutoff, color="black", linestyle="--", linewidth=1)
 
        if len(res_before) > 1:
            ax_bot.axhspan(-std_before, std_before, xmin=0,
                            xmax=(cutoff - time_sat.iloc[0]) / (time_sat.iloc[-1] - time_sat.iloc[0]),
                            color=color_before, alpha=0.12)
        if len(res_after) > 1:
            ax_bot.axhspan(-std_after, std_after,
                            xmin=(cutoff - time_sat.iloc[0]) / (time_sat.iloc[-1] - time_sat.iloc[0]),
                            xmax=1,
                            color=color_after, alpha=0.12)
 
        annotation = (
            f"Levene W = {stat:.2f}\n"
            f"p = {pval:.2e}\n"
            f"$\\sigma$ before = {std_before:.2f}\n"
            f"$\\sigma$ after = {std_after:.2f}\n"
            f"ratio var. = {var_ratio:.1f}x"
        ) if not np.isnan(stat) else "Données insuffisantes"
 
        ax_bot.text(0.02, 0.97, annotation, transform=ax_bot.transAxes,
                     fontsize=7, va="top", ha="left",
                     bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="grey"))
 
        ax_bot.set_xlabel("Years", fontsize=8)
        ax_bot.set_ylabel("Residues (m)", fontsize=8)
        ax_bot.tick_params(axis="both", labelsize=7)
        ax_bot.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax_bot.legend(fontsize=6, loc="lower left")
        ax_bot.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)
 
    # masquer les sous-graphiques inutilisés
    total_slots = NROWS_STATIONS * NCOLS
    for j in range(N, total_slots):
        row_pair = (j // NCOLS) * 2
        col = j % NCOLS
        axes[row_pair, col].set_visible(False)
        axes[row_pair + 1, col].set_visible(False)
 
    fig.suptitle(fig_name, fontsize=15, fontweight="bold")
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
    print(f"Figure sauvegardée : {OUTPUT_FILE}")
    plt.show()
 
    summary_df = pd.DataFrame(results_summary)
    print("\nRésumé du test de Levene par station :")
    print(summary_df.to_string(index=False))
    return summary_df


summary = plot_variance_change_test(
    AWS_data.STATIONS_accumulation,
    "Variance change test for Nilsson and Gardner, 2026 dataset - before and after 2010-11-01",
    satellite="Nilsson and Gardner, 2026",
    variable="dh",
    cutoff_date="2010-11-01",
)
 
def plot_trend_comparison(STATIONS, fig_name, list_satellites, variables, colors,
                           lowess_frac=0.15, output_filename=None,
                           start_date=None, end_date=None):

    """Pour chaque station, trace la série brute (fine, transparente) et sa
    tendance LOWESS (trait plein) pour un ou plusieurs satellites, afin de
    comparer visuellement leurs tendances respectives.

    STATIONS : dict{dict}, chaque entrée doit contenir 'file'
    fig_name : str, titre de la figure
    list_satellites : list[str], noms des satellites à passer à satellite_on_aws
    variables : list[str], même longueur que list_satellites, variable à tracer
        pour chaque satellite
    colors : list[str], même longueur que list_satellites
    lowess_frac : float, fraction de points utilisée pour le lissage LOWESS
        (0.15 = 15% des points dans chaque fenêtre locale)
    output_filename : str ou None, nom de fichier court pour la sauvegarde
    start_date / end_date : bornes temporelles du plot (str ou None)
    """

    if output_filename is None:
        safe_name = "".join(c for c in fig_name if c.isalnum() or c in " _-")
        safe_name = "_".join(safe_name.split())[:60]
    else:
        safe_name = output_filename
    OUTPUT_FILE = safe_name + ".png"

    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()

    N = len(STATIONS)
    NCOLS = 4
    NROWS = (N + NCOLS - 1) // NCOLS

    fig, axes = plt.subplots(
        NROWS, NCOLS,
        figsize=(14, NROWS * 3.5),
        constrained_layout=True,
    )
    axes_flat = axes.flatten()

    for ax, (station, info) in zip(axes_flat, STATIONS.items()):
        legend_labels = []

        for k in range(len(list_satellites)):
            Satellite = interpolation_altimetry_AWS.satellite_on_aws(info['file'], list_satellites[k], variables[k])
            series = Satellite[variables[k]].dropna()
            series1 = Satellite.iloc[series.index.tolist()]
            time_sat = pd.to_datetime(series1['time_sat'])
            values = series.values

            # tri chronologique
            order = np.argsort(time_sat.values)
            time_sat = time_sat.iloc[order]
            values = values[order]

            # filtrage temporel
            mask = pd.Series(True, index=range(len(time_sat)))
            if start_bound is not None:
                mask &= (time_sat >= start_bound).values
            if end_bound is not None:
                mask &= (time_sat <= end_bound).values
            time_sat = time_sat[mask.values]
            values = values[mask.values]

            if len(values) < 3:
                continue

            # ── Tendance LOWESS ──────────────────────────────────────
            t_numeric = (time_sat - time_sat.iloc[0]).dt.days.values.astype(float)
            smoothed = lowess(values, t_numeric, frac=lowess_frac, return_sorted=False)

            ax.plot(time_sat, values, linewidth=0.7, color=colors[k], alpha=0.35)
            ax.plot(time_sat, smoothed, linewidth=1.8, color=colors[k])
            legend_labels.append(list_satellites[k])

        ax.legend(legend_labels, fontsize=7)
        ax.set_title(station, fontsize=12, fontweight="bold")
        ax.set_xlabel("Years", fontsize=9)
        ax.set_ylabel('ice surface elevation changes (m)', fontsize=9)
        ax.tick_params(axis="both", labelsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_major_locator(mdates.YearLocator(5))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    for ax in axes_flat[N:]:
        ax.set_visible(False)

    fig.suptitle(fig_name, fontsize=15, fontweight="bold")
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
    print(f"Figure sauvegardée : {OUTPUT_FILE}")
    plt.show()

# plot_trend_comparison(
#     AWS_data.STATIONS_accumulation,
#     'Comparaison des tendances Nilsson vs Zhang',
#     ['Nilsson', 'Zhang'],
#     ['dh', 'elev_interp'],
#     ["#1f77b4", "#d62728"],
#     lowess_frac=0.15,
# )