from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import AWS_data
import interpolation_altimetry_AWS
from plot_trend_seasonality_residues import decompose_trend_seasonal
from paths import FIGURES_DIR as _BASE_FIGURES_DIR


# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Detrending ───────────────────────────────────────────────────────────
"""
def _apply_detrend_inplace(df, time_col, value_col, poly_degree, n_harmonics) -- removes the polynomial trend from a DataFrame column via decompose_trend_seasonal (keeping the seasonal component) and modifies df in place.
"""


# ── Loading data ──────────────────────────────────────────────────────────
"""
def _load_satellite_df(info, satellite_name, variable) -- loads the satellite DataFrame already matched to the AWS measurements and identifies its time column.

def _load_aws_daily(info, aws_variable) -- loads the full daily AWS series and computes its anomaly using a single baseline (the series' first non-NaN value).
"""


# ── Pairing / alignment ───────────────────────────────────────────────────
"""
def _pair_two_satellites(df_x, time_col_x, var_x, df_y, time_col_y, var_y, freq) -- pairs two satellite series by averaging over a common frequency (monthly by default), keeping only the periods where both have a value.

def _pair_aws_satellite(aws_daily_df, aws_variable, df_sat, time_col_sat, var_sat, freq) -- pairs the full daily AWS series to a satellite series by averaging over a common frequency, following the same logic as _pair_two_satellites.
"""


# ── Plotting ──────────────────────────────────────────────────────────────
"""
def _scatter_with_stats(ax, x, y, color, point_size, point_alpha, show_1to1_line, stat_fontsize, xlabel, ylabel) -- plots an x/y scatter with an optional 1:1 line and the r/p (Pearson) statistics shown in the legend.

def plot_station_correlation_matrix(STATIONS, fig_name, satellite_name, variable, satellite_labels, aws_variable, column_colors, time_agg_freq, detrend, detrend_poly_degree, detrend_n_harmonics, point_size, point_alpha, show_1to1_line, stat_fontsize, subplot_size, save, output_dir) -- plots, for each station, a grid of satellite-satellite and satellite-AWS correlations with scatter plots and r/p statistics.
"""


FIGURES_DIR = _BASE_FIGURES_DIR
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Colonne AWS "brute" associée à chaque variable d'anomalie -- comme dans
# tes fonctions locales gps_alt_anomaly / z_surf_combined_anomaly.
AWS_BASE_COLUMN = {
    "z_surf_combined_anomaly": "z_surf_combined",
    "gps_alt_anomaly": "gps_alt",
}

# Couleurs par défaut, une par colonne (5 satellites ancre + 1 colonne AWS).
# Modifie ce dict pour changer une couleur globalement, ou passe
# `column_colors=` à l'appel pour ne changer que pour un appel donné.
DEFAULT_COLUMN_COLORS = {
    "Copernicus_Climate_Data_Store": "#1f77b4",
    "Nilsson and Gardner, 2026": "#ff7f0e",
    "Andersen et al., 2025": "#2ca02c",
    "Khan et al., 2025": "#d62728",
    "Zhang et al., 2022": "#9467bd",
    "AWS": "black",
}


def _apply_detrend_inplace(df, time_col, value_col, poly_degree, n_harmonics):
    """
    Retire la tendance polynomiale de df[value_col] (la saisonnalité est
    conservée), via decompose_trend_seasonal (plot_trend_seasonality_residues.py) -- même
    principe que ta fonction detrented_time_series.
    Modifie df en place et le renvoie. Les positions où la tendance n'a pas
    pu être ajustée (pas assez de points valides) deviennent NaN, et seront
    naturellement écartées par les dropna() en aval.
    """
    trend, seasonal, residual = decompose_trend_seasonal(
        df[time_col], df[value_col].values,
        poly_degree=poly_degree, n_harmonics=n_harmonics,
    )
    df[value_col] = df[value_col].values - trend
    return df


def _load_satellite_df(info, satellite_name, variable):
    """Charge le DataFrame satellite (déjà apparié aux mesures AWS par
    interpolation_altimetry_AWS), et repère la colonne de temps satellite utilisée."""
    df = interpolation_altimetry_AWS.satellite_on_aws(info["file"], satellite_name, variable)
    time_col = "time_sat" if "time_sat" in df.columns else "time"
    return df, time_col


def _load_aws_daily(info, aws_variable):
    """
    Charge la série AWS QUOTIDIENNE COMPLÈTE (pas seulement les dates où un
    satellite est passé) directement depuis le csv AWS (info['file']), et
    calcule l'anomalie `aws_variable` avec UNE SEULE baseline (1re valeur
    non-NaN de toute la série quotidienne -- contrairement à avant où la
    baseline était recalculée sur le sous-ensemble apparié à chaque
    satellite, ce qui limitait aussi la quantité de données disponibles).
    Renvoie un DataFrame ["time", aws_variable] (sans NaN), ou None.
    """
    base_col = AWS_BASE_COLUMN.get(aws_variable)
    if base_col is None:
        raise ValueError(f"aws_variable={aws_variable!r} inconnu (voir AWS_BASE_COLUMN).")

    df = pd.read_csv(info["file"], parse_dates=["time"])
    if base_col not in df.columns:
        return None

    valid = df[base_col].dropna()
    if valid.empty:
        return None

    df[aws_variable] = df[base_col] - valid.iloc[0]
    return df[["time", aws_variable]].dropna()


def _pair_two_satellites(df_x, time_col_x, var_x, df_y, time_col_y, var_y, freq="M"):
    """
    Apparie deux séries satellite différentes en les regroupant par mois
    (moyenne mensuelle de chaque satellite), puis en ne gardant que les mois
    où LES DEUX satellites ont au moins une valeur -- plutôt qu'une
    tolérance en jours, qui perdait trop de points quand les fréquences de
    survol diffèrent fortement entre satellites.

    freq : fréquence de regroupement pandas ("M" = mensuel, "W" = hebdo,
        "Q" = trimestriel...).
    Renvoie (x_values, y_values) alignés (une valeur par mois en commun),
    ou (None, None) si aucun mois en commun.
    """
    dx = df_x[[time_col_x, var_x]].dropna().copy()
    dy = df_y[[time_col_y, var_y]].dropna().copy()
    if dx.empty or dy.empty:
        return None, None

    dx["period"] = pd.to_datetime(dx[time_col_x]).dt.to_period(freq)
    dy["period"] = pd.to_datetime(dy[time_col_y]).dt.to_period(freq)

    # moyenne par mois (par satellite), au cas où plusieurs passages le même mois.
    # renommées en "x"/"y" (uniques) au cas où var_x == var_y (ex. deux
    # satellites qui utilisent tous les deux "dh") -- sinon pd.concat créerait
    # deux colonnes de même nom et merged[var_x] renverrait un DataFrame à 2
    # colonnes au lieu d'une Series.
    dx_agg = dx.groupby("period")[var_x].mean().rename("x")
    dy_agg = dy.groupby("period")[var_y].mean().rename("y")

    merged = pd.concat([dx_agg, dy_agg], axis=1, join="inner").dropna()
    if merged.empty:
        return None, None
    return merged["x"].values, merged["y"].values


def _pair_aws_satellite(aws_daily_df, aws_variable, df_sat, time_col_sat, var_sat, freq="M"):
    """
    Apparie la série AWS quotidienne complète (aws_daily_df) à une série
    satellite, par moyenne mensuelle des deux séries puis jointure sur les
    mois en commun -- même logique que _pair_two_satellites, mais côté AWS
    on part des données quotidiennes (beaucoup plus denses que le
    sous-ensemble apparié aux dates de passage satellite).
    """
    if aws_daily_df is None or aws_daily_df.empty:
        return None, None

    dx = aws_daily_df.copy()
    dy = df_sat[[time_col_sat, var_sat]].dropna().copy()
    if dx.empty or dy.empty:
        return None, None

    dx["period"] = pd.to_datetime(dx["time"]).dt.to_period(freq)
    dy["period"] = pd.to_datetime(dy[time_col_sat]).dt.to_period(freq)

    dx_agg = dx.groupby("period")[aws_variable].mean().rename("x")
    dy_agg = dy.groupby("period")[var_sat].mean().rename("y")

    merged = pd.concat([dx_agg, dy_agg], axis=1, join="inner").dropna()
    if merged.empty:
        return None, None
    return merged["x"].values, merged["y"].values


def _scatter_with_stats(ax, x, y, color, point_size, point_alpha,
                          show_1to1_line, stat_fontsize, xlabel, ylabel):
    """Trace un scatter x/y avec droite 1:1 optionnelle et r/p en légende."""
    if x is None or y is None or len(x) < 2:
        ax.text(0.5, 0.5, "Pas assez\nde données", ha="center", va="center",
                fontsize=9, style="italic", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        return

    r, p = pearsonr(x, y)
    ax.scatter(x, y, color=color, s=point_size, alpha=point_alpha,
               label=f"r = {r:.2f}\np = {p:.2g}\nn = {len(x)}")

    if show_1to1_line:
        lims = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(lims, lims, linestyle="--", linewidth=1, color="grey", zorder=1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

    ax.legend(fontsize=stat_fontsize, loc="best", framealpha=0.7, handlelength=0)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(True, linestyle="--", alpha=0.4)


def plot_station_correlation_matrix(
        STATIONS, fig_name,
        satellite_name, variable, satellite_labels=None,
        aws_variable="z_surf_combined_anomaly",
        column_colors=None,            # liste [c0..c4, c_aws] ou dict {label: couleur} ; None -> DEFAULT_COLUMN_COLORS
        time_agg_freq="M",             # fréquence de regroupement pour apparier 2 satellites ("M"=mensuel, "W"=hebdo...)
        detrend=False,                  # True -> retire la tendance polynomiale de chaque série (satellites + AWS) avant corrélation
        detrend_poly_degree=1,
        detrend_n_harmonics=1,
        point_size=12, point_alpha=0.7,
        show_1to1_line=True,
        stat_fontsize=7,
        subplot_size=(2.6, 2.6),
        save=True,
        output_dir=None,
        ):
    """
    Pour chaque station de STATIONS, trace une grille de corrélations :
    - colonnes 0..n-1 : chaque satellite comme "ancre" (axe x), corrélé aux
      n-1 AUTRES satellites (axe y, un par ligne) ; la case diagonale
      (satellite vs lui-même) est laissée vide.
    - colonne n (dernière) : AWS (`aws_variable`, axe x) corrélé à chacun
      des n satellites (axe y, un par ligne).
    Chaque case affiche un nuage de points, la droite 1:1, et r/p (Pearson)
    en légende.

    satellite_name, variable : listes de même longueur n (ordre = ordre des
        colonnes ET des lignes pour la partie satellite-satellite).
    satellite_labels : labels d'affichage ; None -> reprend satellite_name.
    column_colors : couleur par colonne. Liste de longueur n+1 (dans l'ordre
        satellite_name + AWS en dernier), ou dict {label: couleur} (les clés
        non trouvées retombent sur DEFAULT_COLUMN_COLORS). None -> défauts
        du module.
    time_agg_freq : fréquence de regroupement pandas ("M" par défaut =
        mensuel) utilisée pour apparier deux satellites différents : chaque
        satellite est moyenné par période, puis on ne garde que les
        périodes où les deux ont une valeur.
    detrend : si True, retire la tendance polynomiale (garde la
        saisonnalité) de chaque série -- satellites ET AWS -- avant de
        calculer les corrélations, via decompose_trend_seasonal
        (plot_trend_seasonality_residues.py). Si False (défaut), les séries brutes sont utilisées.
    detrend_poly_degree, detrend_n_harmonics : paramètres transmis à
        decompose_trend_seasonal si detrend=True (voir plot_trend_seasonality_residues.py).
    """
    n = len(satellite_name)
    if len(variable) != n:
        raise ValueError("satellite_name et variable doivent avoir la même longueur.")
    if satellite_labels is None:
        satellite_labels = list(satellite_name)

    # ── résolution des couleurs de colonne ──────────────────────────────
    resolved_colors = dict(DEFAULT_COLUMN_COLORS)
    if isinstance(column_colors, dict):
        resolved_colors.update(column_colors)
    elif column_colors is not None:
        if len(column_colors) != n + 1:
            raise ValueError("column_colors (liste) doit avoir n+1 éléments (n satellites + AWS).")
        keys = list(satellite_name) + ["AWS"]
        resolved_colors.update(dict(zip(keys, column_colors)))

    col_colors = [resolved_colors.get(sat, "tab:blue") for sat in satellite_name]
    aws_color = resolved_colors.get("AWS", "black")

    out_dir = Path(output_dir) if output_dir is not None else FIGURES_DIR

    for station, info in STATIONS.items():

        # ── chargement : un DataFrame par satellite, + AWS quotidien à part ──
        sat_dfs, sat_time_cols = [], []
        for sat_name, var_name in zip(satellite_name, variable):
            df, time_col = _load_satellite_df(info, sat_name, var_name)
            if detrend:
                df = _apply_detrend_inplace(
                    df, time_col, var_name, detrend_poly_degree, detrend_n_harmonics
                )
            sat_dfs.append(df)
            sat_time_cols.append(time_col)

        aws_daily_df = _load_aws_daily(info, aws_variable)
        if detrend and aws_daily_df is not None:
            aws_daily_df = _apply_detrend_inplace(
                aws_daily_df, "time", aws_variable, detrend_poly_degree, detrend_n_harmonics
            )

        n_cols = n + 1
        n_rows = n
        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(subplot_size[0] * n_cols, subplot_size[1] * n_rows),
            constrained_layout=True,
        )
        axes = np.atleast_2d(axes)

        # ── colonnes satellite-satellite (diagonale vide) ────────────────
        for j in range(n):  # colonne = satellite ancre
            for i in range(n):  # ligne = autre satellite
                ax = axes[i, j]
                if i == j:
                    ax.axis("off")
                    continue

                x, y = _pair_two_satellites(
                    sat_dfs[j], sat_time_cols[j], variable[j],
                    sat_dfs[i], sat_time_cols[i], variable[i],
                    freq=time_agg_freq,
                )
                _scatter_with_stats(
                    ax, x, y, col_colors[j], point_size, point_alpha,
                    show_1to1_line, stat_fontsize,
                    xlabel=satellite_labels[j], ylabel=satellite_labels[i],
                )

        # ── dernière colonne : AWS vs chaque satellite ────────────────────
        for i in range(n):
            ax = axes[i, n]
            x, y = _pair_aws_satellite(
                aws_daily_df, aws_variable,
                sat_dfs[i], sat_time_cols[i], variable[i],
                freq=time_agg_freq,
            )
            _scatter_with_stats(
                ax, x, y, aws_color, point_size, point_alpha,
                show_1to1_line, stat_fontsize,
                xlabel=f"AWS {aws_variable}", ylabel=satellite_labels[i],
            )

        # ── en-têtes de colonnes / lignes ─────────────────────────────────
        for j in range(n):
            axes[0, j].set_title(satellite_labels[j], fontsize=10, fontweight="bold")
        axes[0, n].set_title("AWS", fontsize=10, fontweight="bold")

        fig.suptitle(
            f"{fig_name} — {station}" + (" (détrendé)" if detrend else ""),
            fontsize=14, fontweight="bold",
        )

        if save:
            out_path = out_dir / f"{fig_name}_{station}.png"
            fig.savefig(out_path, dpi=200, bbox_inches="tight")
            print(f"Figure sauvegardée : {out_path}")

        plt.show()


# ── Exemple d'appel ──────────────────────────────────────────────────────


plot_station_correlation_matrix(
    AWS_data.STATION_UPE_L,
    "Correlation matrix",
    satellite_name=['Copernicus_Climate_Data_Store', 'Nilsson and Gardner, 2026',
                    'Andersen et al., 2025', 'Khan et al., 2025', 'Zhang et al., 2022'],
    variable=['dh', 'dh', 'ZZ', 'dh_vol', 'elev_interp'],
    aws_variable="z_surf_combined_anomaly",
    time_agg_freq="M",
    detrend=False,             # ou False pour garder les séries brutes
    detrend_poly_degree=2,
    detrend_n_harmonics=1,
)