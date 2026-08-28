import os
from pathlib import Path

import AWS_data
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plot_time_series_altimetry_AWS
import interpolation_altimetry_AWS
import scipy.optimize
import sklearn
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.tsa.seasonal import STL, seasonal_decompose
from statsmodels.tsa.stattools import adfuller

from paths import FIGURES_DIR as _BASE_FIGURES_DIR

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Trend/seasonal decomposition ────────────────────────────────────────────
"""
def decompose_trend_seasonal(time, values, poly_degree, n_harmonics) -- simultaneously fits a polynomial trend and a harmonic (sin/cos) seasonality to a time series by least squares, returning (trend, seasonal, residual).

def fit_sin(tt, yy, t0) -- fits a single sinusoid to a time series and returns its amplitude/frequency/phase/offset plus the date of its yearly maximum.
"""

# ── Small metrics helpers ───────────────────────────────────────────────────
"""
def mae(y_true, y_pred) -- mean absolute error (wraps sklearn.metrics.mean_absolute_error).

def autocorrelation(x) -- home-made ACF via np.correlate (positive lags only); kept for compatibility/ad hoc use -- for figures, prefer statsmodels.graphics.tsaplots.plot_acf, which provides white-noise confidence bands.
"""

# ── Detrended residuals and AWS anomaly trend ───────────────────────────────
"""
def compute_detrended_residuals(station_file, satellite_name, variable, poly_degree, n_harmonics, ref_date, ...) -- computes the detrended residuals of a satellite series matched to a station, using decompose_trend_seasonal.

def compute_aws_anomaly_trend(station_file, poly_degree, n_harmonics, ref_date, start_date, end_date, ...) -- computes the trend/seasonal/residual decomposition of the AWS's own elevation anomaly (rather than a satellite series).
"""

# ── Per-station trend-removed plots ─────────────────────────────────────────
"""
def plot_variable_trend_removed(station_name, station_file, satellite_name, variable, color, fig_name, ...) -- plots one satellite variable at a station with its polynomial trend removed.

def plot_station_residuals_scatter(station_name, station_file, list_satellites, poly_degree, n_harmonics, ref_date, ...) -- scatter of detrended residuals across the requested satellites for one station.

def subplots_detrending_station(station_name, station_file, list_satellites, poly_degree, n_harmonics, ref_date, ...) -- grid of subplots summarizing the trend/seasonal/residual decomposition of every requested satellite at one station.
"""

# ── Fit summary table across stations ───────────────────────────────────────
"""
def build_fit_summary_table(stations_dict, list_satellites, poly_degree, n_harmonics, ref_date, start_date, ...) -- builds and optionally saves a CSV table summarizing the trend fit (and per-station figures) across every station in stations_dict and every requested satellite.
"""

# ── Paths ─────────────────────────────────────────────────────────────────
"""
FIGURES_DIR_ROOT -- output folder for this module's figures, a subfolder of paths.py's shared FIGURES_DIR.
"""


print(np.__version__)

FIGURES_DIR_ROOT = _BASE_FIGURES_DIR / "detrending"


def decompose_trend_seasonal(time, values, poly_degree=2, n_harmonics=2):
    """Ajuste simultanément une tendance polynomiale et une saisonnalité
    harmonique (sin/cos) sur une série temporelle, par moindres carrés.

    time : array-like de datetime64/Timestamp
    values : array-like de float (peut contenir des NaN, ignorés)
    poly_degree : degré du polynôme de tendance (1=linéaire, 2=quadratique, ...)
    n_harmonics : nombre d'harmoniques saisonnières (1=cycle annuel seul,
        2=annuel + semi-annuel, etc.)

    Retourne (trend, seasonal, residual), trois arrays de même longueur que
    values (avec NaN aux positions où values était NaN).
    """
    time = pd.DatetimeIndex(time)
    values = np.asarray(values, dtype="float64")

    valid = np.isfinite(values)
    if valid.sum() < (poly_degree + 1 + 2 * n_harmonics + 1):
        nan_arr = np.full_like(values, np.nan)
        return nan_arr, nan_arr, nan_arr

    t0 = time[valid][0]
    t_years = (time - t0).total_seconds().values / (365.25 * 24 * 3600)

    cols = []
    for d in range(poly_degree + 1):
        cols.append(t_years**d)

    for k in range(1, n_harmonics + 1):
        cols.append(np.sin(2 * np.pi * k * t_years))
        cols.append(np.cos(2 * np.pi * k * t_years))
    X = np.column_stack(cols)

    coeffs, *_ = np.linalg.lstsq(X[valid], values[valid], rcond=None)

    n_poly = poly_degree + 1
    trend = X[:, :n_poly] @ coeffs[:n_poly]
    seasonal = X[:, n_poly:] @ coeffs[n_poly:]

    trend = np.where(valid, trend, np.nan)
    seasonal = np.where(valid, seasonal, np.nan)
    residual = values - trend - seasonal

    return trend, seasonal, residual


import pandas as pd

def fit_sin(tt, yy, t0=None):
    """
    Fit sin to the input time sequence, and return fitting parameters:
    amp, omega, phase, offset, freq, period, fitfunc, date_max, day_of_year_max

    t0 : datetime-like optionnel, date réelle correspondant à tt=0.
         Nécessaire pour calculer date_max / day_of_year_max quand tt
         est un tableau de floats (ex: années écoulées depuis t0),
         plutôt qu'un tableau datetime64.
    """

    tt = np.array(tt)
    yy = np.array(yy, dtype=float)

    is_datetime = np.issubdtype(tt.dtype, np.datetime64)

    if is_datetime:
        tt_num = (tt - tt[0]) / np.timedelta64(1, "s")
        t0_real = pd.Timestamp(tt[0])
        unit_seconds = 1.0  # tt_num est déjà en secondes
    elif t0 is not None:
        # tt est en années écoulées depuis t0
        t0_real = pd.Timestamp(t0)
        tt_num = tt  # on garde tel quel (unités arbitraires cohérentes avec p0)
        unit_seconds = None
    else:
        try:
            tt_num = np.array([(t - tt[0]).total_seconds() for t in tt], dtype=float)
            t0_real = None
        except Exception:
            tt_num = tt.astype(float)
            t0_real = None
        unit_seconds = None

    dt = tt_num[1] - tt_num[0]

    ff = np.fft.fftfreq(len(tt_num), dt)
    Fyy = np.abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:]) + 1])

    guess_amp = np.std(yy) * np.sqrt(2)
    guess_offset = np.mean(yy)

    def sinfunc(t, A, w, p, c):
        return A * np.sin(w * t + p) + c

    popt, pcov = scipy.optimize.curve_fit(
        sinfunc, tt_num, yy,
        p0=[0.3, 2 * np.pi, -2.5, 0],
        maxfev=10000
    )

    A, w, p, c = popt
    f = w / (2.0 * np.pi)
    T = 1.0 / f if f != 0 else np.inf

    fitfunc = lambda t: A * np.sin(w * t + p) + c

    if A < 0:
        A_label = np.abs(A)
        p_label = p % np.pi
    else:
        A_label = A
        p_label = p % np.pi

    # --- Jour du maximum de la sinusoïde (formule directe via T et φ) ---
    sign_A = 1.0 if A >= 0 else -1.0
    t_max = (T / (2 * np.pi)) * (sign_A * np.pi / 2 - p)
    t_max = t_max % T  # dans la même unité que tt_num

    date_max = None
    day_of_year_max = None

    if is_datetime:
        date_max = t0_real + pd.Timedelta(seconds=t_max)
        day_of_year_max = date_max.dayofyear
    elif t0_real is not None:
        # tt_num est en années -> convertir t_max (années) en jours
        date_max = t0_real + pd.Timedelta(days=t_max * 365.25)
        day_of_year_max = date_max.dayofyear

    return {
        "A_label": A_label,
        "p_label": p_label,
        "A": A,
        "ω": w,
        "φ": p,
        "offset": c,
        "freq": f,
        "T": T,
        "fitfunc": fitfunc,
        "maxcov": np.max(pcov),
        "rawres": (popt, pcov),
        "t_max": t_max,
        "date_max": date_max,
        "day_of_year_max": day_of_year_max,
    }



def mae(y_true, y_pred):
    return sklearn.metrics.mean_absolute_error(y_true, y_pred)


def autocorrelation(x):
    """ACF maison par np.correlate (lags positifs uniquement).
    Conservée pour compatibilité / usage ponctuel ; pour les figures,
    préférer statsmodels.graphics.tsaplots.plot_acf qui fournit les
    bandes de confiance sous H0 de bruit blanc."""
    n = len(x)
    mean = np.mean(x)
    var = np.var(x)
    result = np.correlate(x - mean, x - mean, mode="full") / (var * n)
    return result[result.size // 2:]


def compute_detrended_residuals(
    station_file,
    satellite_name,
    variable,
    poly_degree=2,
    n_harmonics=1,
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
):
    """Calcule tendance, signal détrendé, fit sinusoïdal et résidus pour un
    couple station/satellite, sans rien tracer.

    Retourne un dict avec les clés :
        time_sat, values_aligned, trend, valid_trend,
        time (restreint à valid_trend), detrended, season, residuals,
        fit_result
    ou None si les données sont absentes, la décomposition impossible, ou
    le fit sinusoïdal échoue.
    """
    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()

    Satellite = interpolation_altimetry_AWS.satellite_on_aws(station_file, satellite_name, variable)
    series = Satellite[variable].dropna()
    series1 = Satellite.iloc[series.index.tolist()]
    time_sat = pd.to_datetime(series1["time_sat"])
    values = series.values

    mask = pd.Series(True, index=range(len(time_sat)))
    if start_bound is not None:
        mask &= (time_sat >= start_bound).values
    if end_bound is not None:
        mask &= (time_sat <= end_bound).values
    time_sat = time_sat[mask.values]
    values = values[mask.values]

    if len(values) == 0:
        return None

    if ref_date is not None:
        ref = pd.to_datetime(ref_date)
        closest_idx = (time_sat - ref).abs().argmin()
        baseline = values[closest_idx]
    else:
        baseline = values[0]
    values_aligned = values - baseline

    trend, seasonal, residual = decompose_trend_seasonal(
        time_sat, values_aligned, poly_degree=poly_degree, n_harmonics=n_harmonics
    )
    valid_trend = np.isfinite(trend)
    if not valid_trend.any():
        return None

    detrended = values_aligned - trend  # NaN hors de valid_trend
    time_valid = time_sat[valid_trend]

    t0 = time_valid.iloc[0]
    t_years = ((time_valid - t0).dt.total_seconds() / (365.25 * 24 * 3600)).to_numpy()

    try:
        fit_result = fit_sin(t_years, detrended[valid_trend], t0=t0)
    except Exception as e:
        print(f"Échec de l'ajustement sinusoïdal pour {satellite_name} : {e}")
        return None

    season = fit_result["fitfunc"](t_years)
    residuals = detrended[valid_trend] - season

    return {
        "time_sat": time_sat,
        "values_aligned": values_aligned,
        "trend": trend,
        "valid_trend": valid_trend,
        "time": time_valid,
        "detrended": detrended,
        "season": season,
        "residuals": residuals,
        "fit_result": fit_result,
    }


def compute_aws_anomaly_trend(
    station_file,
    poly_degree=2,
    n_harmonics=1,
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
):
    """Calcule l'anomalie AWS z_surf_combined (via plot_time_series_altimetry_AWS.z_surf_combined_anomaly)
    et sa tendance polynomiale, pour superposition dans le subplot
    'Raw signal + trend'.

    station_file : même chemin que celui passé aux fonctions satellite
        (ex: info_aws["file"] dans AWS_data.STATIONS_ablation).

    Retourne un dict {time, values, trend, valid_trend} ou None si les
    données AWS sont absentes / la décomposition impossible.
    """
    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()

    aws_df = plot_time_series_altimetry_AWS.z_surf_combined_anomaly(station_file)
    aws_series = aws_df["z_surf_combined_anomaly"].dropna()
    time_aws = pd.to_datetime(aws_df.loc[aws_series.index, "time"])
    values_aws = aws_series.values

    mask = pd.Series(True, index=time_aws.index)
    if start_bound is not None:
        mask &= time_aws >= start_bound
    if end_bound is not None:
        mask &= time_aws <= end_bound
    time_aws = time_aws[mask]
    values_aws = values_aws[mask.values]

    if len(values_aws) == 0:
        return None

    if ref_date is not None:
        ref = pd.to_datetime(ref_date)
        closest_idx = np.argmin(np.abs((time_aws - ref).values))
        baseline_aws = values_aws[closest_idx]
    else:
        baseline_aws = values_aws[0]
    values_aws_aligned = values_aws - baseline_aws

    time_aws_dt = pd.DatetimeIndex(time_aws)
    trend, _, _ = decompose_trend_seasonal(
        time_aws_dt, values_aws_aligned, poly_degree=poly_degree, n_harmonics=n_harmonics
    )
    valid_trend = np.isfinite(trend)
    if not valid_trend.any():
        return None

    return {
        "time": time_aws_dt,
        "values": values_aws_aligned,
        "trend": trend,
        "valid_trend": valid_trend,
    }


def plot_variable_trend_removed(
    station_name,
    station_file,
    satellite_name,
    variable,
    color,
    fig_name=None,
    poly_degree=2,
    n_harmonics=1,
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
    sat_band_alpha=0.15,
    axes=None,
    save=True,
    figures_dir=None,
    plot_aws=True,
    aws_color="grey",
):
    """Produit 3 subplots (raw+trend / détrendé+fit / résidus) pour un
    couple station/satellite. Si axes est fourni (3 axes existants), dessine
    dedans sans créer de figure ni sauvegarder — utile pour composer une
    grille.

    plot_aws : bool, si True superpose dans le subplot 1 l'anomalie AWS
        z_surf_combined (via compute_aws_anomaly_trend) et sa propre
        tendance polynomiale.
    aws_color : couleur utilisée pour la courbe AWS et sa tendance.
    """
    standalone = axes is None
    if standalone:
        fig, axes = plt.subplots(
            3, 1, figsize=(12, 5), gridspec_kw={"hspace": 0.4}, sharex=True
        )
    else:
        fig = axes[0].figure

    result = compute_detrended_residuals(
        station_file,
        satellite_name,
        variable,
        poly_degree=poly_degree,
        n_harmonics=n_harmonics,
        ref_date=ref_date,
        start_date=start_date,
        end_date=end_date,
    )
    if result is None:
        print(f"Pas de données / fit impossible pour {station_name}/{satellite_name}")
        return None

    time_sat = result["time_sat"]
    values_aligned = result["values_aligned"]
    trend = result["trend"]
    valid_trend = result["valid_trend"]
    time_valid = result["time"]
    detrended = result["detrended"]
    season = result["season"]
    residuals = result["residuals"]
    fit_result = result["fit_result"]

    # --- Subplot 1 : signal brut + tendance (satellite) ---
    axes[0].plot(
        time_sat, values_aligned, linewidth=1.2, color=color,
        label=satellite_name, zorder=3,
    )
    axes[0].plot(
        time_sat[valid_trend], trend[valid_trend], linewidth=1.2,
        color="black", zorder=3,
    )

    # --- AWS z_surf_combined anomaly + sa propre tendance ---
    if plot_aws:
        aws_result = compute_aws_anomaly_trend(
            station_file,
            poly_degree=poly_degree,
            n_harmonics=n_harmonics,
            ref_date=ref_date,
            start_date=start_date,
            end_date=end_date,
        )
        if aws_result is not None:
            axes[0].plot(
                aws_result["time"], aws_result["values"], linewidth=1.0,
                color=aws_color, alpha=0.8, label="AWS z_surf_combined",
                zorder=2,
            )
            axes[0].plot(
                aws_result["time"][aws_result["valid_trend"]],
                aws_result["trend"][aws_result["valid_trend"]],
                linewidth=1.2, color=aws_color, linestyle="--", zorder=2,
            )

    axes[0].set_title("Raw signal + trend", y=1)
    axes[0].set_ylabel("SEC (m)", fontsize=12)
    axes[0].grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    # --- Subplot 2 : détrendé + fit sinusoïdal ---
    axes[1].plot(
        time_valid, detrended[valid_trend], "o", color=color,
        label=satellite_name, zorder=3,
    )
    axes[1].plot(time_valid, season, linewidth=1.2, color="black", zorder=3)
    fit_label = (
        f"A={fit_result['A_label']:.3f} m\n"
        f"T={fit_result['T']:.2f} yr\n"
        f"φ={fit_result['p_label']:.2f}\n"
        f"Day of maximum={fit_result['day_of_year_max']}"
    )
    axes[1].text(
        0.02, 0.95, fit_label,
        transform=axes[1].transAxes, verticalalignment="top", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
    )
    axes[1].set_title("Raw - trend (seasonal + residual)", y=1)
    axes[1].set_ylabel("SEC (m)", fontsize=12)
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axes[1].xaxis.set_major_locator(mdates.YearLocator(1))
    plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=45, ha="right")
    axes[1].grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    # --- Subplot 3 : résidus ---
    mae_residuals = mae(np.zeros(len(residuals)), residuals)
    axes[2].plot(time_valid, residuals, "o", color=color, label=satellite_name, zorder=3)
    axes[2].fill_between(
        time_valid,
        residuals - np.std(residuals),
        residuals + np.std(residuals),
        color=color,
        alpha=sat_band_alpha,
        linewidth=0,
        zorder=1,
    )
    axes[2].axhline(0, color="black", linewidth=0.8, linestyle="--")
    axes[2].set_title("Residuals (detrended - fitted sinusoid)", y=1)
    axes[2].set_ylabel("SEC (m)", fontsize=12)
    axes[2].text(
        0.02, 0.95, f"MAE on residuals: {mae_residuals:.4f} m",
        transform=axes[2].transAxes, verticalalignment="top", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
    )
    axes[2].grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axes[2].xaxis.set_major_locator(mdates.YearLocator(1))
    plt.setp(axes[2].xaxis.get_majorticklabels(), rotation=45, ha="right")

    if standalone:
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(
            handles, labels,
            loc="lower center", ncol=len(labels),
            bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=11,
        )
        fig.subplots_adjust(bottom=0.12)

        out_dir = figures_dir or FIGURES_DIR_ROOT
        out_dir.mkdir(parents=True, exist_ok=True)
        OUTPUT_FILE = f"Time series trend removed {station_name} {satellite_name}.png"
        fig.suptitle(f"{station_name} - {satellite_name}", fontsize=15, fontweight="bold")
        plt.savefig(out_dir / OUTPUT_FILE, dpi=150, bbox_inches="tight")
        print(f"Figure sauvegardée : {OUTPUT_FILE}")

    return fit_result


def plot_station_residuals_scatter(
    station_name,
    station_file,
    list_satellites,
    poly_degree=2,
    n_harmonics=1,
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
    lag=1,
    figsize_per_axes=(4, 4),
    figures_dir=None,
    save=True,
):
    """Une figure par station : scatter résidu(t) vs résidu(t-lag) du fit
    sinusoïdal, une colonne par satellite dans list_satellites.

    Un nuage aligné sur la diagonale y=x indique l'absence d'autocorrélation
    au lag choisi ; un nuage étiré le long d'une droite indique une
    autocorrélation."""
    out_dir = figures_dir or (FIGURES_DIR_ROOT / station_name)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_sats = len(list_satellites)
    fig, axes = plt.subplots(
        1, n_sats,
        figsize=(figsize_per_axes[0] * n_sats, figsize_per_axes[1]),
        squeeze=False,
    )
    axes = axes[0]

    for j, satellite_name in enumerate(list_satellites):
        ax = axes[j]
        info_sat = interpolation_altimetry_AWS.SATELLITE[satellite_name]

        result = compute_detrended_residuals(
            station_file,
            satellite_name,
            info_sat["var"],
            poly_degree=poly_degree,
            n_harmonics=n_harmonics,
            ref_date=ref_date,
            start_date=start_date,
            end_date=end_date,
        )
        if result is None:
            ax.set_visible(False)
            continue

        residuals = result["residuals"]
        if len(residuals) <= lag:
            ax.text(
                0.5, 0.5, "Trop peu de points", ha="center", va="center",
                transform=ax.transAxes,
            )
            continue

        x = residuals[:-lag]
        y = residuals[lag:]
        ax.scatter(x, y, s=14, color=info_sat["color"], alpha=0.6)

        lims = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(lims, lims, color="black", linewidth=0.8, linestyle="--", zorder=1)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        ax.set_title(satellite_name, fontsize=11, fontweight="bold")
        ax.set_xlabel(f"résidu(t-{lag})", fontsize=10)
        if j == 0:
            ax.set_ylabel("résidu(t)", fontsize=10)

    fig.suptitle(
        f"{station_name} — scatter lag-{lag} des résidus", fontsize=15, fontweight="bold"
    )
    fig.tight_layout()

    if save:
        OUTPUT_FILE = out_dir / f"Residuals lag scatter {station_name}.png"
        plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
        print(f"Figure sauvegardée : {OUTPUT_FILE}")

    # plt.show()
    return fig, axes


def subplots_detrending_station(
    station_name,
    station_file,
    list_satellites,
    poly_degree=2,
    n_harmonics=1,
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
):
    """Une figure par station, avec 3 lignes (raw+trend / détrendé+fit /
    résidus) x N colonnes (une par satellite dans list_satellites)."""
    figures_dir = FIGURES_DIR_ROOT / station_name
    figures_dir.mkdir(parents=True, exist_ok=True)

    NROWS = 3
    NCOLS = len(list_satellites)

    fig, axes = plt.subplots(
        NROWS, NCOLS,
        figsize=(NCOLS * 5, NROWS * 3),
        sharex="col",
        sharey="row",
        gridspec_kw={"hspace": 0.5, "wspace": 0.35},
    )
    if NCOLS == 1:
        axes = axes.reshape(NROWS, 1)

    station_fit_results = {}
    legend_handles, legend_labels = [], []
    for k, satellite_name in enumerate(list_satellites):
        info_sat = interpolation_altimetry_AWS.SATELLITE[satellite_name]
        fit_result = plot_variable_trend_removed(
            station_name=station_name,
            station_file=station_file,
            satellite_name=satellite_name,
            variable=info_sat["var"],
            color=info_sat["color"],
            poly_degree=poly_degree,
            n_harmonics=n_harmonics,
            ref_date=ref_date,
            start_date=start_date,
            end_date=end_date,
            axes=axes[:, k],
            save=False,
        )
        station_fit_results[satellite_name] = fit_result
        axes[0, k].set_title(satellite_name, fontsize=11, fontweight="bold")

        handles, labels = axes[0, k].get_legend_handles_labels()
        for h, l in zip(handles, labels):
            if l not in legend_labels:
                legend_handles.append(h)
                legend_labels.append(l)

    amplitudes = [r["A"] for r in station_fit_results.values() if r is not None]
    periods = [r["T"] for r in station_fit_results.values() if r is not None]
    mean_amplitude = np.mean(amplitudes) if amplitudes else np.nan
    mean_period = np.mean(periods) if periods else np.nan

    fig.suptitle(
        f"{station_name}\n"
        f"Mean Amplitude (all datasets) = {mean_amplitude:.3f} m   |   "
        f"Mean Period = {mean_period:.2f} yr",
        fontsize=14,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.85, bottom=0.12)

    fig.legend(
        legend_handles, legend_labels,
        loc="lower center", ncol=len(legend_labels),
        bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=11,
    )

    OUTPUT_FILE = figures_dir / f"Time series decomposition {station_name}.png"
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
    print(f"Figure sauvegardée : {OUTPUT_FILE}")
    # plt.show()


def build_fit_summary_table(
    stations_dict,
    list_satellites,
    poly_degree=2,
    n_harmonics=1,
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
    save=True,
    output_path=None,
):
    """Tableau récapitulatif : pour chaque satellite, moyenne (et écart-type)
    de l'amplitude et de la période du fit sinusoïdal, moyennées sur toutes
    les stations de stations_dict.

    Retourne un DataFrame indexé par satellite avec les colonnes :
        n_stations, mean_amplitude_m, std_amplitude_m,
        mean_period_yr, std_period_yr
    """
    records = {satellite_name: {"A": [], "T": []} for satellite_name in list_satellites}

    for station_name, info_aws in stations_dict.items():
        for satellite_name in list_satellites:
            info_sat = interpolation_altimetry_AWS.SATELLITE[satellite_name]
            result = compute_detrended_residuals(
                info_aws["file"],
                satellite_name,
                info_sat["var"],
                poly_degree=poly_degree,
                n_harmonics=n_harmonics,
                ref_date=ref_date,
                start_date=start_date,
                end_date=end_date,
            )
            if result is None:
                continue
            records[satellite_name]["A"].append(result["fit_result"]["A"])
            records[satellite_name]["T"].append(result["fit_result"]["T"])

    rows = []
    for satellite_name, vals in records.items():
        n = len(vals["A"])
        rows.append(
            {
                "satellite": satellite_name,
                "n_stations": n,
                "mean_amplitude_m": np.mean(vals["A"]) if n else np.nan,
                "std_amplitude_m": np.std(vals["A"]) if n else np.nan,
                "mean_period_yr": np.mean(vals["T"]) if n else np.nan,
                "std_period_yr": np.std(vals["T"]) if n else np.nan,
            }
        )

    summary_df = pd.DataFrame(rows).set_index("satellite")

    if save:
        out_path = output_path or (FIGURES_DIR_ROOT / "fit_summary_by_satellite.csv")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(out_path)
        print(f"Tableau sauvegardé : {out_path}")

    return summary_df


if __name__ == "__main__":
    # --- Figures détaillées par station / satellite ---
    fitted_params = {}
    for station, info_aws in AWS_data.STATIONS_ablation.items():
        station_figures_dir = FIGURES_DIR_ROOT / station
        station_figures_dir.mkdir(parents=True, exist_ok=True)

        fitted_params[station] = {}
        for satellite, info_sat in interpolation_altimetry_AWS.SATELLITE.items():
            fitted_params[station][satellite] = plot_variable_trend_removed(
                station_name=station,
                station_file=info_aws["file"],
                satellite_name=satellite,
                variable=info_sat["var"],
                color=info_sat["color"],
                poly_degree=2,
                n_harmonics=1,
                ref_date=None,
                start_date="2010-11-01",
                end_date=None,
                figures_dir=station_figures_dir,
            )

        # --- Scatter résidu(t) vs résidu(t-lag) pour cette station ---
        # plot_station_residuals_scatter(
        #     station_name=station,
        #     station_file=info_aws["file"],
        #     list_satellites=list(interpolation_altimetry_AWS.SATELLITE.keys()),
        #     poly_degree=2,
        #     n_harmonics=1,
        #     start_date="2010-11-01",
        #     lag=1,
        #     figures_dir=station_figures_dir,
        # )

    # --- Grille de subplots par station (raw+trend / détrendé+fit / résidus) ---
    # for station, info_aws in AWS_data.STATIONS_ablation.items():
    #     subplots_detrending_station
    #         station_name=station,
    #         station_file=info_aws["file"],
    #         list_satellites=list(interpolation_altimetry_AWS.SATELLITE.keys()),
    #         poly_degree=2,
    #         n_harmonics=1,
    #         ref_date=None,
    #         start_date="2010-11-01",
    #         end_date=None,
    #     )

    # --- Scatter résidu(t) vs résidu(t-lag) : maintenant appelé par station
    # directement dans la boucle ci-dessus (voir plot_station_residuals_scatter).

    # --- Tableau récapitulatif amplitude/période moyennes par satellite ---
    # summary_df = build_fit_summary_table(
    #     stations_dict=AWS_data.STATIONS_ablation,
    #     list_satellites=list(interpolation_altimetry_AWS.SATELLITE.keys()),
    #     poly_degree=2,
    #     n_harmonics=1,
    #     start_date="2010-11-01",
    # )
    # print(summary_df)