import os
from pathlib import Path

import AWS_data
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import interpolation_altimetry_AWS
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from paths import PRODEM_DATA_FILE, FIGURES_DIR as _BASE_FIGURES_DIR

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── AWS anomaly / small helpers ──────────────────────────────────────────
"""
def monthly_std_hourly(hourly_path, columns, min_valid_fraction) -- monthly standard deviation of one or more hourly AWS columns, grouped by real calendar year-month.

def z_surf_combined_anomaly(aws_path) -- loads an AWS CSV and adds a column with the combined surface elevation anomaly relative to its first valid value.

def gps_alt_anomaly(aws_path) -- loads an AWS CSV and adds a column with the GPS altitude anomaly relative to its first valid value.

def snow_height(aws_path) -- loads an AWS CSV as-is (passthrough, for callers that just need the raw snow_height column alongside the other AWS variables).
"""

# ── Per-station constants ─────────────────────────────────────────────────
"""
AWS_BASE_COLUMN -- maps each anomaly variable name (e.g. "z_surf_combined_anomaly") to the raw AWS column it's derived from.

DATES -- per-station date window (start/end) used where a station's Delta z_surf_combined or plotting range needs a station-specific default instead of the dataset's full extent.
"""

# ── Multi-station overview plots ─────────────────────────────────────────
"""
def plots_variables(STATIONS, fig_name, list_satellites, variables, colors, ...) -- grid of subplots, one per station, plotting the requested satellite variable(s) and matching AWS variable(s).

def plots_variables_with_satellite_uncertainty(STATIONS, fig_name, list_satellites, variables, colors, ...) -- same as plots_variables, adding each satellite's own reported uncertainty as a band around its curve.

def _multi_satellite_std_band(info, list_satellites, variables, ...) -- computes the inter-product mean and standard deviation band across the requested satellites at one station, used to draw a "spread across datasets as uncertainty" band.

def plots_variables_full_uncertainty(STATIONS, fig_name, list_satellites, variables, colors, ...) -- grid of subplots combining the AWS series, every requested satellite, and (optionally) the multi-satellite mean +/- inter-product std band from _multi_satellite_std_band.
"""

# ── GC-Net precise GPS locations ──────────────────────────────────────────
"""
def load_gcnet_precise_locations(csv_path) -- loads the GC-Net precise GPS measurements CSV (column "site" = normalized station code, to match against STATIONS).

def _gcnet_points_for_station(gcnet_df, station) -- returns the GPS points sorted by date for a given station, or None if no measurement exists for it.

def _gcnet_baseline(t_gps, h_gps, ref_date, use_first_point) -- determines the reference (baseline) value used to align the GC-Net GPS points, either the first available point or the one nearest ref_date.
"""

# ── Inter-dataset statistics ──────────────────────────────────────────────
"""
def inter_dataset_std_before_after(STATIONS, station_name, list_satellites, ...) -- compares the inter-product standard deviation across the requested satellites before vs after a cutoff, for one station.
"""

# ── Per-station AWS uncertainty plots ─────────────────────────────────────
"""
def _plot_aws_variable_on_ax(ax, station, info, ...) -- draws one AWS variable (with its uncertainty) for one station on a given matplotlib axis; shared by plot_aws_variable_uncertainty and plot_aws_variable_uncertainty_single.

def plot_aws_variable_uncertainty(STATIONS, fig_name, aws_variable, ...) -- grid of subplots, one per station, of _plot_aws_variable_on_ax.

def plot_aws_variable_uncertainty_single(station, info, fig_name, ...) -- single-station version of plot_aws_variable_uncertainty, one figure for one station.
"""

# ── Paths ─────────────────────────────────────────────────────────────────
"""
FIGURES_DIR -- output folder for this module's figures, imported from paths.py's shared FIGURES_DIR.
"""


FIGURES_DIR = _BASE_FIGURES_DIR
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def monthly_std_hourly(
    hourly_path,
    columns=("z_surf_combined", "gps_alt", "snow_height"),
    min_valid_fraction=0.0,
):
    """Computes the monthly standard deviation of several daily series from
    an AWS file, grouped by actual year-month (as opposed to a simple
    calendar-month grouping, which would merge January 2015 and January
    2016 together, for instance).

    input:
        hourly_path : str or Path, path to the hourly CSV (must contain a
            'time' column plus the columns listed in `columns`)
        columns : tuple[str], names of the columns to process
            (default: 'z_surf_combined', 'gps_alt', 'snow_height')
        min_valid_fraction : float between 0 and 1, minimum fraction of
            non-NaN values required in a month to compute a std for it
            (NaN otherwise). 0.0 = no filtering (default behaviour).

    output:
        pd.DataFrame indexed by month (Timestamp, start of month), with one
        standard-deviation column per requested variable (e.g.
        'z_surf_combined', 'gps_alt'). NaNs are ignored in the computation
        (np.nanstd), unless the month is empty or insufficiently filled.
    """
    daily_df = pd.read_csv(hourly_path, parse_dates=["time"])
    daily_df = daily_df.set_index("time")
    daily_df = daily_df.sort_index()

    def _std_or_nan(group):
        n_valid = group.notna().sum()
        n_total = len(group)
        if n_total == 0:
            return np.nan
        if min_valid_fraction > 0.0 and (n_valid / n_total) < min_valid_fraction:
            return np.nan
        return np.nanstd(group.values)

    results = {}
    for col in columns:
        if col not in daily_df.columns:
            raise KeyError(f"Column '{col}' missing from {hourly_path}")
        results[col] = daily_df[col].resample("MS").apply(_std_or_nan)

    monthly = pd.DataFrame(results)
    return monthly


def z_surf_combined_anomaly(aws_path):
    aws = pd.read_csv(aws_path, parse_dates=["time"])
    aws["z_surf_combined_anomaly"] = (
        aws["z_surf_combined"] - aws["z_surf_combined"].dropna().iloc[0]
    )
    return aws


def gps_alt_anomaly(aws_path):
    aws = pd.read_csv(aws_path, parse_dates=["time"])
    aws["gps_alt_anomaly"] = aws["gps_alt"] - aws["gps_alt"].dropna().iloc[0]
    return aws


def snow_height(aws_path):
    aws = pd.read_csv(aws_path, parse_dates=["time"])
    return aws


AWS_BASE_COLUMN = {
    "z_surf_combined_anomaly": "z_surf_combined",
    "gps_alt_anomaly": "gps_alt",
    "snow_height": "snow_height",
}

# Per-station date window used to compute Delta z_surf_combined (baseline =
# value at/after "start", endpoint = value at/before "end"). Passed to
# plot_aws_variable_uncertainty() via its `dates_by_station` argument.
# Stations not listed here fall back to that call's global start_date/end_date.
DATES = {
    "JAR": {"start": "2010-09-01", "end": "2022-07-01"},
    "KAN_L": {"start": "2010-09-01", "end": "2023-07-01"},
    "KAN_M": {"start": "2010-09-01", "end": "2023-08-01"},
    "KPC_L": {"start": "2012-07-01", "end": "2023-08-01"},
    "KPC_U": {"start": "2010-09-01", "end": "2023-08-01"},
    "NUK_L": {"start": "2014-08-01", "end": "2023-08-01"},
    "NUK_U": {"start": "2013-08-01", "end": "2023-08-01"},
    "QAS_A": {"start": "2012-08-01", "end": "2015-08-01"},
    "QAS_L": {"start": "2010-11-01", "end": "2023-08-01"},
    "QAS_M": {"start": "2016-08-01", "end": "2022-08-01"},
    "QAS_U": {"start": "2010-10-01", "end": "2023-08-01"},
    "SWC": {"start": "2010-09-01", "end": "2022-07-01"},
    "TAS_A": {"start": "2013-08-01", "end": "2023-08-01"},
    "TAS_L": {"start": "2010-09-01", "end": "2023-08-01"},
    "TAS_U": {"start": "2010-09-01", "end": "2015-08-01"},
    "THU_L": {"start": "2010-09-01", "end": "2023-08-01"},
    "UPE_L": {"start": "2010-09-01", "end": "2023-08-01"},
    "UPE_U": {"start": "2010-09-01", "end": "2023-08-01"},
}


# with ProgressBar():
def plots_variables(
    STATIONS,
    fig_name,
    list_satellites,
    variables,
    colors,
    aws_variables=None,
    aws_colors=None,
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
    show_uncertainty=True,
    aws_hourly_file_key="hourly_data",
    aws_error_style="band",
    aws_min_valid_fraction=0.5,
    aws_errorbar_step=1,
):
    """STATIONS : dict{dict} , fig_name : str , list_satellites : list[str] ,
    variables : list[str] , colors : list[str] , ref_date : str or None
    aws_variables : list[str] or None, among 'z_surf_combined_anomaly'
    and/or 'gps_alt_anomaly'
    aws_colors : list[str] or None, colors matching aws_variables (same order)
    If ref_date is given (e.g. '2010-01-15'), every series (satellite AND
    AWS) is re-based to equal 0 at the date nearest to ref_date. If None,
    each series starts at 0 at its own first point.
    start_date / end_date : time bounds of the plot (str or None).
    Defaults to start_date='2010-11-01' and end_date=None (= today).

    show_uncertainty : bool, if True shows the uncertainty (real monthly
        standard deviation, computed via monthly_std on the hourly data) as
        a band around each AWS curve. Requires
        STATIONS[station][aws_hourly_file_key] to exist; if missing for a
        station, the uncertainty is simply omitted for it.
    aws_hourly_file_key : str, key in info (STATIONS dict entry) pointing to
        the hourly file used for monthly_std.
    aws_error_style : str or dict, 'band' (shaded area), 'errorbar' (spaced
        error bars), or 'both'. Can be a single global str or a dict
        {aws_var: style}.
    aws_min_valid_fraction : float between 0 and 1, passed to monthly_std.
    aws_errorbar_step : int, plot one error bar every N (days)."""

    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
    N = len(STATIONS)
    NCOLS = 3
    NROWS = (N + NCOLS - 1) // NCOLS

    # ── Time bounds ──────────────────────────────────────────────────────
    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()

    # ── Normalize the AWS error style into a dict {aws_var: style} ────────
    def _style_for(aws_var):
        if isinstance(aws_error_style, dict):
            return aws_error_style.get(aws_var, None)
        return aws_error_style

    fig, axes = plt.subplots(
        NROWS,
        NCOLS,
        figsize=(14, NROWS * 3.5),
        constrained_layout=True,
    )
    axes_flat = axes.flatten()
    for ax, (station, info) in zip(axes_flat, STATIONS.items()):
        # ── Satellites ──────────────────────────────────────────────────
        mean_signal_for_gcnet = None

        for k in range(len(list_satellites)):
            Satellite = interpolation_altimetry_AWS.satellite_on_aws(
                info["file"], list_satellites[k], variables[k]
            )
            series = Satellite[variables[k]].dropna()
            series1 = Satellite.iloc[series.index.tolist()]
            time_sat = pd.to_datetime(series1["time_sat"])
            values = series.values

            # time filtering
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
            ax.plot(
                time_sat,
                values_aligned,
                linewidth=1.4,
                color=colors[k],
                label=list_satellites[k],
            )

        # ── AWS (z_surf_combined_anomaly / gps_alt_anomaly) ──────────────
        if aws_variables is not None:
            aws_funcs = {
                "z_surf_combined_anomaly": z_surf_combined_anomaly,
                "gps_alt_anomaly": gps_alt_anomaly,
                "snow_height": snow_height,
            }
            for j, aws_var in enumerate(aws_variables):
                aws_df = aws_funcs[aws_var](info["file"])
                aws_series = aws_df[aws_var].dropna()
                time_aws = aws_df.loc[aws_series.index, "time"]
                values_aws = aws_series.values

                # time filtering
                mask_aws = pd.Series(True, index=time_aws.index)
                if start_bound is not None:
                    mask_aws &= time_aws >= start_bound
                if end_bound is not None:
                    mask_aws &= time_aws <= end_bound
                time_aws = time_aws[mask_aws]
                values_aws = values_aws[mask_aws.values]

                if len(values_aws) == 0:
                    continue

                if ref_date is not None:
                    ref = pd.to_datetime(ref_date)
                    closest_idx = (time_aws - ref).abs().argmin()
                    baseline_aws = (
                        values_aws[time_aws.index.get_loc(closest_idx)]
                        if closest_idx in time_aws.index
                        else values_aws[0]
                    )
                else:
                    baseline_aws = values_aws[0]
                values_aws_aligned = values_aws - baseline_aws
                color_aws = aws_colors[j] if aws_colors is not None else None

                # ── push gps_alt_anomaly into the background ────────────
                is_background = aws_var == "gps_alt_anomaly"
                line_zorder = 0 if is_background else 3
                line_alpha = 0.5 if is_background else 1.0
                band_zorder = -1 if is_background else 1

                ax.plot(
                    time_aws,
                    values_aws_aligned,
                    linewidth=1.4,
                    linestyle="--",
                    color=color_aws,
                    alpha=line_alpha,
                    zorder=line_zorder,
                    label=aws_var,
                )

                # ── Dynamic uncertainty (real monthly std) ───────────────
                hourly_path = info.get(aws_hourly_file_key)
                if show_uncertainty and hourly_path is not None:
                    base_col = AWS_BASE_COLUMN.get(aws_var)
                    std_df = monthly_std_hourly(
                        hourly_path,
                        columns=(base_col,),
                        min_valid_fraction=aws_min_valid_fraction,
                    )
                    sigma_lookup = std_df[base_col]
                    sigma_lookup.index = sigma_lookup.index.to_period("M")

                    day_periods = pd.DatetimeIndex(time_aws).to_period("M")
                    sigma_values = sigma_lookup.reindex(day_periods).values

                    valid_sigma = ~pd.isna(sigma_values)
                    if valid_sigma.any():
                        style = _style_for(aws_var)
                        t_arr = pd.DatetimeIndex(time_aws)[valid_sigma]
                        v_arr = np.asarray(values_aws_aligned)[valid_sigma]
                        s_arr = sigma_values[valid_sigma].astype(float)

                        if style in ("band", "both"):
                            ax.fill_between(
                                t_arr,
                                v_arr - s_arr,
                                v_arr + s_arr,
                                color=color_aws,
                                alpha=0.15,
                                linewidth=0,
                                zorder=band_zorder,
                            )

                        if style in ("errorbar", "both"):
                            step = max(1, aws_errorbar_step)
                            ax.errorbar(
                                t_arr[::step],
                                v_arr[::step],
                                yerr=s_arr[::step],
                                fmt="none",
                                ecolor=color_aws,
                                elinewidth=0.8,
                                capsize=2,
                                alpha=0.7,
                                zorder=line_zorder,
                            )

        ax.legend(fontsize=7)
        ax.set_title(station, fontsize=12, fontweight="bold")
        ax.set_xlabel("Years", fontsize=9)
        ax.set_ylabel("Surface elevation changes (m)", fontsize=9)
        ax.tick_params(axis="both", labelsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_major_locator(mdates.YearLocator(1))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    for ax in axes_flat[N:]:
        ax.set_visible(False)
    fig.suptitle(fig_name, fontsize=15, fontweight="bold")
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
    print(f"Figure saved: {OUTPUT_FILE}")
    plt.show()


# plots_variables(
#     AWS_data.STATIONS_ablation, 'Surface Elevation Change for ablation areas from 2010-11-01 to now with z_surf_combined_anomaly (AWS)',
#     ['Copernicus_Climate_Data_Store', 'Nilsson and Gardner, 2026', 'Andersen et al., 2025', 'Khan et al., 2025', 'Zhang et al., 2022'],
#     ['dh', 'dh', 'ZZ', 'dh_vol', 'elev_interp'],
#     ["#EB2F25", "#17EE17", "#EDE20D", "#D20DF0", "#03FEF1"],
#     aws_variables=['z_surf_combined_anomaly'],
#     aws_colors=["#000000"]
# )


# _________________________________________________________________________________________________________________________________________________________


def plots_variables_with_satellite_uncertainty(
    STATIONS,
    fig_name,
    list_satellites,
    variables,
    colors,
    uncertainty_variables=None,
    ref_date=None,
    start_date=None,
    end_date=None,
    sat_band_alpha=0.15,
):
    """Plots, for each station, the satellite series (list_satellites/
    variables/colors) with an uncertainty band (value +/- uncertainty)
    around each curve.

    STATIONS : dict{dict}, each entry must contain 'file'
    fig_name : str, figure title
    list_satellites : list[str]
    variables : list[str], same length as list_satellites
    colors : list[str], same length as list_satellites
    uncertainty_variables : list[str or None] or None, same length as
        list_satellites. Name of the uncertainty variable to use for each
        satellite (e.g. 'rms' for Nilsson), or None if no uncertainty is
        available for that satellite (e.g. Khan).
    ref_date : str or None. If given (e.g. '2010-01-15'), every series is
        re-based to equal 0 at the date nearest to ref_date. If None, each
        series starts at 0 at its own first point.
    start_date / end_date : time bounds of the plot (str or None).
        Defaults to start_date='2010-11-01' and end_date=None (= today).
    sat_band_alpha : float, transparency of the uncertainty bands."""

    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
    N = len(STATIONS)
    NCOLS = 2
    NROWS = (N + NCOLS - 1) // NCOLS

    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()

    fig, axes = plt.subplots(
        NROWS,
        NCOLS,
        figsize=(14, NROWS * 3.5),
        constrained_layout=True,
    )
    axes_flat = np.atleast_1d(axes).ravel()
    for ax, (station, info) in zip(axes_flat, STATIONS.items()):
        for k in range(len(list_satellites)):
            Satellite = interpolation_altimetry_AWS.satellite_on_aws(
                info["file"], list_satellites[k], variables[k]
            )
            series = Satellite[variables[k]].dropna()
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
                continue

            if ref_date is not None:
                ref = pd.to_datetime(ref_date)
                closest_idx = (time_sat - ref).abs().argmin()
                baseline = values[closest_idx]
            else:
                baseline = values[0]
            values_aligned = values - baseline

            ax.plot(
                time_sat,
                values_aligned,
                linewidth=1.4,
                linestyle="-",
                alpha=1,
                color=colors[k],
                label=list_satellites[k],
                zorder=3,
            )

            # ── Satellite uncertainty band ───────────────────────────────
            unc_var = (
                uncertainty_variables[k] if uncertainty_variables is not None else None
            )
            if unc_var is not None:

                def _to_scalar(v):
                    try:
                        return float(v)
                    except (TypeError, ValueError):
                        return np.nan

                Satellite_unc = interpolation_altimetry_AWS.satellite_on_aws(
                    info["file"], list_satellites[k], unc_var
                )
                unc_raw = Satellite_unc[unc_var]
                unc_values = np.array(
                    [_to_scalar(v) for v in unc_raw.values], dtype="float64"
                )
                unc_time = pd.to_datetime(Satellite_unc["time_sat"])

                unc_df = pd.DataFrame(
                    {"time": unc_time.values, "unc": unc_values}
                ).dropna()
                unc_df["time"] = unc_df["time"].astype("datetime64[ns]")

                main_df = pd.DataFrame(
                    {"time": time_sat.values, "value": values_aligned}
                )
                main_df["time"] = main_df["time"].astype("datetime64[ns]")

                merged = pd.merge(main_df, unc_df, on="time", how="inner")

                if not merged.empty:
                    ax.fill_between(
                        merged["time"],
                        merged["value"] - merged["unc"],
                        merged["value"] + merged["unc"],
                        color=colors[k],
                        alpha=sat_band_alpha,
                        linewidth=0,
                        zorder=1,
                    )

        ax.legend(fontsize=7)
        ax.set_title(station, fontsize=12, fontweight="bold")
        ax.set_xlabel("Years", fontsize=9)
        ax.set_ylabel("Surface elevation changes (m)", fontsize=9)
        ax.tick_params(axis="both", labelsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_major_locator(mdates.YearLocator(1))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    for ax in axes_flat[N:]:
        ax.set_visible(False)
    fig.suptitle(fig_name, fontsize=15, fontweight="bold")
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
    print(f"Figure saved: {OUTPUT_FILE}")
    plt.show()


# plots_variables_with_satellite_uncertainty(
#     AWS_data.STATION_Nilsson,
#     'Surface Elevation Change according to Nilsson and Garnder, 2026 - uncertainty bands',
#     ['Nilsson and Gardner, 2026'], #'Copernicus_Climate_Data_Store', , 'Andersen et al., 2025', 'Khan et al., 2025', 'Zhang et al., 2022'
#     ['dh'], #'dh', , 'ZZ', 'dh_vol', 'elev_interp'
#     ["#17EE17"], #"#EB2F25", , "#EDE20D", "#D20DF0", "#03FEF1"
#     uncertainty_variables=['rms'], #'dh_uncert', , 'ZZer', None, 'elev_uncer_interp'
# )


def _multi_satellite_std_band(
    info,
    list_satellites,
    variables,
    start_bound,
    end_bound,
    ref_date,
    resample_freq="MS",
):
    """Computes, for a station, the mean and inter-satellite standard
    deviation of the aligned series, on a common time grid (monthly by
    default). Only keeps points where at least 2 satellites are present.

    Returns (mean_signal, std_signal) as pd.Series indexed by time, or None
    if fewer than 2 satellites have usable data."""
    series_dict = {}
    for k in range(len(list_satellites)):
        Satellite = interpolation_altimetry_AWS.satellite_on_aws(
            info["file"], list_satellites[k], variables[k]
        )
        series = Satellite[variables[k]].dropna()
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
            continue

        if ref_date is not None:
            ref = pd.to_datetime(ref_date)
            closest_idx = (time_sat - ref).abs().argmin()
            baseline = values[closest_idx]
        else:
            baseline = values[0]
        values_aligned = values - baseline

        s = pd.Series(values_aligned, index=pd.DatetimeIndex(time_sat)).sort_index()
        series_dict[list_satellites[k]] = s.resample(resample_freq).mean()

    if len(series_dict) < 2:
        return None

    combined = pd.DataFrame(series_dict)
    n_valid = combined.notna().sum(axis=1)
    valid_mask = n_valid >= 2

    mean_signal = combined.mean(axis=1, skipna=True)[valid_mask]
    std_signal = combined.std(axis=1, skipna=True, ddof=1)[valid_mask]

    if mean_signal.empty:
        return None
    return mean_signal, std_signal


# ──────────────────────────────────────────────────────────────────────
# GC-Net precise GPS points
# ──────────────────────────────────────────────────────────────────────
def load_gcnet_precise_locations(csv_path):
    """Loads the GC-Net precise GPS measurements CSV.
    Column 'site' = normalized station code (to match against STATIONS)."""
    df = pd.read_csv(csv_path)
    df["date"] = pd.to_datetime(df["date"], format="%d-%m-%Y")
    return df


def _gcnet_points_for_station(gcnet_df, station):
    """Returns the GPS points sorted by date for a given station, or None
    if no measurement exists for it."""
    pts = gcnet_df[gcnet_df["site"] == station].copy()
    if pts.empty:
        return None
    return pts.sort_values("date")


def _gcnet_baseline(t_gps, h_gps, ref_date, use_first_point=True):
    """Determines the reference (baseline) value used to align the GPS
    points.
      - if use_first_point=True (default): always the first available GPS
        measurement within the filtered range, regardless of ref_date.
        This is the desired behaviour to compute the "clean" GNSS anomaly
        relative to its own first point.
      - if use_first_point=False and ref_date is given: height of the GPS
        point nearest to ref_date (legacy behaviour, to visually align the
        GPS points onto the other series).
      - otherwise: first available GPS measurement.
    """
    if not use_first_point and ref_date is not None:
        ref = pd.to_datetime(ref_date)
        closest_idx = (t_gps - ref).abs().values.argmin()
        return h_gps[closest_idx]
    return h_gps[0]


def plots_variables_full_uncertainty(
    STATIONS,
    fig_name,
    list_satellites,
    variables,
    colors,
    uncertainty_variables=None,
    aws_variables=None,
    aws_colors=None,
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
    sat_band_alpha=0.25,
    show_uncertainty=True,
    aws_hourly_file_key="hourly_data",
    aws_error_style="band",
    aws_min_valid_fraction=0.5,
    aws_errorbar_step=1,
    show_multi_sat_std_band=True,
    multi_sat_std_resample="MS",
    multi_sat_std_color="#527bff",
    multi_sat_sigma_alphas=(0.30, 0.15),
    multi_sat_legend_alpha=0.55,
    show_multi_sat_mean=True,
    multi_sat_mean_color="#1000f3",
    gcnet_csv_path=PRODEM_DATA_FILE,
    gcnet_color="black",
    gcnet_marker="o",
    gcnet_markersize=30,
    gcnet_baseline_first_point=True,
    gcnet_align_to_mean=True,
    multi_sat_sigma_label_stat="mean",
    multi_sat_sigma_unit="m",
    multi_sat_sigma_legend_loc="upper left",
):
    """Plots, for each station, the satellite AND AWS series, with an
    uncertainty band around EVERY curve (satellite and AWS), an
    inter-satellite dispersion band (1/2 sigma) with its mean, and, overlaid,
    the GC-Net precise GPS points (when available for the station).

    -- Satellite parameters --
    STATIONS : dict{dict}, each entry must contain 'file'
    fig_name : str, figure title
    list_satellites : list[str]
    variables : list[str], same length as list_satellites
    colors : list[str], same length as list_satellites
    uncertainty_variables : list[str or None] or None, same length as
        list_satellites.

    -- AWS parameters --
    aws_variables : list[str] or None, among 'z_surf_combined_anomaly'
        and/or 'gps_alt_anomaly'
    aws_colors : list[str] or None, colors matching aws_variables (same order)
    show_uncertainty : bool, shows the AWS uncertainty (real monthly std via
        monthly_std_hourly).
    aws_hourly_file_key : str, key in info pointing to the hourly file.
    aws_error_style : str or dict, 'band', 'errorbar', or 'both'.
    aws_min_valid_fraction : float between 0 and 1, passed to
        monthly_std_hourly.
    aws_errorbar_step : int, plot one error bar every N points.

    -- Inter-satellite dispersion --
    show_multi_sat_std_band : bool, shows the 1/2 sigma band.
    multi_sat_std_resample : str, common resampling frequency.
    multi_sat_std_color : band color.
    multi_sat_sigma_alphas : tuple (alpha_1sigma, alpha_2sigma).
    show_multi_sat_mean : bool, shows the inter-satellite mean.
    multi_sat_mean_color : color of that mean.

    -- GC-Net precise GPS points --
    gcnet_csv_path : str or Path or None. If given, loads the precise GPS
        measurements CSV and overlays, for each station present both in
        STATIONS and in the CSV ('site' column), points (marker).
    gcnet_color, gcnet_marker, gcnet_markersize : GPS point style.
    gcnet_baseline_first_point : bool, if True (default) the GNSS anomaly is
        computed relative to the station's first available GPS point,
        independently of ref_date. If False, uses ref_date (legacy
        behaviour, visually consistent with satellites/AWS).

    -- Common parameters --
    ref_date : str or None. If given, satellites and AWS are re-based to
        equal 0 at the date nearest to ref_date. If None, each series starts
        at 0 at its own first point.
    start_date / end_date : time bounds of the plot (str or None).
    sat_band_alpha : float, transparency of the satellite uncertainty bands.
    """

    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
    N = len(STATIONS)
    NCOLS = 2
    NROWS = (N + NCOLS - 1) // NCOLS

    # ── Time bounds ──────────────────────────────────────────────────────
    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()

    # ── Load the GC-Net CSV once (outside the station loop) ────────────
    gcnet_df = (
        load_gcnet_precise_locations(gcnet_csv_path)
        if gcnet_csv_path is not None
        else None
    )

    # ── Normalize the AWS error style into a dict {aws_var: style} ──────
    def _style_for(aws_var):
        if isinstance(aws_error_style, dict):
            return aws_error_style.get(aws_var, None)
        return aws_error_style

    def _sigma_repr_value(std_signal, stat="mean"):
        """Scalar representative sigma value for the local legend."""
        vals = std_signal.dropna().values
        if len(vals) == 0:
            return np.nan
        if stat == "median":
            return float(np.nanmedian(vals))
        if stat == "last":
            return float(vals[-1])
        return float(np.nanmean(vals))

    def _to_scalar(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan

    fig, axes = plt.subplots(
        NROWS,
        NCOLS,
        sharey=False,
        figsize=(14, NROWS * 3.5),
        constrained_layout=True,
    )
    axes_flat = np.atleast_1d(axes).ravel()
    for ax, (station, info) in zip(axes_flat, STATIONS.items()):
        # ── Inter-satellite standard-deviation band (background) ────────
        mean_signal_for_gcnet = None
        need_multi_sat = (
            (show_multi_sat_std_band or gcnet_align_to_mean)
            and list_satellites is not None
            and len(list_satellites) >= 2
        )
        if need_multi_sat:
            band_result = _multi_satellite_std_band(
                info,
                list_satellites,
                variables,
                start_bound,
                end_bound,
                ref_date,
                resample_freq=multi_sat_std_resample,
            )
            if band_result is not None:
                mean_signal, std_signal = band_result
                mean_signal_for_gcnet = mean_signal
                t_band = mean_signal.index

                if show_multi_sat_std_band:
                    for sigma_mult, alpha in zip(
                        (2, 1), reversed(multi_sat_sigma_alphas)
                    ):
                        ax.fill_between(
                            t_band,
                            mean_signal - sigma_mult * std_signal,
                            mean_signal + sigma_mult * std_signal,
                            color=multi_sat_std_color,
                            alpha=alpha,
                            linewidth=0,
                            zorder=-10 + sigma_mult,
                            label=None,
                        )
                        ax.plot(
                            t_band,
                            mean_signal - sigma_mult * std_signal,
                            color=multi_sat_std_color,
                            linewidth=0.7,
                            linestyle=(0, (3, 2)),
                            alpha=0.6,
                            zorder=-6,
                        )
                        ax.plot(
                            t_band,
                            mean_signal + sigma_mult * std_signal,
                            color=multi_sat_std_color,
                            linewidth=0.7,
                            linestyle=(0, (3, 2)),
                            alpha=0.6,
                            zorder=-6,
                        )

                    if show_multi_sat_mean:
                        ax.plot(
                            t_band,
                            mean_signal,
                            color=multi_sat_mean_color,
                            linewidth=2,
                            linestyle="-",
                            alpha=0.9,
                            zorder=10,
                            label="Inter-datasets average",
                        )

                    sigma_repr = _sigma_repr_value(
                        std_signal, stat=multi_sat_sigma_label_stat
                    )
                    if not np.isnan(sigma_repr):
                        sigma_handles = [
                            Patch(
                                facecolor=multi_sat_std_color,
                                alpha=a,
                                edgecolor=multi_sat_std_color,
                                linewidth=0.7,
                                label=f"{s}σ ≈ {s * sigma_repr:.3f} {multi_sat_sigma_unit}",
                            )
                            for s, a in zip((1, 2), multi_sat_sigma_alphas)
                        ]
                        sigma_legend = ax.legend(
                            handles=sigma_handles,
                            loc=multi_sat_sigma_legend_loc,
                            fontsize=6,
                            framealpha=multi_sat_legend_alpha,
                            title="Inter-sat. dispersion",
                            title_fontsize=6,
                            borderpad=0.4,
                            handlelength=1.2,
                            labelspacing=0.3,
                        )
                        ax.add_artist(sigma_legend)

        # ── Satellites ────────────────────────────────────────────────
        if list_satellites is not None:
            for k in range(len(list_satellites)):
                Satellite = interpolation_altimetry_AWS.satellite_on_aws(
                    info["file"], list_satellites[k], variables[k]
                )
                series = Satellite[variables[k]].dropna()
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
                    continue

                if ref_date is not None:
                    ref = pd.to_datetime(ref_date)
                    closest_idx = (time_sat - ref).abs().argmin()
                    baseline = values[closest_idx]
                else:
                    baseline = values[0]
                values_aligned = values - baseline

                ax.plot(
                    time_sat,
                    values_aligned,
                    linewidth=1.4,
                    linestyle="-",
                    alpha=0.8,
                    color=colors[k],
                    label=list_satellites[k],
                    zorder=3,
                )

                # ── Satellite uncertainty band ──────────────────────────
                unc_var = (
                    uncertainty_variables[k]
                    if uncertainty_variables is not None
                    else None
                )
                if unc_var is not None:
                    Satellite_unc = interpolation_altimetry_AWS.satellite_on_aws(
                        info["file"], list_satellites[k], unc_var
                    )
                    unc_raw = Satellite_unc[unc_var]
                    unc_values = np.array(
                        [_to_scalar(v) for v in unc_raw.values], dtype="float64"
                    )
                    unc_time = pd.to_datetime(Satellite_unc["time_sat"])

                    unc_df = pd.DataFrame(
                        {"time": unc_time.values, "unc": unc_values}
                    ).dropna()
                    unc_df["time"] = unc_df["time"].astype("datetime64[ns]")

                    main_df = pd.DataFrame(
                        {"time": time_sat.values, "value": values_aligned}
                    )
                    main_df["time"] = main_df["time"].astype("datetime64[ns]")

                    merged = pd.merge(main_df, unc_df, on="time", how="inner")

                    if not merged.empty:
                        ax.fill_between(
                            merged["time"],
                            merged["value"] - merged["unc"],
                            merged["value"] + merged["unc"],
                            color=colors[k],
                            alpha=sat_band_alpha,
                            linewidth=0,
                            zorder=1,
                        )

        # ── AWS (z_surf_combined_anomaly / gps_alt_anomaly) ─────────────
        if aws_variables is not None:
            aws_funcs = {
                "z_surf_combined_anomaly": z_surf_combined_anomaly,
                "gps_alt_anomaly": gps_alt_anomaly,
                "snow_height": snow_height,
            }
            for j, aws_var in enumerate(aws_variables):
                aws_df = aws_funcs[aws_var](info["file"])
                aws_series = aws_df[aws_var].dropna()
                time_aws = aws_df.loc[aws_series.index, "time"]
                values_aws = aws_series.values

                mask_aws = pd.Series(True, index=time_aws.index)
                if start_bound is not None:
                    mask_aws &= time_aws >= start_bound
                if end_bound is not None:
                    mask_aws &= time_aws <= end_bound
                time_aws = time_aws[mask_aws]
                values_aws = values_aws[mask_aws.values]

                if len(values_aws) == 0:
                    continue

                if ref_date is not None:
                    ref = pd.to_datetime(ref_date)
                    closest_idx = (time_aws - ref).abs().argmin()
                    baseline_aws = (
                        values_aws[time_aws.index.get_loc(closest_idx)]
                        if closest_idx in time_aws.index
                        else values_aws[0]
                    )
                else:
                    baseline_aws = values_aws[0]
                values_aws_aligned = values_aws - baseline_aws
                color_aws = aws_colors[j] if aws_colors is not None else None

                is_background = aws_var == "gps_alt_anomaly"
                line_zorder = 0 if is_background else 3
                line_alpha = 0.5 if is_background else 1.0
                band_zorder = -1 if is_background else 1

                ax.plot(
                    time_aws,
                    values_aws_aligned,
                    linewidth=1.6,
                    linestyle="--",
                    color=color_aws,
                    alpha=line_alpha,
                    zorder=20,
                    label=aws_var,
                )

                # ── Dynamic uncertainty (real monthly std) ───────────────
                hourly_path = info.get(aws_hourly_file_key)
                if show_uncertainty and hourly_path is not None:
                    base_col = AWS_BASE_COLUMN.get(aws_var)
                    std_df = monthly_std_hourly(
                        hourly_path,
                        columns=(base_col,),
                        min_valid_fraction=aws_min_valid_fraction,
                    )
                    sigma_lookup = std_df[base_col]
                    sigma_lookup.index = sigma_lookup.index.to_period("M")

                    day_periods = pd.DatetimeIndex(time_aws).to_period("M")
                    sigma_values = sigma_lookup.reindex(day_periods).values

                    valid_sigma = ~pd.isna(sigma_values)
                    if valid_sigma.any():
                        style = _style_for(aws_var)
                        t_arr = pd.DatetimeIndex(time_aws)[valid_sigma]
                        v_arr = np.asarray(values_aws_aligned)[valid_sigma]
                        s_arr = sigma_values[valid_sigma].astype(float)

                        if style in ("band", "both"):
                            ax.fill_between(
                                t_arr,
                                v_arr - s_arr,
                                v_arr + s_arr,
                                color=color_aws,
                                alpha=0.15,
                                linewidth=0,
                                zorder=band_zorder,
                            )

                        if style in ("errorbar", "both"):
                            step = max(1, aws_errorbar_step)
                            ax.errorbar(
                                t_arr[::step],
                                v_arr[::step],
                                yerr=s_arr[::step],
                                fmt="none",
                                ecolor=color_aws,
                                elinewidth=0.8,
                                capsize=2,
                                alpha=0.7,
                                zorder=line_zorder,
                            )

        # ── GC-Net precise GPS points (anomaly / first point or ref_date) ─
        if gcnet_df is not None:
            pts = _gcnet_points_for_station(gcnet_df, station)
            if pts is not None:
                t_gps = pts["date"]
                h_gps = pts["orthometric_height_m"].values

                mask_gps = pd.Series(True, index=range(len(t_gps)))
                if start_bound is not None:
                    mask_gps &= (t_gps >= start_bound).values
                if end_bound is not None:
                    mask_gps &= (t_gps <= end_bound).values
                t_gps_f = pd.DatetimeIndex(t_gps[mask_gps.values])
                h_gps_f = h_gps[mask_gps.values]

                if len(h_gps_f) > 0:
                    baseline_gps = _gcnet_baseline(
                        t_gps_f,
                        h_gps_f,
                        ref_date,
                        use_first_point=gcnet_baseline_first_point,
                    )
                    h_gps_aligned = h_gps_f - baseline_gps

                    # ── Align the first GNSS point onto the inter-sat mean ──
                    if (
                        gcnet_align_to_mean
                        and mean_signal_for_gcnet is not None
                        and len(mean_signal_for_gcnet) > 0
                    ):
                        t0 = t_gps_f[0]
                        nearest_pos = mean_signal_for_gcnet.index.get_indexer(
                            [t0], method="nearest"
                        )[0]
                        if nearest_pos != -1:
                            mean_val_at_t0 = mean_signal_for_gcnet.iloc[nearest_pos]
                            h_gps_aligned = h_gps_aligned + mean_val_at_t0

                    ax.scatter(
                        t_gps_f,
                        h_gps_aligned,
                        marker=gcnet_marker,
                        s=gcnet_markersize,
                        color=gcnet_color,
                        edgecolor="white",
                        linewidth=0.6,
                        zorder=15,
                        label="GC-Net GPS",
                    )

        ax.set_title(station, fontsize=12, fontweight="bold")
        ax.set_xlabel("Years", fontsize=9)
        ax.set_ylabel("Surface elevation changes (m)", fontsize=9)
        ax.tick_params(axis="both", labelsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_major_locator(mdates.YearLocator(1))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    # ── END of the station loop: runs ONCE ────────────────────────────────

    # ── Single deduplicated legend for the whole figure ───────────────────
    all_handles, all_labels = [], []
    for ax in axes_flat[:N]:
        h, l = ax.get_legend_handles_labels()
        for hh, ll in zip(h, l):
            if ll not in all_labels:
                all_handles.append(hh)
                all_labels.append(ll)

    if all_handles:
        fig.legend(
            all_handles,
            all_labels,
            loc="lower center",
            ncol=min(len(all_labels), 6),
            fontsize=8,
            bbox_to_anchor=(0.5, -0.03),
            frameon=True,
        )

    for ax in axes_flat[N:]:
        ax.set_visible(False)
    fig.suptitle(fig_name, fontsize=15, fontweight="bold")
    plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
    print(f"Figure saved: {OUTPUT_FILE}")
    plt.show()


plots_variables_full_uncertainty(
    AWS_data.STATION_use,
    "Surface Elevation Change - SDL and NSE - from 2011-01-01 - multi sat. stantard deviation bands and mean",
    ['Copernicus_Climate_Data_Store', 'Andersen et al., 2025', 'Nilsson and Gardner, 2026', 'Khan et al., 2025', 'Zhang et al., 2022'], # 
    ['dh', 'ZZ', 'dh', 'dh_vol', 'elev_interp'], #
    ["#EB2F25", "#EDE20D", "#17EE17","#D20DF0", "#03FEF1"], #
    uncertainty_variables=None,    #,['dh_uncert', 'ZZer', 'rms', None, 'elev_uncer_interp']
    aws_variables=None, #['gps_alt_anomaly', 'z_surf_combined_anomaly']
    aws_colors=None, #['blue', 'black']
    show_multi_sat_mean=True,
    show_multi_sat_std_band=True,
)

# ___________________________________________________________________________________________________________________________________________________________________


def inter_dataset_std_before_after(
    STATIONS,
    station_name,
    list_satellites,
    variables,
    split_date,
    start_date="2010-08-01",
    end_date=None,
    ref_date=None,
    resample_freq="MS",
):
    """
    Computes the mean inter-dataset standard deviation (dispersion between
    satellites) before and after `split_date`, for a given station.

    Reuses `_multi_satellite_std_band` (the same function used in the
    background of `plots_variables_full_uncertainty`), so the values
    returned here are consistent with the band shown on the figures.

    Returns a dict with the before/after stats (mean, median, number of
    points, start/end dates of each sub-period).
    """
    info = STATIONS[station_name]
    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()
    split = pd.to_datetime(split_date)

    band_result = _multi_satellite_std_band(
        info,
        list_satellites,
        variables,
        start_bound,
        end_bound,
        ref_date,
        resample_freq=resample_freq,
    )
    if band_result is None:
        raise ValueError(
            f"Not enough overlap between satellites to compute inter-dataset "
            f"dispersion for {station_name}."
        )

    mean_signal, std_signal = band_result
    std_signal = std_signal.dropna()

    before = std_signal[std_signal.index < split]
    after = std_signal[std_signal.index >= split]

    def _stats(s):
        return {
            "mean": float(s.mean()) if len(s) else np.nan,
            "median": float(s.median()) if len(s) else np.nan,
            "n_points": int(len(s)),
            "date_min": s.index.min() if len(s) else None,
            "date_max": s.index.max() if len(s) else None,
        }

    result = {"before": _stats(before), "after": _stats(after)}

    print(f"── {station_name} — inter-dataset std, split at {split.date()} ──")
    for period in ("before", "after"):
        st = result[period]
        print(
            f"  {period:6s}: mean={st['mean']:.4f}  median={st['median']:.4f}  "
            f"n={st['n_points']}  range=[{st['date_min']}, {st['date_max']}]"
        )

    return result


# Example usage for SDL, with the same datasets as the existing figure:
# inter_dataset_std_before_after(
#     AWS_data.STATION_CEN,
#     "CEN",
#     ['Copernicus_Climate_Data_Store', 'Nilsson and Gardner, 2026', 'Andersen et al., 2025', 'Khan et al., 2025', 'Zhang et al., 2022'], #
#     ['dh','dh', 'ZZ', 'dh_vol', 'elev_interp'], #
#     split_date="2017-06-01",
# )


# ___________________________________________________________________________________________________________________________________________________________________


def _plot_aws_variable_on_ax(
    ax,
    station,
    info,
    aws_variable,
    aws_funcs,
    base_col,
    color,
    ref_date,
    start_bound_i,
    end_bound_i,
    aws_hourly_file_key,
    error_style,
    min_valid_fraction,
    errorbar_step,
):
    """Draws one station's AWS anomaly curve, its hourly-uncertainty band/
    error bars, and the Delta z_surf legend line onto `ax`. Factored out so
    plot_aws_variable_uncertainty (grid of every station) and
    plot_aws_variable_uncertainty_single (one station, one standalone
    figure) share the exact same per-station drawing logic and stay in
    sync -- see either function's docstring for what each parameter means.

    Returns True if data was plotted, False if the station had no data in
    [start_bound_i, end_bound_i] (the title is still set to
    "{station} (no data)" in that case, so the panel/figure isn't left
    blank without explanation).
    """
    aws_df = aws_funcs[aws_variable](info["file"])
    aws_series = aws_df[aws_variable].dropna()
    time_aws = aws_df.loc[aws_series.index, "time"]
    values_aws = aws_series.values

    # ── time filtering ──────────────────────────────────────────────────
    mask_aws = pd.Series(True, index=time_aws.index)
    if start_bound_i is not None:
        mask_aws &= time_aws >= start_bound_i
    if end_bound_i is not None:
        mask_aws &= time_aws <= end_bound_i
    time_aws = time_aws[mask_aws]
    values_aws = values_aws[mask_aws.values]

    if len(values_aws) == 0:
        ax.set_title(f"{station} (no data)", fontsize=12, fontweight="bold")
        return False

    # ── reset to zero ────────────────────────────────────────────────────
    if ref_date is not None:
        ref = pd.to_datetime(ref_date)
        closest_idx = (time_aws - ref).abs().argmin()
        baseline = (
            values_aws[time_aws.index.get_loc(closest_idx)]
            if closest_idx in time_aws.index
            else values_aws[0]
        )
    else:
        baseline = values_aws[0]
    values_aligned = values_aws - baseline

    ax.plot(
        time_aws,
        values_aligned,
        linewidth=1.6,
        color=color,
        zorder=3,
        label=aws_variable,
    )

    # ── Hourly uncertainty (real monthly std) ────────────────────────────
    hourly_path = info.get(aws_hourly_file_key)
    sigma_values = None
    if hourly_path is None:
        print(f"[{station}] no '{aws_hourly_file_key}' file -- uncertainty omitted.")
    else:
        std_df = monthly_std_hourly(
            hourly_path,
            columns=(base_col,),
            min_valid_fraction=min_valid_fraction,
        )
        sigma_lookup = std_df[base_col]
        sigma_lookup.index = sigma_lookup.index.to_period("M")

        day_periods = pd.DatetimeIndex(time_aws).to_period("M")
        sigma_values = sigma_lookup.reindex(day_periods).values

        valid_sigma = ~pd.isna(sigma_values)
        if valid_sigma.any():
            t_arr = pd.DatetimeIndex(time_aws)[valid_sigma]
            v_arr = np.asarray(values_aligned)[valid_sigma]
            s_arr = sigma_values[valid_sigma].astype(float)

            if error_style in ("band", "both"):
                ax.fill_between(
                    t_arr,
                    v_arr - s_arr,
                    v_arr + s_arr,
                    color=color,
                    alpha=0.2,
                    linewidth=0,
                    zorder=1,
                    label="± hourly std (monthly)",
                )
            if error_style in ("errorbar", "both"):
                step = max(1, errorbar_step)
                ax.errorbar(
                    t_arr[::step],
                    v_arr[::step],
                    yerr=s_arr[::step],
                    fmt="none",
                    ecolor=color,
                    elinewidth=0.8,
                    capsize=2,
                    alpha=0.7,
                    zorder=2,
                )

    # ── Total balance (delta z_surf) over the displayed period ───────────
    delta_z = float(values_aligned[-1])
    if sigma_values is not None:
        sigma_first = sigma_values[0]
        sigma_last = sigma_values[-1]
        terms = [s**2 for s in (sigma_first, sigma_last) if not np.isnan(s)]
        delta_unc = np.sqrt(sum(terms)) if terms else np.nan
    else:
        delta_unc = np.nan

    window_str = f"{time_aws.iloc[0].date()} → {time_aws.iloc[-1].date()}"
    if np.isnan(delta_unc):
        delta_label = f"Δz_surf ({window_str}) = {delta_z:+.2f} m (unc. n/a)"
    else:
        delta_label = f"Δz_surf ({window_str}) = {delta_z:+.2f} ± {delta_unc:.2f} m"
    delta_handle = Line2D([], [], color="none", label=delta_label)

    ax.set_title(station, fontsize=12, fontweight="bold")
    ax.set_xlabel("Years", fontsize=9)
    ax.set_ylabel(f"{aws_variable} (m)", fontsize=9)
    ax.tick_params(axis="both", labelsize=8)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_major_locator(mdates.YearLocator(1))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    handles, labels = ax.get_legend_handles_labels()
    handles.append(delta_handle)
    labels.append(delta_label)
    ax.legend(handles, labels, fontsize=7, loc="best")

    return True


def plot_aws_variable_uncertainty(
    STATIONS,
    fig_name,
    aws_variable="z_surf_combined_anomaly",
    color="black",
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
    dates_by_station=None,
    aws_hourly_file_key="hourly_data",
    error_style="band",
    min_valid_fraction=0.5,
    errorbar_step=1,
    save=True,
    show=True,
):
    """Lightweight version of plots_variables_full_uncertainty for when you
    only want the AWS -- a single variable (default
    z_surf_combined_anomaly) with its uncertainty, no satellites, no
    inter-dataset dispersion, no GC-Net points. Avoids having to call
    plots_variables_full_uncertainty with list_satellites/variables/colors
    all set to None.

    STATIONS : dict{dict}, each entry must contain 'file' and, if
        error_style uses it, `aws_hourly_file_key`.
    fig_name : str, figure title (and name of the saved PNG).
    aws_variable : str, one of 'z_surf_combined_anomaly',
        'gps_alt_anomaly', 'snow_height' (default: 'z_surf_combined_anomaly').
    color : str, color of the curve and of the uncertainty band/error bars.
    ref_date : str or None. If given, the series is re-based to equal 0 at
        the date nearest to ref_date. If None (default), it starts at 0 at
        its own first point.
    start_date / end_date : global time bounds of the plot (str or None),
        used for every station unless overridden by dates_by_station.
        Defaults to start_date='2010-11-01' and end_date=None (= today).
    dates_by_station : dict{station: {"start": str, "end": str}} or None.
        Per-station override of start_date/end_date, used to restrict both
        the plotted period AND the Delta z_surf endpoints for that station
        specifically: the baseline is the first point on/after "start" and
        the delta endpoint is the last point on/before "end" (same logic as
        the global start_date/end_date, just applied per station). Stations
        missing from this dict fall back to the global start_date/end_date.
        Handy when different stations have different usable record periods
        (e.g. a station only covering 2012-2015) -- pass a per-station date
        dict here rather than calling this function once per station with a
        different start_date/end_date each time.
    aws_hourly_file_key : str, key in info pointing to the hourly file used
        by monthly_std_hourly to compute the uncertainty (real monthly
        standard deviation, computed from HOURLY data, not daily). If
        missing for a station, the uncertainty is simply omitted for it (a
        message is printed, no error).
    error_style : 'band' (shaded area, default), 'errorbar' (spaced error
        bars), or 'both'.
    min_valid_fraction : float between 0 and 1, passed to
        monthly_std_hourly (minimum fraction of non-NaN hourly values
        required in a month to compute a standard deviation for it, NaN
        otherwise -- so no band that month).
    errorbar_step : int, plot one error bar every N points (ignored if
        error_style='band').
    save : bool, saves the PNG in FIGURES_DIR if True (default).
    show : bool, calls plt.show() if True (default).
    """
    aws_funcs = {
        "z_surf_combined_anomaly": z_surf_combined_anomaly,
        "gps_alt_anomaly": gps_alt_anomaly,
        "snow_height": snow_height,
    }
    if aws_variable not in aws_funcs:
        raise KeyError(
            f"aws_variable={aws_variable!r} is unknown -- must be one of "
            f"{list(aws_funcs)}."
        )

    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
    N = len(STATIONS)
    NCOLS = 3
    NROWS = (N + NCOLS - 1) // NCOLS

    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()
    base_col = AWS_BASE_COLUMN.get(aws_variable, aws_variable)

    fig, axes = plt.subplots(
        NROWS,
        NCOLS,
        figsize=(14, NROWS * 3.5),
        constrained_layout=True,
    )
    axes_flat = np.atleast_1d(axes).ravel()

    for ax, (station, info) in zip(axes_flat, STATIONS.items()):
        # ── Per-station date window override (falls back to the global
        # start_date/end_date when the station isn't in dates_by_station) ──
        station_dates = (dates_by_station or {}).get(station)
        if station_dates is not None:
            start_bound_i = pd.to_datetime(station_dates["start"])
            end_bound_i = pd.to_datetime(station_dates["end"])
        else:
            start_bound_i = start_bound
            end_bound_i = end_bound

        _plot_aws_variable_on_ax(
            ax,
            station,
            info,
            aws_variable,
            aws_funcs,
            base_col,
            color,
            ref_date,
            start_bound_i,
            end_bound_i,
            aws_hourly_file_key,
            error_style,
            min_valid_fraction,
            errorbar_step,
        )

    for ax in axes_flat[N:]:
        ax.set_visible(False)

    fig.suptitle(fig_name, fontsize=15, fontweight="bold")

    if save:
        plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
        print(f"Figure saved: {OUTPUT_FILE}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, axes_flat


def plot_aws_variable_uncertainty_single(
    station,
    info,
    fig_name=None,
    aws_variable="z_surf_combined_anomaly",
    color="black",
    ref_date=None,
    start_date="2010-11-01",
    end_date=None,
    dates_by_station=None,
    aws_hourly_file_key="hourly_data",
    error_style="band",
    min_valid_fraction=0.5,
    errorbar_step=1,
    save=True,
    show=True,
    figsize=(8, 4.5),
):
    """One-station, one-figure version of plot_aws_variable_uncertainty --
    for when you want a separate, individually-sized PNG per station (e.g.
    to drop each one into a report) instead of one big subplot grid. Same
    drawing logic (curve, hourly-uncertainty band, Delta z_surf legend
    line) via the shared _plot_aws_variable_on_ax() helper, so a single
    station's panel looks identical whether it comes from this function or
    from the grid produced by plot_aws_variable_uncertainty.

    Call it in a loop, once per station -- note that STATIONS.items()
    already gives you both the station name AND its info dict, so pass
    both through:

        for station, info in AWS_data.STATIONS_ablation.items():
            plot_aws_variable_uncertainty_single(
                station, info,
                fig_name=f"z_surf_combined_anomaly with hourly uncertainty - {station}",
                dates_by_station=DATES,
                show=False,
            )

    station : the station's name/key (used for the title/legend and to look
        up dates_by_station -- it is NOT used to fetch the data, `info` is).
    info : that station's entry from a STATIONS_* dict, e.g.
        AWS_data.STATIONS_ablation[station] (must contain 'file' and, if
        error_style needs it, aws_hourly_file_key). This is why the loop
        must be `for station, info in STATIONS.items():` rather than just
        looping over the station names.
    fig_name : str or None, figure title and saved-PNG name (under
        FIGURES_DIR). Defaults to "{aws_variable} with hourly uncertainty
        - {station}" when None.
    figsize : matplotlib figsize for the single-panel figure.

    All other parameters (aws_variable, color, ref_date, start_date,
    end_date, dates_by_station, aws_hourly_file_key, error_style,
    min_valid_fraction, errorbar_step, save, show) have the exact same
    meaning as in plot_aws_variable_uncertainty. In particular
    dates_by_station still applies the per-station start/end override,
    keyed by `station` -- pass the same DATES dict you'd pass to the grid
    version.
    """
    aws_funcs = {
        "z_surf_combined_anomaly": z_surf_combined_anomaly,
        "gps_alt_anomaly": gps_alt_anomaly,
        "snow_height": snow_height,
    }
    if aws_variable not in aws_funcs:
        raise KeyError(
            f"aws_variable={aws_variable!r} is unknown -- must be one of "
            f"{list(aws_funcs)}."
        )

    if fig_name is None:
        fig_name = f"{aws_variable} with hourly uncertainty - {station}"
    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
    base_col = AWS_BASE_COLUMN.get(aws_variable, aws_variable)

    start_bound = pd.to_datetime(start_date) if start_date is not None else None
    end_bound = pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()
    station_dates = (dates_by_station or {}).get(station)
    if station_dates is not None:
        start_bound_i = pd.to_datetime(station_dates["start"])
        end_bound_i = pd.to_datetime(station_dates["end"])
    else:
        start_bound_i = start_bound
        end_bound_i = end_bound

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    _plot_aws_variable_on_ax(
        ax,
        station,
        info,
        aws_variable,
        aws_funcs,
        base_col,
        color,
        ref_date,
        start_bound_i,
        end_bound_i,
        aws_hourly_file_key,
        error_style,
        min_valid_fraction,
        errorbar_step,
    )

    fig.suptitle(fig_name, fontsize=13, fontweight="bold")

    if save:
        fig.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
        print(f"Figure saved: {OUTPUT_FILE}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, ax


# Overview grid: all stations on one figure.
# plot_aws_variable_uncertainty(
#     AWS_data.STATIONS_ablation,
#     "z_surf_combined_anomaly with hourly uncertainty",
#     dates_by_station=DATES,
# )

# # One individual figure per station (e.g. for dropping each into a report).
# for station, info in AWS_data.STATIONS_ablation.items():
#     plot_aws_variable_uncertainty_single(
#         station,
#         info,
#         fig_name=f"z_surf_combined_anomaly with hourly uncertainty - {station}",
#         dates_by_station=DATES,
#         show=False,
#     )


# for station, info in AWS_data.STATIONS_ablation.items():
#     plots_variables_with_satellite_uncertainty(
#         info,
#         list_satellites=['Nilsson and Gardner, 2026'],
#         fig_name=f"Surface Elevation Change according to  Nilsson and Garnder, 2026 - satellite uncertainty bands - {station}",
#         variables=['dh'],
#         colors=["#17EE17"],
#         uncertainty_variables=['rms'],
#     )


# plots_variables_with_satellite_uncertainty(
#     AWS_data.STATIONS_ablation,
#     'Surface Elevation Change according to  Nilsson and Garnder, 2026 - satellite uncertainty bands, ablation areas',
#     ['Nilsson and Gardner, 2026'], #'Copernicus_Climate_Data_Store', , 'Andersen et al., 2025', 'Khan et al., 2025', 'Zhang et al., 2022'
#     ['dh'], #'dh', , 'ZZ', 'dh_vol', 'elev_interp'
#     ["#17EE17"], #"#EB2F25", , "#EDE20D", "#D20DF0", "#03FEF1"
#     uncertainty_variables=['rms'], #'dh_uncert', , 'ZZer', None, 'elev_uncer_interp'
# )