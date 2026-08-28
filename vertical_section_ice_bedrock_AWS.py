"""
cross_section.py

Vertical cross-section of the Greenland Ice Sheet along the drift direction
of a single AWS station: the transect line runs from the station's first
valid recorded position to its last valid recorded position (station drift
over the ice, common for PROMICE/GC-Net stations), extended by EXTEND_KM on
each end (upstream of the first point, downstream of the last point).

The figure shows:
    - bed topography (BedMachine Greenland v6) with its uncertainty band (errbed)
    - reference surface (BedMachine, nominal_year)
    - mean multi-satellite surface elevation, one curve per hydrological year,
      averaged over Khan et al. 2025, Nilsson and Gardner 2026,
      Andersen et al. 2025 and Copernicus Climate Data Store
    - the station's first and last valid positions, each labelled with the
      station name + date above -- as two separate callouts when the two
      points are far enough apart, or as a single merged
      "STATION\nDATE1 -> DATE2" callout when they're too close together for
      two callouts to stay legible (near-stationary stations, tens to a
      couple hundred metres of net drift). (No dashed vertical ice column
      is drawn at these points anymore -- it was removed on request, since
      a per-point static-BedMachine thickness display, drawn at the two
      dates, risked being read as a change between them even though it
      wasn't one; see "Modelling choice" below.)
    - in the legend: the satellite-derived altitude difference between the
      two marked points (mean +/- inter-product std of the requested
      satellites at each point, propagated into the difference -- the two
      individual per-point altitudes are no longer shown on their own, only
      this difference), the reference-surface slope between the two points
      (how much BedMachine's static surface alone differs between the two
      positions -- the purely geometric part of the altitude difference,
      see "Modelling choice" below), the REAL, time-resolved ice-thickness
      difference from two individually-dated ArcticDEM strips (Delta h_ice,
      temporal -- see note below), and the straight-line horizontal
      (planimetric) displacement of the station between the two points.
      (There is no BedMachine-only ice-thickness-difference line: since
      BedMachine's "thickness" field has a single fixed epoch, no time
      dimension, differencing it between two spatial positions does not
      measure a real change over time -- see "Modelling choice" below --
      so that line was removed rather than left in a way that reads as a
      temporal result.)

Modelling choice :
    Satellite products only provide dh, an elevation change relative to a
    per-product reference epoch -- not an absolute surface elevation. To
    reconstruct an absolute surface usable in a cross-section, each product's
    dh is re-referenced to BED_REF_DATE (by default the BedMachine
    nominal_year) and added on top of the BedMachine reference surface:

        surface(t) = bed_surface_ref + [dh(t) - dh(BED_REF_DATE)]

    The four requested products are then averaged together for each
    hydrological year (Sept 1 of year Y-1 -> Aug 31 of year Y, labelled Y).

    A consequence worth spelling out: the altitude difference between the
    station's two marked points (Delta altitude in the legend) is NOT a
    pure temporal signal at a fixed location -- the station itself moved
    (its horizontal displacement, also in the legend) between the two
    dates, so the two points sit at two different places on BedMachine's
    static reference surface. Delta altitude therefore bundles (a) the
    real dh change measured by the satellites and (b) the purely
    geometric difference in the static reference surface between the two
    positions -- shown separately in the legend as "Reference-surface
    slope". Both z_surf_combined (from the AWS itself) and Delta altitude
    include this geometric term in the same way (both are true altitudes
    at wherever the station physically is), so comparing the two remains
    valid -- but it means the "dynamics" signal inferred from their
    difference isn't only vertical flux-divergence thickening, it also
    includes horizontal advection along the local slope (the u.grad(z_s)
    term in the classic dh/dt budget). The "Reference-surface slope" line
    lets a reader gauge how large that contribution is.

    Note on BedMachine's "thickness" field (surface - bed): it has a single
    ~nominal_year epoch, no time dimension, so differencing it between the
    two positions would only be a SPATIAL gradient (how much the
    2007-2008-era reference ice thickness varies over that stretch of the
    flowline), not a real change over the ~decade separating the two
    dates -- it would not tell you whether the ice column actually
    thickened or thinned between the two dates. For that reason this
    module does not compute or display that spatial difference at all.
    Instead, Delta h_ice (temporal) samples an absolute surface
    elevation from two individually-dated ArcticDEM strips (PGC STAC
    collection 'arcticdem-strips-s2s041-2m'), one near each date, and
    differences them against the bed (assumed constant through time --
    reasonable on a ~decade timescale, unlike the ice surface): see
    sample_arcticdem_strip_at_point() and
    real_delta_thickness_at_station()'s docstrings for the method and its
    main caveat (individual strips are not vertically co-registered here,
    so treat the result as indicative rather than a bias-corrected value).

Station-point altitude uncertainty (legend):
    At each of the two marked points, every requested satellite (per
    list_satellites) is independently sampled at that exact (x, y) location
    and at its own acquisition nearest to that point's date, and
    re-anchored to bed_ref_date the same way as everywhere else in this
    module (see "Modelling choice" above). The legend reports the mean and
    the inter-product standard deviation (ddof=1) of these per-satellite
    values -- the same "spread across datasets as uncertainty" convention
    used for the inter-dataset std band in plot_time_series_altimetry_AWS.py. A point where
    only one product has coverage shows its value with "unc. n/a" (a single
    number has no spread to measure); a point with no coverage at all shows
    "n/a". The altitude difference between the two points is end - start,
    with its uncertainty propagated as sqrt(unc_start^2 + unc_end^2)
    (standard error propagation for a difference of two independent
    quantities; missing sides count as 0 in the sum rather than voiding the
    whole result).

AWS relocations (see find_position_jumps, plot_cross_section's
relocation_handling parameter):
    A station's position record can contain a physical relocation, not
    just ice-flow drift -- PROMICE/GC-Net stations are periodically moved
    upstream when their current site is at risk of melting out near the
    ice margin. Treating the whole record as one continuous drift across
    such a jump builds a transect between two different physical
    locations, which can point the "wrong" way (e.g. the later fix
    landing closer to the accumulation zone than the earlier one, because
    the AWS was moved there -- not because the ice carried it backwards).
    By default (relocation_handling="longest_segment"), any jump larger
    than jump_threshold_m between two consecutive fixes is detected and
    printed as a warning, and the record is split into continuous
    deployment segments at each jump; the transect + the two marked points
    are then restricted to the longest segment (by time span) -- so a
    short recent relocation doesn't discard a longer, better-covered older
    deployment (e.g. a station moved once in 2020 after ten uninterrupted
    years, 2011-2020, keeps that whole decade rather than being cut down
    to 2020-present). Pass relocation_handling="latest_segment" to always
    use the station's current site instead (even if shorter), or
    relocation_handling="full_record" to ignore relocations entirely (not
    recommended for interpretation, but available for comparison/
    debugging).

Adjust BED_TOPOGRAPHY_DATA_FILE below to the actual path of the BedMachine NetCDF file.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from matplotlib.lines import Line2D
from pystac_client import Client
from pyproj import Transformer
from rasterio.windows import Window

import AWS_data
import interpolation_altimetry_AWS
from paths import BED_TOPOGRAPHY_DATA_FILE, VELOCITY_DATA_DIR, FIGURES_DIR as _BASE_FIGURES_DIR


# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Station coordinates ──────────────────────────────────────────────────
"""
def find_position_jumps(station_file, jump_threshold_m) -- flags abrupt jumps between consecutive valid position fixes, almost always an AWS relocation rather than ice-flow drift.

def find_continuous_segments(station_file, jump_threshold_m) -- splits a station's position record into continuous deployment segments at each detected relocation jump.

def longest_continuous_segment(station_file, jump_threshold_m) -- returns the start/end of the longest continuous deployment segment in a station's record.

def latest_continuous_segment_start(station_file, jump_threshold_m) -- returns the timestamp of the first fix right after the station's most recent detected relocation.

def get_station_first_last(station_file, min_time, max_time) -- returns the station's first and last valid recorded position, optionally restricted to a time window.

def get_station_position_near(station_file, target_date, period_start, period_end) -- returns the station's valid position closest to a given date, optionally restricted to a period.

def build_transect_from_station(station_file, extend_km, n_points, min_time, max_time) -- builds the transect line running through a station's drift direction, extended past its first/last positions.

def project_km(x_pt, y_pt, x_start, y_start, ux, uy) -- projects an arbitrary point's distance (km) along the transect direction, measured from its start.
"""

# ── Transect geometry ────────────────────────────────────────────────────
"""
def build_transect(xy_start, xy_end, n_points) -- builds a straight-line transect between two points and returns the sampled coordinates with their cumulative distance.

def sample_along_line(da, x_line, y_line) -- samples a 2D gridded array at the transect points, nearest-neighbour, lazily on dask-backed data.
"""

# ── Transect orientation ─────────────────────────────────────────────────
"""
def compute_bearing_deg(lon1, lat1, lon2, lat2) -- computes the true geographic bearing (0-360°) from one point to another.

def bearing_to_compass(bearing_deg) -- converts a bearing in degrees to a 16-point compass label (e.g. "NNE").
"""

# ── Hydrological years ───────────────────────────────────────────────────
"""
def hydrological_year(timestamp) -- returns the hydrological year label (Sept 1 -> Aug 31) for a given timestamp.
"""

# ── Ice surface velocity (Sentinel-1 offset-tracking, one file per date) ──
"""
def list_velocity_files(velocity_dir=VELOCITY_DATA_DIR, pattern="*.nc") -- lists every velocity NetCDF file found in a directory.

def get_velocity_time(path) -- reads a velocity file's acquisition date from its own time coordinate.

def sample_velocity_magnitude(path, x_line, y_line, variable="land_ice_surface_velocity_magnitude") -- samples the ice velocity magnitude of a single file along the transect.

def build_velocity_by_year(
    x_line, y_line, years, velocity_dir=VELOCITY_DATA_DIR, pattern="*.nc",
    variable="land_ice_surface_velocity_magnitude",
) -- averages the velocity magnitude along the transect over every acquisition of each requested hydrological year.

"""

# ── Point altitude from satellites, at the station's drift endpoints ──────
"""
def satellite_surface_at_point(list_satellites, x_pt, y_pt, t_target, ref_surface_val, bed_ref_ts) -- reconstructs the absolute surface elevation at a single point from each requested satellite's dh.

def mean_std_across_products(values, ddof=1) -- computes the mean and inter-product standard deviation of a list of per-satellite values.

def format_altitude_legend(label, alt, unc, n_products) -- formats one legend line showing an altitude with its uncertainty, gracefully degraded when data is missing.
"""

# ── Real temporal ice-thickness change (ArcticDEM strips) ─────────────────
"""
def sample_arcticdem_strip_at_point(x_pt, y_pt, target_date, search_window_days, collection) -- finds the individually-dated ArcticDEM strip closest in time to target_date covering a point, and returns the raw elevation sampled there.

def real_delta_thickness_at_station(x_first, y_first, t_first, x_last, y_last, t_last, bed_first, bed_last, search_window_days) -- real, time-resolved ice-thickness difference at the station's start/end position, from two ArcticDEM strips minus the (assumed static) BedMachine bed.

def format_real_delta_thickness_legend(result) -- formats one legend line summarizing real_delta_thickness_at_station's result, gracefully degraded when no strip is found.
"""

# ── Main figure ───────────────────────────────────────────────────────────
"""
def plot_cross_section(...) -- plots the full vertical cross-section figure for one station: bed topography, reference surface, yearly multi-satellite surface curves, the two marked drift endpoints, and an optional velocity panel. Internal steps:
    - Station drift transect
    - Bed topography (BedMachine)
    - Per-satellite dh sampled along the line, re-anchored to bed_ref_date
    - Determine which hydrological years are common to all satellites
    - Mean multi-satellite surface per hydrological year
    - Station position at the start/end of the plotted period
    - Satellite-derived altitude at the two marked points
    - Ice surface velocity, mean per hydrological year
    - Plot
"""

# ── Batch plotting over every station of a STATIONS_* dict ────────────────
"""
def plot_all_stations(stations_dict, stations, stop_on_error, **kwargs) -- calls plot_cross_section once per station, skipping (by default) any station whose plot raises an error.
"""

# ── Paths ─────────────────────────────────────────────────────────────────
"""
FIGURES_DIR, DEFAULT_SATELLITES, _LATLON_TO_3413, _XY3413_TO_LATLON, ARCTICDEM_STAC_URL, ARCTICDEM_STRIPS_COLLECTION -- output folder, default satellite set used for the multi-satellite mean, the shared lon/lat <-> EPSG:3413 coordinate transformers (both directions), and the PGC STAC endpoint/collection used to fetch individually-dated ArcticDEM strips (see sample_arcticdem_strip_at_point below).
"""


FIGURES_DIR = _BASE_FIGURES_DIR / "cross_sections"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Default set of products used for the multi-satellite mean (Zhang et al.
# 2022 intentionally excluded here, per request).
DEFAULT_SATELLITES = [
    "Khan et al., 2025",
    "Nilsson and Gardner, 2026",
    "Andersen et al., 2025",
    "Copernicus_Climate_Data_Store",
]

_LATLON_TO_3413 = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
_XY3413_TO_LATLON = Transformer.from_crs("EPSG:3413", "EPSG:4326", always_xy=True)

# PGC STAC catalog, same endpoint used for the ArcticDEM *mosaic* in
# plot_maps.py -- here queried for the *strips* collection instead (each
# strip has its own genuine acquisition date, unlike the mosaic which
# blends multiple epochs together and can't be used for real temporal
# differencing). See sample_arcticdem_strip_at_point().
ARCTICDEM_STAC_URL = "https://stac.pgc.umn.edu/api/v1/"
ARCTICDEM_STRIPS_COLLECTION = "arcticdem-strips-s2s041-2m"


# ── Station coordinates ──────────────────────────────────────────────────
def find_position_jumps(station_file, jump_threshold_m=150.0):
    """Detect abrupt jumps between consecutive valid (lat, lon) fixes that
    are too large to be explained by continuous ice-flow drift -- almost
    certainly a physical relocation of the AWS itself (PROMICE/GC-Net
    stations are periodically moved upstream when the current site is at
    risk of melting out near the ice margin), not a sudden acceleration of
    the ice. build_transect_from_station and plot_cross_section otherwise
    silently treat the *entire* position record as one continuous drift --
    across an undetected relocation that produces a transect built between
    two physically different locations, which can even point the "wrong"
    way (e.g. a later fix landing closer to the accumulation zone than an
    earlier one, because the AWS was moved there, not because the ice
    carried it backwards).

    jump_threshold_m : any consecutive-fix displacement above this is
        flagged. 150 m is deliberately far above plausible drift between
        two fixes (even a fast-flowing outlet-glacier AWS at a very high
        ~5 m/day only covers ~150 m in a month, and most AWS/GC-Net sites
        are much slower than that -- see the velocity panel), and well
        below a real relocation, which is normally several hundred metres
        to a few kilometres. This check works regardless of the file's
        sampling cadence: a real relocation jumps within a single
        sampling interval no matter how fine that interval is, while
        genuine drift per interval only gets smaller as the cadence gets
        finer.

    Returns a list of dicts (empty if no jump found), one per detected
    jump, each with "time_before", "time_after" (the two straddling
    fixes) and "jump_m" (their straight-line distance in metres).
    """
    aws = pd.read_csv(station_file)
    aws["time"] = pd.to_datetime(aws["time"])
    valid = aws.dropna(subset=["lat", "lon"]).sort_values("time")
    if len(valid) < 2:
        return []

    x, y = _LATLON_TO_3413.transform(valid["lon"].to_numpy(), valid["lat"].to_numpy())
    times = valid["time"].to_numpy()

    jumps = []
    for i in range(1, len(valid)):
        d = float(np.hypot(x[i] - x[i - 1], y[i] - y[i - 1]))
        if d > jump_threshold_m:
            jumps.append({
                "time_before": pd.Timestamp(times[i - 1]),
                "time_after": pd.Timestamp(times[i]),
                "jump_m": d,
            })
    return jumps


def find_continuous_segments(station_file, jump_threshold_m=150.0):
    """Split a station's position record into continuous deployment
    segments at each relocation jump detected by find_position_jumps.

    Returns a list of dicts, ordered chronologically, one per segment,
    each with "start", "end" (timestamps of the segment's first and last
    valid fix) and "duration_days" ((end - start).days). A record with no
    detected jump returns a single segment spanning the whole record; an
    empty/single-fix record returns an empty list.
    """
    aws = pd.read_csv(station_file)
    aws["time"] = pd.to_datetime(aws["time"])
    valid = aws.dropna(subset=["lat", "lon"]).sort_values("time")
    if valid.empty:
        return []

    times = [pd.Timestamp(t) for t in valid["time"]]
    jump_after_times = {
        j["time_after"] for j in find_position_jumps(station_file, jump_threshold_m=jump_threshold_m)
    }

    segments = []
    seg_start_prev = times[0]
    for i in range(1, len(times)):
        if times[i] in jump_after_times:
            seg_end = times[i - 1]
            segments.append({
                "start": seg_start_prev, "end": seg_end,
                "duration_days": (seg_end - seg_start_prev).days,
            })
            seg_start_prev = times[i]
    segments.append({
        "start": seg_start_prev, "end": times[-1],
        "duration_days": (times[-1] - seg_start_prev).days,
    })
    return segments


def longest_continuous_segment(station_file, jump_threshold_m=150.0):
    """(segment_start, segment_end) of the longest continuous deployment
    segment (by time span) in a station's position record, split at each
    relocation jump detected by find_position_jumps (see
    find_continuous_segments). A record with no detected jump returns the
    whole record's (first, last) fix, same as not restricting at all. Ties
    are broken in favour of the more recent segment. Returns (None, None)
    for an empty/single-fix record.

    Used by plot_cross_section (relocation_handling="longest_segment", the
    default) so a short, recent relocation doesn't discard a longer,
    better-covered earlier deployment -- e.g. a station moved once in 2020
    after ten uninterrupted years (2011-2020) keeps that whole decade
    instead of being cut down to just 2020-present.
    """
    segments = find_continuous_segments(station_file, jump_threshold_m=jump_threshold_m)
    if not segments:
        return None, None
    longest = max(segments, key=lambda s: (s["duration_days"], s["start"]))
    return longest["start"], longest["end"]


def latest_continuous_segment_start(station_file, jump_threshold_m=150.0):
    """Timestamp of the first valid fix *after* the most recent relocation
    jump detected by find_position_jumps, or None if no jump was found
    (the whole record is then treated as one continuous deployment). Used
    by plot_cross_section (relocation_handling="latest_segment") to
    restrict the transect and marked points to the station's *current*
    site, even if a longer-but-older segment exists -- see
    longest_continuous_segment for the default behaviour instead.
    """
    jumps = find_position_jumps(station_file, jump_threshold_m=jump_threshold_m)
    if not jumps:
        return None
    return jumps[-1]["time_after"]


def get_station_first_last(station_file, min_time=None, max_time=None):
    """First and last valid (time, lon, lat) recorded for a station,
    sorted chronologically. Rows with missing lat/lon are dropped first.

    min_time, max_time : if given, fixes outside [min_time, max_time] are
        excluded before picking first/last -- used to restrict to a single
        continuous deployment segment, after a detected relocation (see
        find_position_jumps / longest_continuous_segment /
        latest_continuous_segment_start)."""
    aws = pd.read_csv(station_file)
    aws["time"] = pd.to_datetime(aws["time"])
    valid = aws.dropna(subset=["lat", "lon"]).sort_values("time")
    if min_time is not None:
        valid = valid[valid["time"] >= min_time]
    if max_time is not None:
        valid = valid[valid["time"] <= max_time]
    if valid.empty:
        bounds = []
        if min_time is not None:
            bounds.append(f"on/after {pd.Timestamp(min_time).date()}")
        if max_time is not None:
            bounds.append(f"on/before {pd.Timestamp(max_time).date()}")
        suffix = f" {' and '.join(bounds)}" if bounds else ""
        raise ValueError(f"No valid lat/lon found in {station_file}{suffix}.")
    first = valid.iloc[0]
    last = valid.iloc[-1]
    return (
        (first["time"], first["lon"], first["lat"]),
        (last["time"], last["lon"], last["lat"]),
    )


def get_station_position_near(station_file, target_date, period_start=None, period_end=None):
    """Station's valid (time, lon, lat) closest to target_date. Used to mark
    the station's position at the start/end of the period actually covered
    by the plotted hydrological years, rather than the absolute first/last
    valid fix of the whole record (which can fall outside that period).

    If period_start/period_end are given, candidate rows are restricted to
    that window before searching for the nearest one -- this both enforces
    "not before nor after" the plotted period, and protects against a
    spurious/erroneous timestamp far outside the real record (e.g. a stray
    row dated decades before AWS deployment) from being picked as
    "nearest" simply because of a large data gap elsewhere.
    """
    aws = pd.read_csv(station_file)
    aws["time"] = pd.to_datetime(aws["time"])
    valid = aws.dropna(subset=["lat", "lon"])
    if period_start is not None:
        valid = valid[valid["time"] >= period_start]
    if period_end is not None:
        valid = valid[valid["time"] <= period_end]
    if valid.empty:
        raise ValueError(
            f"No valid lat/lon found in {station_file} "
            f"within [{period_start}, {period_end}]."
        )
    idx = (valid["time"] - target_date).abs().idxmin()
    row = valid.loc[idx]
    return row["time"], row["lon"], row["lat"]


def build_transect_from_station(station_file, extend_km=5, n_points=500, min_time=None, max_time=None):
    """Transect running through a station's drift direction: from its first
    valid position to its last valid position, extended by extend_km on
    each end (upstream of the first point, downstream of the last point).

    min_time, max_time : forwarded to get_station_first_last -- restricts
        "first"/"last" to fixes within [min_time, max_time], used to build
        the transect from a single continuous deployment segment only,
        after a detected relocation (see find_position_jumps /
        longest_continuous_segment in plot_cross_section).

    Returns a dict with:
        x_line, y_line, distance_km : the sampled transect (see build_transect)
        dist_first_km, dist_last_km : position of the station's first/last
            valid point along distance_km
        t_first, t_last : timestamps of the first/last valid point
    """
    (t_first, lon_first, lat_first), (t_last, lon_last, lat_last) = (
        get_station_first_last(station_file, min_time=min_time, max_time=max_time)
    )
    x_first, y_first = _LATLON_TO_3413.transform(lon_first, lat_first)
    x_last, y_last = _LATLON_TO_3413.transform(lon_last, lat_last)

    dx, dy = x_last - x_first, y_last - y_first
    length_m = np.hypot(dx, dy)
    if length_m == 0:
        raise ValueError(
            f"Station has an identical first and last valid position "
            f"({station_file}) -- cannot define a transect direction."
        )
    ux, uy = dx / length_m, dy / length_m

    extend_m = extend_km * 1000.0
    xy_start = (x_first - ux * extend_m, y_first - uy * extend_m)
    xy_end = (x_last + ux * extend_m, y_last + uy * extend_m)

    x_line, y_line, distance_km = build_transect(xy_start, xy_end, n_points=n_points)
    bearing_deg = compute_bearing_deg(lon_first, lat_first, lon_last, lat_last)

    return {
        "x_line": x_line,
        "y_line": y_line,
        "distance_km": distance_km,
        "dist_first_km": extend_km,
        "dist_last_km": extend_km + length_m / 1000.0,
        "t_first": t_first,
        "t_last": t_last,
        "bearing_deg": bearing_deg,
        "lon_first": lon_first,
        "lat_first": lat_first,
        "lon_last": lon_last,
        "lat_last": lat_last,
        "x_start": xy_start[0],
        "y_start": xy_start[1],
        "ux": ux,
        "uy": uy,
    }


def project_km(x_pt, y_pt, x_start, y_start, ux, uy):
    """Distance (km) of an arbitrary point along the transect direction,
    projected from the transect's (extended) start point."""
    return ((x_pt - x_start) * ux + (y_pt - y_start) * uy) / 1000.0


# ── Transect geometry ────────────────────────────────────────────────────
def build_transect(xy_start, xy_end, n_points=500):
    """Straight-line transect between two EPSG:3413 points.

    Returns (x_line, y_line, distance_km), all length n_points, with
    distance_km the cumulative distance from xy_start (0 at the start
    station, positive towards the end station).
    """
    x0, y0 = xy_start
    x1, y1 = xy_end
    t = np.linspace(0, 1, n_points)
    x_line = x0 + t * (x1 - x0)
    y_line = y0 + t * (y1 - y0)
    total_km = np.hypot(x1 - x0, y1 - y0) / 1000.0
    distance_km = t * total_km
    return x_line, y_line, distance_km


def sample_along_line(da, x_line, y_line):
    """Nearest-neighbour vectorized sampling of a 2D (y, x) DataArray at the
    transect points. Works lazily on dask-backed arrays: only the selected
    points are actually loaded into memory."""
    x_idx = xr.DataArray(x_line, dims="points")
    y_idx = xr.DataArray(y_line, dims="points")
    sampled = da.sel(x=x_idx, y=y_idx, method="nearest")
    return np.asarray(sampled.values, dtype="float64")


# ── Transect orientation ─────────────────────────────────────────────────
def compute_bearing_deg(lon1, lat1, lon2, lat2):
    """True geographic bearing (degrees clockwise from North, 0-360) from
    point 1 to point 2, computed from lat/lon (not from the projected
    EPSG:3413 coordinates, which distort angles this close to the pole)."""
    lat1_r, lat2_r = np.radians(lat1), np.radians(lat2)
    dlon_r = np.radians(lon2 - lon1)
    x = np.sin(dlon_r) * np.cos(lat2_r)
    y = np.cos(lat1_r) * np.sin(lat2_r) - np.sin(lat1_r) * np.cos(lat2_r) * np.cos(dlon_r)
    return (np.degrees(np.arctan2(x, y)) + 360) % 360


def bearing_to_compass(bearing_deg):
    """16-point compass label for a bearing in degrees."""
    directions = [
        "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
        "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW",
    ]
    idx = int((bearing_deg + 11.25) // 22.5) % 16
    return directions[idx]


# ── Hydrological years ───────────────────────────────────────────────────
def hydrological_year(timestamp):
    """Hydrological year label for a timestamp: Sept 1 (Y-1) -> Aug 31 (Y)
    is labelled Y. Adjust here if your convention differs (e.g. calendar
    year, or a different start month)."""
    ts = pd.Timestamp(timestamp)
    return ts.year + 1 if ts.month >= 9 else ts.year


# ── Ice surface velocity (Sentinel-1 offset-tracking, one file per date) ──
def list_velocity_files(velocity_dir=VELOCITY_DATA_DIR, pattern="*.nc"):
    """All velocity NetCDF files found in velocity_dir (non-recursive)."""
    files = sorted(Path(velocity_dir).glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No velocity file found in {velocity_dir} (pattern={pattern!r})."
        )
    return files


def get_velocity_time(path):
    """Acquisition date of a velocity file, read from its own 'time'
    coordinate (more robust than parsing the filename)."""
    with xr.open_dataset(path, chunks={}) as ds:
        return pd.Timestamp(np.asarray(ds["time"].values).ravel()[0])


def sample_velocity_magnitude(path, x_line, y_line, variable="land_ice_surface_velocity_magnitude"):
    """Sample the velocity magnitude of a single velocity file along the
    transect."""
    with xr.open_dataset(path, chunks={}) as ds:
        da = ds[variable]
        if "time" in da.dims:
            da = da.isel(time=0)
        return sample_along_line(da, x_line, y_line)


def build_velocity_by_year(
    x_line, y_line, years, velocity_dir=VELOCITY_DATA_DIR, pattern="*.nc",
    variable="land_ice_surface_velocity_magnitude",
):
    """For each requested hydrological year, average the velocity magnitude
    (sampled along the transect) over all acquisitions falling in that year.
    Years without any acquisition are simply absent from the returned dict.

    Returns (velocity_by_year, units) where units is read from the
    variable's `units` attribute of the first file opened (empty string if
    absent).
    """
    files = list_velocity_files(velocity_dir, pattern)
    file_times = [(f, get_velocity_time(f)) for f in files]

    units = ""
    with xr.open_dataset(file_times[0][0], chunks={}) as ds0:
        units = ds0[variable].attrs.get("units", "")

    velocity_by_year = {}
    for year in years:
        files_this_year = [f for f, t in file_times if hydrological_year(t) == year]
        if not files_this_year:
            continue
        samples = [
            sample_velocity_magnitude(f, x_line, y_line, variable=variable)
            for f in files_this_year
        ]
        velocity_by_year[year] = np.nanmean(np.vstack(samples), axis=0)

    return velocity_by_year, units


# ── Point altitude from satellites, at the station's drift endpoints ──────
def satellite_surface_at_point(list_satellites, x_pt, y_pt, t_target, ref_surface_val, bed_ref_ts):
    """For each requested satellite, sample its dh at the single point
    (x_pt, y_pt) -- not along the transect line -- at the acquisition
    nearest to t_target and at the acquisition nearest to bed_ref_ts, then
    reconstruct an absolute surface elevation from that product alone:

        surface_this_sat = ref_surface_val + [dh(t_target) - dh(bed_ref_ts)]

    (same re-anchoring convention as the rest of the module, see the module
    docstring). A product with no data at that exact point (all-NaN pixel,
    e.g. outside its coverage or masked) is silently skipped.

    Returns (values, used): values is the list of per-satellite surface
    elevations (float) that were actually computed, used is the matching
    list of satellite names -- so the caller knows which/how many products
    contributed.
    """
    values, used = [], []
    t_target = pd.to_datetime(t_target)
    for sat_name in list_satellites:
        info_sat = interpolation_altimetry_AWS.SATELLITE[sat_name]
        ds = interpolation_altimetry_AWS.satellite_opening(info_sat["file"])
        var = ds[info_sat["var"]]

        point = var.sel(x=x_pt, y=y_pt, method="nearest")
        dh_ref_val = float(point.sel(time=bed_ref_ts, method="nearest").values)
        dh_target_val = float(point.sel(time=t_target, method="nearest").values)

        if np.isnan(dh_ref_val) or np.isnan(dh_target_val):
            continue

        values.append(ref_surface_val + (dh_target_val - dh_ref_val))
        used.append(sat_name)

    return values, used


def mean_std_across_products(values, ddof=1):
    """Mean and inter-product standard deviation of a list of per-satellite
    values -- same convention as _multi_satellite_std_band in
    plot_time_series_altimetry_AWS.py (spread across products used as the uncertainty
    estimate). Returns (mean, std): std is NaN if fewer than 2 values are
    available (a single product has no spread to measure), and both are NaN
    if the list is empty."""
    if len(values) == 0:
        return np.nan, np.nan
    arr = np.asarray(values, dtype="float64")
    if len(arr) == 1:
        return float(arr[0]), np.nan
    return float(np.mean(arr)), float(np.std(arr, ddof=ddof))


def format_altitude_legend(label, alt, unc, n_products):
    """One legend line: 'label: alt ± unc m', gracefully degraded when
    altitude or uncertainty can't be computed (no product coverage at that
    point/date, or only one product available so no inter-product spread)."""
    if np.isnan(alt):
        return f"{label}: n/a (no satellite coverage at that point/date)"
    if np.isnan(unc):
        return f"{label}: {alt:.1f} m (unc. n/a, {n_products} product)"
    return f"{label}: {alt:.1f} ± {unc:.1f} m ({n_products} products)"


# ── Real temporal ice-thickness change (ArcticDEM strips) ─────────────────
def sample_arcticdem_strip_at_point(
    x_pt, y_pt, target_date, search_window_days=180,
    collection=ARCTICDEM_STRIPS_COLLECTION,
):
    """Finds, among ArcticDEM *strips* (PGC STAC collection
    'arcticdem-strips-s2s041-2m' -- individually-dated 2 m stereo-derived
    DEMs, NOT the static multi-epoch mosaic used elsewhere in this project
    for the hillshade background in plot_maps.py) covering the point
    (x_pt, y_pt) in EPSG:3413, the one whose acquisition date is closest to
    target_date, within +/- search_window_days. Widen search_window_days if
    no strip is found near a given date/location -- strip coverage is
    opportunistic (tied to commercial stereo acquisitions), not systematic
    like the altimetry products used elsewhere in this module.

    Returns (elevation_m, strip_date, strip_id): the strip's raw elevation
    sampled at that exact point (nearest pixel), its actual acquisition
    date (which will generally NOT fall exactly on target_date -- that's
    why it's returned, so the caller can report the real date used), and
    the STAC item id (for traceability, e.g. to go check that strip's own
    'readme'/'metadata' asset if a result looks suspicious). Returns
    (np.nan, None, None) if no strip is found in the search window.

    IMPORTANT -- no vertical co-registration is applied here. Unlike the
    ArcticDEM *mosaic*, individual strips are NOT guaranteed to share a
    common absolute vertical datum: each one carries its own systematic
    bias from its own bundle adjustment (commonly on the order of a metre,
    sometimes several), on top of the ~0.1-0.3 m relative accuracy within
    a strip itself. A rigorous analysis should co-register each strip
    first -- e.g. against ICESat-2 ATL06 ground tracks, or via a
    stable-terrain slope/aspect regression (Nuth & Kaab, 2011) against the
    ArcticDEM mosaic nearby -- before differencing two strips. That step is
    NOT implemented here: treat the resulting Delta h_ice as indicative
    only, not a bias-corrected measurement.
    """
    target_date = pd.Timestamp(target_date)
    window = pd.Timedelta(days=search_window_days)

    # Small bbox around the point (STAC bbox search needs a box, not a
    # point) -- 500 m on each side is comfortably inside a single strip
    # footprint (strips are several km wide) while keeping the search tight.
    buffer_m = 500.0
    lon_min, lat_min = _XY3413_TO_LATLON.transform(x_pt - buffer_m, y_pt - buffer_m)
    lon_max, lat_max = _XY3413_TO_LATLON.transform(x_pt + buffer_m, y_pt + buffer_m)
    bbox = [
        min(lon_min, lon_max), min(lat_min, lat_max),
        max(lon_min, lon_max), max(lat_min, lat_max),
    ]

    catalog = Client.open(ARCTICDEM_STAC_URL)
    search = catalog.search(
        collections=[collection],
        bbox=bbox,
        datetime=f"{(target_date - window).date()}/{(target_date + window).date()}",
    )
    items = list(search.items())
    if not items:
        print(
            f"[sample_arcticdem_strip_at_point] no ArcticDEM strip found "
            f"within +/-{search_window_days} d of {target_date.date()} at "
            f"({x_pt:.0f}, {y_pt:.0f}) -- try widening search_window_days."
        )
        return np.nan, None, None

    def _item_date(item):
        return pd.Timestamp(item.datetime).tz_localize(None)

    best_item = min(items, key=lambda it: abs(_item_date(it) - target_date))
    strip_date = _item_date(best_item)

    asset = best_item.assets.get("dem")
    if asset is None:
        raise ValueError(
            f"Strip {best_item.id} has no 'dem' asset -- unexpected for "
            f"collection {collection!r}."
        )

    with rasterio.open(asset.href) as src:
        row, col = src.index(x_pt, y_pt)
        window_1px = Window(col, row, 1, 1)
        value = src.read(1, window=window_1px)[0, 0]
        if src.nodata is not None and value == src.nodata:
            value = np.nan

    return float(value), strip_date, best_item.id


def real_delta_thickness_at_station(
    x_first, y_first, t_first, x_last, y_last, t_last,
    bed_first, bed_last, search_window_days=180,
):
    """Real, time-resolved ice-thickness difference under the station
    between its start and end position/date. This is deliberately NOT
    computed from BedMachine's static 'thickness' field (surface - bed):
    since BedMachine has a single reference epoch and no time dimension,
    differencing it between the station's two (different) positions would
    only be the SPATIAL gradient of that static field, not a real change
    over time -- see the module docstring's "Modelling choice" note.

    Here, instead, an absolute surface elevation is sampled at each
    position from an individually-dated ArcticDEM strip closest in time to
    that position's own date (via sample_arcticdem_strip_at_point), and the
    bed is assumed constant through time (BedMachine's own convention --
    bed topography changes on glacial/geological timescales, far slower
    than the ~decade separating t_first and t_last):

        thickness(t) ~= arcticdem_surface(t) - bed   (bed assumed static)

    Returns a dict with: delta_thickness_m (thickness_last - thickness_first,
    NaN if either strip is missing), thickness_first_m, thickness_last_m,
    strip_date_first, strip_date_last (the strips' actual acquisition
    dates -- generally not exactly t_first/t_last), strip_id_first,
    strip_id_last. See sample_arcticdem_strip_at_point's docstring for the
    important caveat that individual strips are not vertically
    co-registered here -- this result should be treated as indicative.
    """
    elev_first, strip_date_first, strip_id_first = sample_arcticdem_strip_at_point(
        x_first, y_first, t_first, search_window_days=search_window_days,
    )
    elev_last, strip_date_last, strip_id_last = sample_arcticdem_strip_at_point(
        x_last, y_last, t_last, search_window_days=search_window_days,
    )

    thickness_first = elev_first - bed_first if not np.isnan(elev_first) else np.nan
    thickness_last = elev_last - bed_last if not np.isnan(elev_last) else np.nan
    delta_thickness = (
        thickness_last - thickness_first
        if not (np.isnan(thickness_first) or np.isnan(thickness_last))
        else np.nan
    )

    return {
        "delta_thickness_m": delta_thickness,
        "thickness_first_m": thickness_first,
        "thickness_last_m": thickness_last,
        "strip_date_first": strip_date_first,
        "strip_date_last": strip_date_last,
        "strip_id_first": strip_id_first,
        "strip_id_last": strip_id_last,
    }


def format_real_delta_thickness_legend(result):
    """One legend line summarizing real_delta_thickness_at_station's
    result, gracefully degraded to 'n/a' with a short reason if either
    ArcticDEM strip is missing, or to a distinct message if the
    computation was skipped altogether (result["skipped_reason"] set --
    see plot_cross_section's compute_real_thickness parameter)."""
    if result.get("skipped_reason"):
        return f"Δh_ice, temporal: not computed ({result['skipped_reason']})"
    if np.isnan(result["delta_thickness_m"]):
        return "Δh_ice, temporal (ArcticDEM strips): n/a (no strip found near start and/or end)"
    return (
        f"Δh_ice, temporal (ArcticDEM strips, "
        f"{result['strip_date_first'].date()} → {result['strip_date_last'].date()}): "
        f"{result['delta_thickness_m']:+.1f} m (not co-registered -- indicative only)"
    )


# ── Main figure ───────────────────────────────────────────────────────────
def plot_cross_section(
    station="KAN_L",
    stations_dict=AWS_data.STATIONS_ablation,
    extend_km=5,
    list_satellites=None,
    n_points=500,
    bed_ref_date=None,
    start_year=None,
    end_year=None,
    cmap_name="viridis",
    include_velocity=True,
    velocity_dir=VELOCITY_DATA_DIR,
    velocity_pattern="*.nc",
    relocation_handling="longest_segment",
    jump_threshold_m=150.0,
    compute_real_thickness=True,
    real_thickness_search_window_days=180,
    save=True,
    output_path=None,
    show=True,
):
    """
    Plot a vertical cross-section of the ice sheet along the drift direction
    of a single AWS station (first valid position -> last valid position),
    extended by extend_km on each end.

    station : key in stations_dict (e.g. AWS_data.STATIONS_ablation)
    extend_km : distance (km) added beyond the station's first and last
        valid positions, on each end of the transect.
    list_satellites : list of satellite names (keys of interpolation_altimetry_AWS.SATELLITE);
        defaults to DEFAULT_SATELLITES (Khan, Nilsson, Andersen, Copernicus).
    n_points : number of points sampled along the straight transect line.
    bed_ref_date : reference date used to re-anchor each satellite's dh onto
        the BedMachine reference surface. Defaults to the BedMachine
        `nominal_year` attribute (Jan 1st of that year) if present in the
        NetCDF metadata, else "2008-01-01".
    start_year, end_year : restrict the hydrological years plotted (inclusive).
        If None, all years with data for every requested satellite are used.
    include_velocity : bool, if True adds a second panel with ice surface
        velocity magnitude (mean per hydrological year, from Sentinel-1
        offset-tracking files in velocity_dir). Years without any
        acquisition are simply skipped in that panel.
    velocity_dir, velocity_pattern : directory and glob pattern used to
        find the per-date velocity NetCDF files.
    relocation_handling : how to handle a detected AWS relocation (an
        abrupt jump in the station's recorded position -- see
        find_position_jumps -- almost always a physical move of the
        instrument to a new site, not the ice suddenly carrying it
        backwards). One of:
            "longest_segment" (default) -- split the record into
                continuous deployment segments at every detected jump,
                and use the longest one (by time span). A short recent
                relocation doesn't discard a longer, better-covered older
                deployment this way -- e.g. a station moved once in 2020
                after ten uninterrupted years (2011-2020) keeps that whole
                decade instead of being cut down to 2020-present.
            "latest_segment" -- always restrict to fixes on/after the
                most recent detected jump, i.e. the station's *current*
                site, even if that segment is shorter than an earlier one.
            "full_record" -- ignore any detected jump and use the whole
                position record as before (kept for comparison/debugging;
                not recommended for interpretation: a transect spanning
                two physically different locations doesn't describe a
                single flowline, and can even point the "wrong" way if a
                relocation moved the AWS upstream).
        Either way, any detected jump is printed as a warning so it's
        never silently hidden.
    jump_threshold_m : forwarded to find_position_jumps -- the minimum
        consecutive-fix displacement (metres) treated as a relocation
        rather than ice-flow drift.
    compute_real_thickness : bool, if True (default) computes the real,
        time-resolved ice-thickness change under the station
        (real_delta_thickness_at_station -- two ArcticDEM strips minus the
        bed) and shows it in the legend. Set False to skip it entirely --
        no PGC STAC query, no strip download -- e.g. for a quick/offline
        run, a batch call over many stations where the network round-trip
        adds up, or simply when that line isn't needed. When False, the
        legend shows "Δh_ice, temporal: not computed (compute_real_thickness=False)"
        instead.
    real_thickness_search_window_days : forwarded to
        real_delta_thickness_at_station (then to
        sample_arcticdem_strip_at_point) as search_window_days -- how many
        days on each side of the station's start/end date to search for the
        closest ArcticDEM strip. Ignored when compute_real_thickness=False.
    show : if True (default), display the figure with plt.show() (blocks
        until the window is closed). Set False for batch use (e.g. from
        plot_all_stations()) so the loop doesn't stall on each figure --
        the figure is then closed right after saving instead.
    """
    if list_satellites is None:
        list_satellites = DEFAULT_SATELLITES
    if relocation_handling not in ("longest_segment", "latest_segment", "full_record"):
        raise ValueError(
            f"relocation_handling={relocation_handling!r} -- must be "
            f"'longest_segment', 'latest_segment' or 'full_record'."
        )

    # ── Station drift transect ───────────────────────────────────────────
    station_file = stations_dict[station]["file"]

    jumps = find_position_jumps(station_file, jump_threshold_m=jump_threshold_m)
    segment_start = segment_end = None
    if jumps:
        jump_desc = "; ".join(
            f"{j['time_before'].date()} -> {j['time_after'].date()} "
            f"({j['jump_m']:,.0f} m)".replace(",", " ")
            for j in jumps
        )
        print(
            f"[{station}] {len(jumps)} position jump(s) detected -- likely "
            f"AWS relocation(s), not ice-flow drift: {jump_desc}."
        )
        if relocation_handling == "longest_segment":
            segment_start, segment_end = longest_continuous_segment(
                station_file, jump_threshold_m=jump_threshold_m
            )
            span_days = (segment_end - segment_start).days
            print(
                f"[{station}] relocation_handling='longest_segment' -- "
                f"restricting to the longest uninterrupted deployment, "
                f"{segment_start.date()} to {segment_end.date()} "
                f"({span_days} days). Pass relocation_handling="
                f"'latest_segment' to always use the station's current "
                f"site instead (even if shorter), or 'full_record' to "
                f"ignore relocations entirely."
            )
        elif relocation_handling == "latest_segment":
            segment_start = jumps[-1]["time_after"]
            print(
                f"[{station}] relocation_handling='latest_segment' -- "
                f"restricting to fixes on/after {segment_start.date()} (the "
                f"station's current site). Pass relocation_handling="
                f"'longest_segment' to use the longest uninterrupted "
                f"deployment instead (even if it's not the most recent "
                f"one), or 'full_record' to ignore relocations entirely."
            )

    transect = build_transect_from_station(
        station_file, extend_km=extend_km, n_points=n_points,
        min_time=segment_start, max_time=segment_end,
    )
    x_line = transect["x_line"]
    y_line = transect["y_line"]
    distance_km = transect["distance_km"]
    bearing_deg = transect["bearing_deg"]
    compass = bearing_to_compass(bearing_deg)

    # ── Bed topography (BedMachine) ────────────────────────────────────
    bed_ds = xr.open_dataset(BED_TOPOGRAPHY_DATA_FILE, chunks={})

    if bed_ref_date is None:
        nominal_year = bed_ds.attrs.get("nominal_year", 2008)
        bed_ref_date = f"{int(nominal_year)}-01-01"
    bed_ref_ts = pd.to_datetime(bed_ref_date)

    bed_line = sample_along_line(bed_ds["bed"], x_line, y_line)
    ref_surface_line = sample_along_line(bed_ds["surface"], x_line, y_line)
    errbed_line = sample_along_line(bed_ds["errbed"], x_line, y_line)

    # ── Per-satellite dh sampled along the line, re-anchored to bed_ref_date ──
    # dh_by_satellite[name] = (time_index, dh_along_line[time, points])
    dh_by_satellite = {}
    for sat_name in list_satellites:
        info_sat = interpolation_altimetry_AWS.SATELLITE[sat_name]
        ds = interpolation_altimetry_AWS.satellite_opening(info_sat["file"])
        var = ds[info_sat["var"]]

        # Sample the full time series along the transect (points dim),
        # keeping time so we can later group by hydrological year.
        x_idx = xr.DataArray(x_line, dims="points")
        y_idx = xr.DataArray(y_line, dims="points")
        var_line = var.sel(x=x_idx, y=y_idx, method="nearest")  # dims: (time, points)

        # Value at the reference date, per point along the line.
        dh_ref = var_line.sel(time=bed_ref_ts, method="nearest")

        dh_by_satellite[sat_name] = {
            "time": pd.to_datetime(var_line["time"].values),
            "dh_reanchored": (var_line - dh_ref),  # xr.DataArray (time, points)
        }

    # ── Determine which hydrological years are common to all satellites ──
    years_per_sat = []
    for sat_name, d in dh_by_satellite.items():
        years = sorted(set(hydrological_year(t) for t in d["time"]))
        years_per_sat.append(set(years))
    common_years = sorted(set.intersection(*years_per_sat))

    if start_year is not None:
        common_years = [y for y in common_years if y >= start_year]
    if end_year is not None:
        common_years = [y for y in common_years if y <= end_year]

    if not common_years:
        raise ValueError(
            "No hydrological year is common to all requested satellites "
            "over the requested start_year/end_year range."
        )

    # ── Mean multi-satellite surface per hydrological year ──────────────
    surfaces_by_year = {}
    for year in common_years:
        per_sat_surfaces = []
        for sat_name, d in dh_by_satellite.items():
            time_index = d["time"]
            mask = np.array([hydrological_year(t) == year for t in time_index])
            if not mask.any():
                continue
            dh_year_mean = d["dh_reanchored"].isel(time=np.where(mask)[0]).mean(
                dim="time"
            ).values
            surface_this_sat = ref_surface_line + dh_year_mean
            per_sat_surfaces.append(surface_this_sat)

        if per_sat_surfaces:
            surfaces_by_year[year] = np.nanmean(np.vstack(per_sat_surfaces), axis=0)

    years_sorted = sorted(surfaces_by_year.keys())

    # ── Station position at the start/end of the plotted period ──────────
    # The two marked points must fall within the period actually covered by
    # the plotted hydrological years (i.e. consistent with the colorbar on
    # the right), not the station's absolute first/last GPS fix -- which
    # can extend well before or after the years shown. They must also stay
    # within [segment_start, segment_end] (if a relocation was detected),
    # for the same reason the transect itself is restricted above --
    # otherwise "start"/"end" could land on the station's old, physically
    # different site, or on a shorter/different segment than the one the
    # transect was actually built from.
    target_start = pd.Timestamp(f"{years_sorted[0] - 1}-09-01")
    target_end = pd.Timestamp(f"{years_sorted[-1]}-08-31")
    period_start = max(target_start, segment_start) if segment_start is not None else target_start
    period_end = min(target_end, segment_end) if segment_end is not None else target_end
    if period_start > period_end:
        raise ValueError(
            f"[{station}] the retained deployment segment "
            f"({segment_start.date()} to {segment_end.date()}) doesn't "
            f"overlap the plotted hydrological years "
            f"({target_start.date()} to {target_end.date()}) -- try a "
            f"different relocation_handling, or start_year/end_year."
        )
    t_mark_first, lon_mark_first, lat_mark_first = get_station_position_near(
        station_file, period_start, period_start=period_start, period_end=period_end
    )
    t_mark_last, lon_mark_last, lat_mark_last = get_station_position_near(
        station_file, period_end, period_start=period_start, period_end=period_end
    )

    x_mark_first, y_mark_first = _LATLON_TO_3413.transform(lon_mark_first, lat_mark_first)
    x_mark_last, y_mark_last = _LATLON_TO_3413.transform(lon_mark_last, lat_mark_last)

    # Straight-line (planimetric) horizontal displacement of the station
    # between the two marked points, in EPSG:3413 (i.e. true ground
    # distance, not the along-transect distance -- the two coincide here
    # since the transect is built along the first->last drift direction,
    # but computed independently from x/y so it stays correct even if that
    # ever changes).
    horizontal_displacement_m = float(
        np.hypot(x_mark_last - x_mark_first, y_mark_last - y_mark_first)
    )

    # Map-consistent orientation: west Greenland (station longitude west of
    # 45°W) displays the terminus (downstream) side on the left; east
    # Greenland (longitude east of 45°W) keeps the natural left-to-right
    # order -- matching how the transect would appear on a standard
    # north-up map (west=left, east=right).
    station_lon = (lon_mark_first + lon_mark_last) / 2
    invert_axis = station_lon < -45

    # Project onto the transect line (it may not fall exactly on one of the
    # n_points samples) to get its position in km along distance_km.
    dist_mark_first_km = project_km(
        x_mark_first, y_mark_first, transect["x_start"], transect["y_start"],
        transect["ux"], transect["uy"],
    )
    dist_mark_last_km = project_km(
        x_mark_last, y_mark_last, transect["x_start"], transect["y_start"],
        transect["ux"], transect["uy"],
    )

    # thickness/errbed are intentionally NOT sampled here anymore: they were
    # only used by the removed per-point dashed ice column (see below) --
    # bed and surface are still needed (bed for real_delta_thickness_at_station,
    # surface for the satellite re-anchoring and the reference-surface slope).
    bed_mark_first = sample_along_line(bed_ds["bed"], [x_mark_first], [y_mark_first])[0]
    surface_mark_first = sample_along_line(bed_ds["surface"], [x_mark_first], [y_mark_first])[0]

    bed_mark_last = sample_along_line(bed_ds["bed"], [x_mark_last], [y_mark_last])[0]
    surface_mark_last = sample_along_line(bed_ds["surface"], [x_mark_last], [y_mark_last])[0]

    # ── Satellite-derived altitude at the two marked points ──────────────
    # Mean +/- inter-product std of the requested satellites (list_satellites),
    # each sampled at the exact point (not along the transect line) and at
    # its acquisition nearest to that point's date -- see
    # satellite_surface_at_point()'s docstring for the re-anchoring formula.
    sat_vals_first, sat_used_first = satellite_surface_at_point(
        list_satellites, x_mark_first, y_mark_first, t_mark_first,
        surface_mark_first, bed_ref_ts,
    )
    alt_first, alt_unc_first = mean_std_across_products(sat_vals_first)

    sat_vals_last, sat_used_last = satellite_surface_at_point(
        list_satellites, x_mark_last, y_mark_last, t_mark_last,
        surface_mark_last, bed_ref_ts,
    )
    alt_last, alt_unc_last = mean_std_across_products(sat_vals_last)

    # Difference between the two points, with standard error propagation
    # for a difference of two (assumed independent) quantities:
    # sqrt(unc_first^2 + unc_last^2). If one side's uncertainty is missing
    # (n/a) it counts as 0 in the sum rather than voiding the whole result;
    # if the altitude itself is missing on either side, the difference is
    # n/a too.
    if np.isnan(alt_first) or np.isnan(alt_last):
        diff_alt, diff_unc = np.nan, np.nan
    else:
        diff_alt = alt_last - alt_first
        diff_terms = [u**2 for u in (alt_unc_first, alt_unc_last) if not np.isnan(u)]
        diff_unc = np.sqrt(sum(diff_terms)) if diff_terms else np.nan

    # ── Reference-surface slope between the two points ────────────────────
    # Purely geometric contribution: how much the static BedMachine
    # reference surface itself differs between the station's start and end
    # positions (744 m apart in the JAR example), independently of any real
    # temporal elevation change. Since diff_alt = ref_surface_diff +
    # (real dh change, re-anchored at each point), this line lets a reader
    # see how much of diff_alt is "just" advection along the local slope.
    ref_surface_diff = float(surface_mark_last - surface_mark_first)

    # ── Real, time-resolved ice-thickness change (ArcticDEM strips) ──────
    # Samples an individually-dated DEM strip at each position/date and
    # differences against the (assumed static) bed -- see
    # real_delta_thickness_at_station's docstring for the co-registration
    # caveat. (No BedMachine-only spatial-gradient version is computed --
    # see the module docstring's "Modelling choice" note for why.)
    # Modulated by compute_real_thickness: this is the only step in the
    # figure that needs network access (a PGC STAC query, then downloading
    # part of a strip) -- skip it (no query at all) by passing
    # compute_real_thickness=False, e.g. for an offline run or a batch call
    # over many stations where the round-trip time adds up.
    if compute_real_thickness:
        real_thickness_result = real_delta_thickness_at_station(
            x_mark_first, y_mark_first, t_mark_first,
            x_mark_last, y_mark_last, t_mark_last,
            bed_mark_first, bed_mark_last,
            search_window_days=real_thickness_search_window_days,
        )
    else:
        real_thickness_result = {
            "delta_thickness_m": np.nan,
            "thickness_first_m": np.nan,
            "thickness_last_m": np.nan,
            "strip_date_first": None,
            "strip_date_last": None,
            "strip_id_first": None,
            "strip_id_last": None,
            "skipped_reason": "compute_real_thickness=False",
        }

    # ── Ice surface velocity, mean per hydrological year ─────────────────
    velocity_by_year, velocity_units = {}, ""
    if include_velocity:
        velocity_by_year, velocity_units = build_velocity_by_year(
            x_line, y_line, years_sorted,
            velocity_dir=velocity_dir, pattern=velocity_pattern,
        )

    # ── Plot ──────────────────────────────────────────────────────────
    # Two stacked panels sharing the x-axis:
    #   - top:    everything except velocity -- bed, bed uncertainty,
    #             reference surface, and the yearly mean multi-satellite
    #             surface curves, all on the same elevation axis.
    #   - bottom: ice surface velocity magnitude, mean per hydrological
    #             year (Sentinel-1 offset-tracking).
    cmap = plt.get_cmap(cmap_name)
    norm = plt.Normalize(vmin=min(years_sorted), vmax=max(years_sorted))

    fig, (ax_full, ax_vel) = plt.subplots(
        2, 1, figsize=(13, 11), sharex=True,
        gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08},
    )

    # Set the x-limits (and their left/right orientation) up front, before
    # any of the collision-avoidance/annotation placement below -- that
    # logic needs each marked point's actual on-screen (pixel) position,
    # which depends on this orientation (see invert_axis above), not just
    # on the raw distance_km values. Data limits don't need anything to be
    # plotted first, so this is safe to do immediately.
    if invert_axis:
        ax_full.set_xlim(distance_km.max(), distance_km.min())
    else:
        ax_full.set_xlim(distance_km.min(), distance_km.max())
    ax_vel.set_xlim(ax_full.get_xlim())

    # -- Top panel: bed, bed uncertainty, reference surface, yearly curves --
    ax_full.fill_between(
        distance_km, np.nanmin(bed_line) - 100, bed_line,
        color="#D4B896", alpha=0.85, zorder=1, label="Bed topography",
    )
    ax_full.fill_between(
        distance_km, bed_line - errbed_line, bed_line + errbed_line,
        color="#5C3A1E", alpha=0.6, zorder=1.5, label="Bed uncertainty (± errbed)",
    )
    ax_full.plot(
        distance_km, ref_surface_line,
        color="0.4", linewidth=1.2, linestyle="--", zorder=2,
        label=f"Reference surface (BedMachine, {bed_ref_ts.year})",
    )
    for year in years_sorted:
        ax_full.plot(
            distance_km, surfaces_by_year[year],
            color=cmap(norm(year)), linewidth=1.6, zorder=3,
        )

    # -- Station position at start/end of the plotted period: name + date label --
    # (No dashed ice column here anymore -- the per-point vertical dotted
    # line, its triangle markers, bed-uncertainty error bar and thickness
    # text label were removed on request, to avoid any reading of "a change
    # between the two points" from a purely per-point, static-BedMachine
    # quantity. Only the station-name/date callouts remain, to mark where
    # along the transect each point falls.)
    #
    # Collision-avoidance: when the start/end points are close together
    # (short drift over the plotted period, or a small extend_km), their
    # labels would otherwise land on top of each other -- the two
    # "station name + date" callouts sit at the same axes-fraction height,
    # and the two rotated thickness labels sit right next to each other.
    # In that case each label is nudged apart from the other by a small
    # on-screen offset (in points) and linked back to its true position
    # with a thin leader line. Far-apart points (the common case) get zero
    # offset, i.e. unchanged behaviour.
    #
    # Both "is it crowded" and "which direction to push each label" are
    # decided from the two points' actual on-screen pixel positions
    # (ax_full.transData, using the x-limits set above), not from
    # dist_mark_first_km/last_km or invert_axis directly. Two earlier
    # versions of this logic got this wrong:
    #   - using a fixed push direction (first always left, last always
    #     right) ignored invert_axis, so on a reversed axis it nudged each
    #     label *towards* the other marker instead of away from it -- this
    #     is what made "JAR 2010" look more downstream than "JAR 2023".
    #   - even after fixing that by reading invert_axis, it still assumed
    #     dist_mark_first_km < dist_mark_last_km always (i.e. that "first"
    #     is always upstream of "last" along the transect). For a station
    #     with very little net drift (tens of metres, e.g. QAS_L or
    #     UPE_L), that projected order isn't reliable -- GPS noise on the
    #     two independently-picked "nearest fix to target date" points can
    #     put them in either order -- so a push direction based on
    #     first/last identity alone could still point the wrong way.
    #   - deciding "crowded" from a fixed 15%-of-span threshold in km also
    #     doesn't account for the actual figure size/DPI/font: a
    #     separation just above that threshold (e.g. NUK_U, ~16% of span)
    #     was still only a few tens of pixels on screen -- not enough room
    #     for the two-line bold "STATION\nYYYY-MM-DD" label -- and got
    #     zero offset, producing overlapping/garbled text.
    # Measuring the actual pixel gap side-steps all three: it's correct
    # regardless of axis direction, regardless of whether "first" or
    # "last" ends up further along the transect, and it reflects what will
    # really overlap on screen instead of a proxy for it.
    px_first = ax_full.transData.transform((dist_mark_first_km, 0))[0]
    px_last = ax_full.transData.transform((dist_mark_last_km, 0))[0]
    sep_pts = abs(px_last - px_first) * 72.0 / fig.dpi  # pixels -> points (DPI-independent)

    min_gap_pts = 90  # ~ the on-screen width needed for the two-line bold label
    is_crowded = sep_pts < min_gap_pts

    for dist_km, t_mark in [
        (dist_mark_first_km, t_mark_first),
        (dist_mark_last_km, t_mark_last),
    ]:
        # When the two points are far enough apart (not is_crowded), each
        # gets its own "STATION\nDATE" callout directly above it. When
        # they're crowded, a single merged callout is drawn once below
        # instead (see after this loop) -- two separate callouts pushed
        # apart by a fixed point offset still have to send their leader
        # lines back down to two nearly-adjacent anchors, which for
        # near-stationary stations (tens to a couple hundred metres of net
        # drift, e.g. QAS_L, UPE_L, KPC_U) reads as two lines awkwardly
        # converging on almost the same spot rather than as two clearly
        # separate pointers -- a single "STATION\nDATE1 -> DATE2" label is
        # both clearer and avoids that visual altogether.
        if not is_crowded:
            ax_full.annotate(
                f"{station}\n{pd.Timestamp(t_mark).date()}",
                xy=(dist_km, 0.96), xycoords=("data", "axes fraction"),
                ha="center", va="top", fontsize=8.5, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5),
                zorder=7,
            )

    if is_crowded:
        mid_dist_km = (dist_mark_first_km + dist_mark_last_km) / 2
        # Order the two dates left-to-right to match how they actually
        # render on screen (px_first/px_last, computed above -- already
        # accounts for invert_axis), not chronological order. On an
        # inverted axis (west Greenland) the station can genuinely drift
        # from its start position on the right to its end position on the
        # left -- terminus is placed on the left there, matching a
        # north-up map, so "later" can render left of "earlier". Writing
        # the merged label as "first -> last" regardless would then read
        # backwards compared to the plot (the earlier date appearing on
        # the label's left even though its point is on the right).
        date_first_str = str(pd.Timestamp(t_mark_first).date())
        date_last_str = str(pd.Timestamp(t_mark_last).date())
        dates_left_to_right = (
            f"{date_first_str} → {date_last_str}" if px_first <= px_last
            else f"{date_last_str} → {date_first_str}"
        )
        ax_full.annotate(
            f"{station}\n{dates_left_to_right}",
            xy=(mid_dist_km, 0.96), xycoords=("data", "axes fraction"),
            ha="center", va="top", fontsize=8.5, fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5),
            zorder=7,
        )
    ax_full.set_ylabel("Elevation (m a.s.l.)")
    ax_full.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)

    # -- Legend: usual entries + the point-to-point comparison lines --
    if np.isnan(diff_alt):
        legend_diff = "Δ altitude (end − start): n/a"
    elif np.isnan(diff_unc):
        legend_diff = f"Δ altitude (end − start): {diff_alt:+.1f} m (unc. n/a)"
    else:
        legend_diff = f"Δ altitude (end − start): {diff_alt:+.1f} ± {diff_unc:.1f} m"

    legend_ref_surface_diff = (
        f"Reference-surface slope (BedMachine, end − start): "
        f"{ref_surface_diff:+.1f} m"
    )

    # Note: plain "," grouping (e.g. "1,611") reads as an English thousands
    # separator but looks like a French decimal number ("1.611") to a
    # French-speaking reader -- use a space instead, which is unambiguous
    # (and matches the SI/French convention for digit grouping).
    legend_displacement = (
        f"Horizontal displacement (start → end): "
        f"{horizontal_displacement_m:,.0f}".replace(",", " ") + " m"
    )
    legend_real_thickness_diff = format_real_delta_thickness_legend(
        real_thickness_result
    )

    extra_handles = [
        Line2D([], [], color="none", label=legend_diff),
        Line2D([], [], color="none", label=legend_ref_surface_diff),
        Line2D([], [], color="none", label=legend_real_thickness_diff),
        Line2D([], [], color="none", label=legend_displacement),
    ]
    handles, labels = ax_full.get_legend_handles_labels()
    handles += extra_handles
    labels += [h.get_label() for h in extra_handles]
    ax_full.legend(handles, labels, loc="lower left", fontsize=9)

    ax_full.set_title(
        f"Vertical cross-section along {station}'s drift direction "
        f"(+{extend_km} km on each end)\n"
        f"mean of {', '.join(list_satellites)}",
        fontsize=12, fontweight="bold", pad=48,
    )

    # -- Flow-direction annotation: accumulation zone <-> terminus --
    # The station drifts with the ice: its position at the start of the
    # plotted period (low-distance end) is upstream, closer to the
    # accumulation zone; its position at the end of the plotted period
    # (high-distance end) is downstream, closer to the glacier
    # terminus/ablation zone. Left/right placement of the two labels
    # follows invert_axis so the figure stays consistent with a north-up
    # map (west on the left, east on the right).
    if invert_axis:
        left_text, right_text = "Glacier terminus (downstream)", "Accumulation zone (upstream)"
    else:
        left_text, right_text = "Accumulation zone (upstream)", "Glacier terminus (downstream)"

    ax_full.annotate(
        "", xy=(0.02, 1.055), xytext=(0.14, 1.055), xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", color="black", lw=1.4),
    )
    ax_full.text(
        0.02, 1.075, left_text,
        transform=ax_full.transAxes, ha="left", va="bottom",
        fontsize=9, fontweight="bold",
    )
    ax_full.annotate(
        "", xy=(0.98, 1.055), xytext=(0.86, 1.055), xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", color="black", lw=1.4),
    )
    ax_full.text(
        0.98, 1.075, right_text,
        transform=ax_full.transAxes, ha="right", va="bottom",
        fontsize=9, fontweight="bold",
    )

    # -- Bottom panel: ice surface velocity magnitude --
    if velocity_by_year:
        for year, vel_line in velocity_by_year.items():
            ax_vel.plot(
                distance_km, vel_line,
                color=cmap(norm(year)), linewidth=1.6, zorder=3,
            )
        missing_years = [y for y in years_sorted if y not in velocity_by_year]
        if missing_years:
            print(
                "No velocity acquisition found for hydrological year(s): "
                f"{missing_years} -- these years are absent from the velocity panel."
            )
    else:
        ax_vel.text(
            0.5, 0.5, "No velocity data available",
            ha="center", va="center", transform=ax_vel.transAxes,
            fontsize=10, style="italic", color="0.4",
        )
    # (x-limits already set earlier, before the collision-avoidance block,
    # since that logic needs each marked point's on-screen pixel position.)
    ax_vel.set_xlabel(
        f"Distance along {station}'s drift-direction transect (km) "
        f"— west on the left, east on the right, as on a map"
    )
    ylabel_vel = "Ice surface velocity magnitude"
    if velocity_units:
        ylabel_vel += f" ({velocity_units})"
    ax_vel.set_ylabel(ylabel_vel)
    ax_vel.grid(True, linestyle="--", linewidth=0.4, alpha=0.5)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[ax_full, ax_vel], pad=0.015)
    cbar.set_label("Hydrological year")

    if save:
        if output_path is None:
            output_path = (
                FIGURES_DIR
                / f"cross_section_{station}_drift_"
                  f"{years_sorted[0]}-{years_sorted[-1]}.png"
            )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Figure saved: {output_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, (ax_full, ax_vel), surfaces_by_year, velocity_by_year


# ── Batch plotting over every station of a STATIONS_* dict ────────────────
def plot_all_stations(
    stations_dict=AWS_data.STATIONS_ablation,
    stations=None,
    stop_on_error=False,
    **kwargs,
):
    """Call plot_cross_section() once per station, so you don't have to
    change the station name by hand every time.

    stations_dict : dict of stations to loop over -- same shape as
        AWS_data.STATIONS_ablation (its default), i.e. {station_name:
        {"file": ..., ...}}. Pass e.g. AWS_data.STATIONS_accumulation to
        loop over a different set.
    stations : optional subset/order of station keys to plot. Defaults to
        every key in stations_dict, in the order they appear there.
    stop_on_error : if False (default), a station whose plot raises an
        exception (e.g. a degenerate first==last GPS fix, no hydrological
        year common to every requested satellite, no velocity file for the
        requested period) is skipped with a printed warning instead of
        aborting the whole loop -- handy since data quality/coverage varies
        a lot from station to station. Set True to stop at the first error
        instead (e.g. while debugging).
    **kwargs : forwarded as-is to plot_cross_section() for every station
        (e.g. extend_km=10, include_velocity=False, start_year=2015,
        cmap_name="plasma"). Do not pass `station` or `stations_dict` here
        -- they are set by this function. `show` defaults to False so the
        loop doesn't block waiting for each figure window to be closed
        one by one; pass show=True to override and inspect each figure as
        it is produced.

    Returns {station_name: result_or_None}, where result is the
    (fig, (ax_full, ax_vel), surfaces_by_year, velocity_by_year) tuple
    returned by plot_cross_section for that station, or None if the
    station was skipped after an error.
    """
    if stations is None:
        stations = list(stations_dict.keys())

    kwargs.setdefault("show", False)

    results = {}
    for station in stations:
        print(f"── {station} " + "─" * max(0, 40 - len(station)))
        try:
            results[station] = plot_cross_section(
                station=station, stations_dict=stations_dict, **kwargs
            )
        except Exception as exc:
            results[station] = None
            print(f"  skipped ({type(exc).__name__}: {exc})")
            if stop_on_error:
                raise

    n_ok = sum(1 for v in results.values() if v is not None)
    print(f"\n{n_ok}/{len(stations)} station(s) plotted successfully.")
    return results


if __name__ == "__main__":
    plot_all_stations(
        stations_dict=AWS_data.STATIONS_ablation,
        extend_km=5,
    )