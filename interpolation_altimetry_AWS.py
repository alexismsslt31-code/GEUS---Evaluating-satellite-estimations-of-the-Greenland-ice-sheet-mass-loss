import pandas as pd
import xarray as xr
import numpy as np
import scipy
from scipy.spatial import cKDTree
from netCDF4 import Dataset
import os
import pyproj
from pathlib import Path

import AWS_data
from paths import SATELLITE_DATA_DIR

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Satellite product file paths ─────────────────────────────────────────
"""
SATELLITE -- per-satellite-product dict (file path, primary variable name, plot color), keyed by product name (Khan et al. 2025, Nilsson and Gardner 2026, Andersen et al. 2025, Copernicus Climate Data Store, Zhang et al. 2022).
"""

# ── Opening data ─────────────────────────────────────────────────────────
"""
def satellite_opening(path) -- opens a satellite product NetCDF as an xarray Dataset.

def aws_opening(path) -- reads an AWS CSV as a pandas DataFrame.
"""

# ── Matching a satellite product to an AWS station's drift ────────────────
"""
def satellite_points_selection(aws_path, satellite_path, variable) -- for each satellite acquisition time, finds the AWS station's position closest in time (nearest-neighbour on a KD-tree built from the satellite grid's lat/lon) and samples the requested variable at that pixel, following the station's drift over time rather than a single fixed point.

def satellite_on_aws(aws_path, satellite_name, variable) -- calls satellite_points_selection for a named product and merges the sampled satellite series with the AWS record on time (outer join).
"""


# ── Chemins vers les fichiers de données satellite ──────────────────────────────────────────────────────────────
SATELLITE = {
    "Khan et al., 2025" : {"file": SATELLITE_DATA_DIR/ "Khan et al., 2025/Greenland_netcdf_1kmgrid_DB/Greenland_dh_icevol_1kmgrid_DB.nc", "var": "dh_vol", "color": "#D20DF0"},
    "Nilsson and Gardner, 2026" : {"file": SATELLITE_DATA_DIR/ "Nilsson and Gardner, 2026/Lat_Lon_Greenland_G1920V01_IceSheetGlacierIceHeight.nc", "var": "dh", "color": "#17EE17"},
    "Andersen et al., 2025" : {"file": SATELLITE_DATA_DIR/ "Andersen et al., 2025/Lat_Lon_CCI_GrIS_RA_dSEC_dh_5km_012011_032025.nc", "var": "ZZ", "color": "#EDE20D"},
    "Copernicus_Climate_Data_Store" : {"file" : SATELLITE_DATA_DIR/ "Copernicus Climate Data Store/Lat_Lon_C3S_GrIS_RA_SEC_25km_Vers6_199108-202601_2026-04-18.nc", "var": "dh", "color": "#EB2F25"},
    "Zhang et al., 2022" : {"file" : SATELLITE_DATA_DIR/ "Zhang et al., 2022/Time_x_y_Surface_Elevation_Anomaly_Greenland_Monthly_5km_Grid.nc", "var": "elev_interp", "color": "#03FEF1"},

    
    # "DEM" :{"file" : SATELLITE_DATA/ "Greenland_DEM_1kmgrid_DB.nc"},
    # "water_equivalent" :{"file" : SATELLITE_DATA/ "Greenland_dhdt_mass_1kmgrid_DB.nc"},
    # "firn" : {"file" : SATELLITE_DATA/"Greenland_dhdt_firn_1kmgrid_DB.nc"}
}
# ── 1. Ouverture des données satellite et aws ────────────────────────────────────────────────────────
def satellite_opening(path) :
    satellite_data = xr.open_dataset(path)
    return satellite_data

#print(satellite_opening(SATELLITE_DATA/ "Lat_Lon_C3S_GrIS_RA_SEC_25km_Vers6_199108-202601_2026-04-18.nc" ))
#print(satellite_opening(SATELLITE_DATA/"Time_x_y_Surface_Elevation_Anomaly_Greenland_Monthly_5km_Grid.nc"))

def aws_opening(path) :
    aws_data = pd.read_csv(path)
    return aws_data

def satellite_points_selection(aws_path, satellite_path,variable):
    aws = aws_opening(aws_path)
    satellite = satellite_opening(satellite_path)
    aws["time"] = pd.to_datetime(aws["time"])

    # --- Coordonnées initiales de la station ---
    first_valid = aws[["lat", "lon"]].dropna().iloc[0]
    lat_0 = first_valid["lat"]
    lon_0 = first_valid["lon"]

    # --- Grille satellite ---
    lat_grid = satellite["lat"].values
    lon_grid = satellite["lon"].values
    lat_flat = lat_grid.ravel()
    lon_flat = lon_grid.ravel()

    def latlon_to_xyz(lat, lon):
        lat_r, lon_r = np.deg2rad(lat), np.deg2rad(lon)
        return np.column_stack([
            np.cos(lat_r) * np.cos(lon_r),
            np.cos(lat_r) * np.sin(lon_r),
            np.sin(lat_r)
        ])

    # --- KDTree ---
    xyz_all = latlon_to_xyz(lat_flat, lon_flat)
    valid_mask = np.isfinite(xyz_all).all(axis=1)
    valid_indices = np.where(valid_mask)[0]
    xyz_valid = xyz_all[valid_mask]
    tree = cKDTree(xyz_valid)

    # --- Pixel initial (lat_0, lon_0) ---
    _, idx_in_valid_0 = tree.query(latlon_to_xyz(np.array([lat_0]), np.array([lon_0])))
    flat_idx_0 = valid_indices[idx_in_valid_0[0]]
    y_idx_0, x_idx_0 = np.unravel_index(flat_idx_0, lat_grid.shape)

    # --- Temps ---
    sat_times = satellite["time"].values
    aws_time_min = aws["time"].min().to_datetime64()

    # --- Boucle principale ---
    last_valid_lat = lat_0
    last_valid_lon = lon_0
    results = []

    for i, sat_time in enumerate(sat_times):

        if sat_time < aws_time_min:
            y_idx, x_idx = y_idx_0, x_idx_0

        else:
            closest_aws_idx = np.argmin(np.abs(aws["time"].values - sat_time))
            lat_aws = aws["lat"].iloc[closest_aws_idx]
            lon_aws = aws["lon"].iloc[closest_aws_idx]

            if not np.isfinite(lat_aws) or not np.isfinite(lon_aws):
                lat_aws = last_valid_lat
                lon_aws = last_valid_lon
            else:
                last_valid_lat = lat_aws
                last_valid_lon = lon_aws

            _, idx_in_valid = tree.query(latlon_to_xyz(np.array([lat_aws]), np.array([lon_aws])))
            flat_idx = valid_indices[idx_in_valid[0]]
            y_idx, x_idx = np.unravel_index(flat_idx, lat_grid.shape)

        val_raw = satellite[variable].isel(time=i, y=y_idx, x=x_idx).values
        try:
            val = float(val_raw)
        except (TypeError, ValueError):
            val = np.nan

        results.append({
            "time_sat": sat_time,
            "y_idx":    y_idx,
            "x_idx":    x_idx,
            "lat_used": lat_grid[y_idx, x_idx],
            "lon_used": lon_grid[y_idx, x_idx],
            variable: val,
        })

    return pd.DataFrame(results)


def satellite_on_aws(aws_path, satellite_name, variable):
    """aws_path : direction vers le fichier aws
    satellite_name : str
    variable : str"""
    aws = aws_opening(aws_path)
    satellite_data = satellite_points_selection(aws_path, SATELLITE[satellite_name]["file"],variable)
    satellite_data["time_sat"] = pd.to_datetime(satellite_data["time_sat"])
    aws["time"] = pd.to_datetime(aws["time"])
    aws_satellite_data = pd.merge(
        satellite_data,
        aws,
        left_on="time_sat",
        right_on="time",
        how="outer"
    )
    return (aws_satellite_data)

# print(satellite_on_aws(AWS_data.STATION_KAN_U['KAN_U']['hourly_data'], 'Copernicus_Climate_Data_Store', 'dh'))