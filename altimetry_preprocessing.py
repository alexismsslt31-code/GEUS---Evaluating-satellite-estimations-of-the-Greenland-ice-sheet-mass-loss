import pandas as pd
import xarray as xr
import numpy as np
import pyproj
from netCDF4 import Dataset
import dask
from pathlib import Path
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt

from paths import DATA

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Paths (legacy pre-processing layout) ───────────────────────────────────
"""
AWS_DATA, SATELLITE_DATA, SATELLITE, AWS_STATIONS -- this script pre-processes the raw satellite/AWS downloads into the final NetCDF/CSV layout paths.py and the rest of the codebase read from. DATA is imported from paths.py, but the subfolder names below ("weather_stations", "satellites") are the RAW, pre-processing-stage layout, distinct from paths.py's final "weather stations"/"satellite products" layout -- kept separate deliberately, not a duplicate to fix.
"""

# ── Adding lat/lon coordinates to each product's dataset ──────────────────
"""
Module-level code below (outside any function) reprojects each satellite product's native grid to lat/lon coordinates (or renames/reformats existing ones) and writes the result back out as a NetCDF ready for the rest of the codebase. Most of these blocks are commented out (already run once); the Zhang et al. block is the one currently active.
"""


# ── Chemins d'accès ──────────────────────────────────────────────────────────────
AWS_DATA = DATA/ "weather_stations/"
SATELLITE_DATA = DATA/ "satellites"

# ── Chemins vers les fichiers de données satellite ──────────────────────────────────────────────────────────────
SATELLITE = {
    "khan" : {"file": SATELLITE_DATA/ "doi_10_5061_dryad_s4mw6m9dh__v20250422/Greenland_netcdf_1kmgrid_DB/Greenland_dhdt_icevol_1kmgrid_DB.nc"},
    "nilsson" : {"file": SATELLITE_DATA/ "ERS_EnviSat_ICESat_CryoSat-2_1992_2023/Greenland_G1920V01_IceSheetGlacierIceHeight.nc"},
    "andersen" : {"file": SATELLITE_DATA/ "CryoSat-2_2011_2025/CCI_GrIS_RA_dSEC_dh_5km_012011_032025.nc"},
    "copernicus" : {"file": SATELLITE_DATA/ "Copernicus_Climate_Data_Store/C3S_GrIS_RA_SEC_25km_Vers6_199108-202601_2026-04-18.nc"},
    "zhang" : {"file": SATELLITE_DATA/ "Zhang/Surface_Elevation_Anomaly_Greenland_Monthly_5km_Grid.nc"}
}

# ── Chemins vers les fichiers de données AWS ──────────────────────────────────────────────────────────────
AWS_STATIONS = {
    #Accumulation
    "CEN": {"file": AWS_DATA/ "CEN_month.csv"},
    "HUM": {"file": AWS_DATA/ "HUM_month.csv"},
    "NEM": {"file": AWS_DATA/ "NEM_month.csv"},
    "TUN": {"file": AWS_DATA/ "TUN_month.csv"},
    "EGP": {"file": AWS_DATA/ "EGP_month.csv"},
    "NAE": {"file": AWS_DATA/ "NAE_month.csv"},
    "ZAC_A": {"file": AWS_DATA/ "ZAC_A_month.csv"},
    "NAU": {"file": AWS_DATA/ "NAU_month.csv"},
    "CP1": {"file": AWS_DATA/ "CP1_month.csv"},
    "DY2": {"file": AWS_DATA/ "DY2_month.csv"},
    "NSE": {"file": AWS_DATA/ "NSE_month.csv"},
    "SDL": {"file": AWS_DATA/ "SDL_month.csv"},
    "SDM": {"file": AWS_DATA/ "SDM_month.csv"},

    #Ablation
    "THU_L": {"file": AWS_DATA/ "THU_L_month.csv"},
    "UPE_L": {"file": AWS_DATA/ "UPE_L_month.csv"},
    "UPE_U": {"file": AWS_DATA/ "UPE_U_month.csv"},
    "SWC": {"file": AWS_DATA/ "SWC_month.csv"},
    "KAN_L": {"file": AWS_DATA/ "KAN_L_month.csv"},
    "KAN_M": {"file": AWS_DATA/ "KAN_M_month.csv"},
    "KAN_U": {"file": AWS_DATA/ "KAN_U_month.csv"},
    "NUK_U": {"file": AWS_DATA/ "NUK_U_month.csv"},
    "QAS_L": {"file": AWS_DATA/ "QAS_L_month.csv"},
    "QAS_M": {"file": AWS_DATA/ "QAS_M_month.csv"},
    "QAS_U": {"file": AWS_DATA/ "QAS_U_month.csv"},
    "TAS_L": {"file": AWS_DATA/ "TAS_L_month.csv"},
    "TAS_U": {"file": AWS_DATA/ "TAS_U_month.csv"},
    "TAS_A": {"file": AWS_DATA/ "TAS_A_month.csv"},
}

# ── Ajout de la latitude et longitude aux coordonnées des datasets si ce n'est pas le cas ────────────────────────────────────────────────────────

#________________________Nilsson______________________________
nilsson = xr.open_dataset(SATELLITE["nilsson"]["file"])

projection = pyproj.Proj("EPSG:3413")

x = nilsson.x.values  # shape (783,)
y = nilsson.y.values  # shape (1421,)

# Créer la grille 2D
xx, yy = np.meshgrid(x, y)  # shape (1421, 783)

# Transformer en lat/lon en une seule opération
lon_grid, lat_grid = projection(xx, yy, inverse=True)  # shape (1421, 783)

# --- Ajouter au dataset sans charger les données ---
nilsson = nilsson.assign_coords(
    lat=(["y", "x"], lat_grid.astype(np.float32)),
    lon=(["y", "x"], lon_grid.astype(np.float32))
)

# Enregistrement du dataset avec lat et lon en coordonnées
nilsson.to_netcdf(
    SATELLITE_DATA / "ERS_EnviSat_ICESat_CryoSat-2_1992_2023/Lat_Lon_Greenland_G1920V01_IceSheetGlacierIceHeight.nc",
    format="NETCDF4"
)

# #________________________Andersen______________________________
andersen = xr.open_dataset(SATELLITE["andersen"]["file"])
andersen = andersen.set_coords(
    ["Lat", "Lon"]
)
andersen = andersen.rename(
    {"Lat": "lat",
    "Lon": "lon"}
)

andersen = andersen.assign_coords(
    lat=(["y", "x"], andersen["lat"].values),
    lon=(["y", "x"], andersen["lon"].values)
)

# Temps sous format <U19 au lieu de datetime : ça va être galère à gérer donc on le repasse en datetime.
andersen = andersen.assign_coords(
    time=pd.to_datetime(andersen["time"].values)
)

# Enregistrement du dataset modifié et plus pratique
andersen.to_netcdf(
  SATELLITE_DATA/ "CryoSat-2_2011_2025/Lat_Lon_CCI_GrIS_RA_dSEC_dh_5km_012011_032025.nc"
)

#________________________Copernicus______________________________
copernicus = xr.open_dataset(SATELLITE["copernicus"]["file"])
copernicus = copernicus.set_coords(
    ["time","lat", "lon"]
)

copernicus = copernicus.rename(
    {"t": "time"}
)
# print(list(copernicus.data_vars))
print(copernicus)

copernicus.to_netcdf(
  SATELLITE_DATA/ "Copernicus_Climate_Data_Store/Lat_Lon_C3S_GrIS_RA_SEC_25km_Vers6_199108-202601_2026-04-18.nc"
)

#________________________Zhang______________________________
zhang = xr.open_dataset(SATELLITE["zhang"]["file"])
projection = pyproj.Proj("EPSG:3413")

lat = zhang.lat.values  # shape (614, 325)
lon = zhang.lon.values  # shape (614, 325)

# lat/lon → x/y projeté
xx, yy = projection(lon, lat, inverse=False)  # shape (614, 325)

# Extraire les vecteurs 1D (si la grille est régulière)
x_1d = xx[0, :].astype(np.float32)   # shape (x,)
y_1d = yy[:, 0].astype(np.float32)   # shape (y,)

# Renommer les dimensions d'abord
zhang = zhang.rename({'lat': 'y', 'lon': 'x'})

# Assigner les nouvelles coordonnées 2D
zhang = zhang.assign_coords(
    y=("y", y_1d),
    x=("x", x_1d),
    lat=(["y", "x"], lat.astype(np.float32)),
    lon=(["y", "x"], lon.astype(np.float32)),
)

zhang = zhang.assign_coords(time=pd.to_datetime(zhang.time.values.astype(int), format='%Y') + pd.to_timedelta((zhang.time.values % 1) * 365.25, unit='D'))

# plt.pcolormesh(zhang.lon, zhang.lat, zhang.elev_interp.isel(time=0))
# plt.colorbar()
# plt.show()

print(zhang)
zhang.to_netcdf(
    SATELLITE_DATA / "Zhang/Time_x_y_Surface_Elevation_Anomaly_Greenland_Monthly_5km_Grid.nc"
)

#________________________Khan______________________________

khan = xr.open_dataset(SATELLITE["khan"]["file"])

ds_check = xr.open_dataset(SATELLITE["khan"]["file"])
print(ds_check["dhdt_vol"].encoding.get("chunksizes"))
ds_check.close()

khan = xr.open_dataset(
    SATELLITE["khan"]["file"],
    chunks={"time": -1, "y": 504, "x": 293}  # time entier (nécessaire pour cumsum), y/x découpés
)

khan["dh_vol"] = khan["dhdt_vol"].cumsum(dim="time", skipna=False)

khan["dh_vol"].attrs = {
    "units": khan["dhdt_vol"].attrs.get("units", "m"),
    "long_name": "cumulative elevation change in ice equivalent",
    "standard_name": "cumulative_elevation_change_in_ice_equivalent", 
    "description": "Cumulative sum of dhdt_vol over time, computed as the running total of monthly elevation change rates from the start of the time series.",
}

encoding = {
    "dh_vol": {"zlib": True, "complevel": 4},
    "dhdt_vol": {"zlib": True, "complevel": 4},
}


with ProgressBar():
    khan.to_netcdf(
        SATELLITE_DATA / "doi_10_5061_dryad_s4mw6m9dh__v20250422/Greenland_netcdf_1kmgrid_DB/Greenland_dh_icevol_1kmgrid_DB.nc",
        encoding=encoding,
    )

print(khan)