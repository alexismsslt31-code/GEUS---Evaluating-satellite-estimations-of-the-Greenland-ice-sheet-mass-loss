from pathlib import Path

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Data directories ─────────────────────────────────────────────────────
"""
DATA, BED_TOPOGRAPHY_DATA_DIR, ICE_MASK_DATA_DIR, PRODEM_DATA_DIR, SATELLITE_DATA_DIR, VELOCITY_DATA_DIR, AWS_DATA_DIR -- single source of truth for every data subfolder used across the codebase, so no other module hardcodes its own copy of the base data path.
"""

# ── Project root and figures output ──────────────────────────────────────
"""
PROJECT_ROOT, FIGURES_DIR -- project root folder and the shared figures output folder every plotting module writes its PNGs under (in its own subfolder, e.g. FIGURES_DIR / "cross_sections").
"""

# ── Data files ────────────────────────────────────────────────────────────
"""
BED_TOPOGRAPHY_DATA_FILE, PRODEM_DATA_FILE, SATELLITE_DATA_FILE -- specific data files (rather than directories) referenced by more than one module: the BedMachine Greenland v6 NetCDF, the GC-Net precise-locations CSV, and the per-satellite-product NetCDF paths keyed by product name.
"""


# Paths leading to the data folder. Data are organised as follow :

# data/
# ├── bed topography/            
# ├── ice mask/
# ├── PRODEM/
# ├── satellite products/
#             ├── Andersen et al., 2025/
#             ├── Copernicus Climate Data Store/
#             ├── Khan et al., 2025/
#             ├── Nilsson and Gardner, 2026/
#             ├── Zhang et al., 2022/
# ├── velocity/
# ├── weather stations/
#             ├── day/
#             ├── hour/
#             ├── month/

#Directories
DATA = Path("C:/Users/alexi/Documents/ENM - Toulouse/Stage Copenhague/GEUS/data/")
BED_TOPOGRAPHY_DATA_DIR = DATA / "bed topography"
ICE_MASK_DATA_DIR = DATA / "ice mask"
PRODEM_DATA_DIR = DATA / "PRODEM"
SATELLITE_DATA_DIR = DATA / "satellite products"
VELOCITY_DATA_DIR = DATA / "velocity"
AWS_DATA_DIR = DATA / "weather stations/"

# Project root and shared figures output folder. Every plotting module in
# this codebase writes its PNGs under FIGURES_DIR (in its own subfolder,
# e.g. FIGURES_DIR / "cross_sections", FIGURES_DIR / "detrending"), instead
# of each hardcoding its own copy of this path.
PROJECT_ROOT = Path(
    "C:/Users/alexi/Documents/ENM - Toulouse/Stage Copenhague/GEUS/"
    "evaluating-satellite-estimations-of-the-greenland-ice-sheet-mass-loss"
)
FIGURES_DIR = PROJECT_ROOT / "figures"

#Files (when necessary)
BED_TOPOGRAPHY_DATA_FILE = BED_TOPOGRAPHY_DATA_DIR / "IDBMG4_6-20260819_071841/BedMachineGreenland-v6.nc"
# ICE_MASK_DATA_FILE = ICE_MASK_DATA_DIR / 
PRODEM_DATA_FILE = PRODEM_DATA_DIR / "GEUS_GC-Net_precise_locations.csv"
SATELLITE_DATA_FILE = {
    "Khan et al., 2025" : {"file": SATELLITE_DATA_DIR/ "Khan et al., 2025/Greenland_netcdf_1kmgrid_DB/Greenland_dh_icevol_1kmgrid_DB.nc"},
    "Nilsson and Gardner, 2026" : {"file": SATELLITE_DATA_DIR/ "Nilsson and Gardner, 2026/Lat_Lon_Greenland_G1920V01_IceSheetGlacierIceHeight.nc"},
    "Andersen et al., 2025" : {"file": SATELLITE_DATA_DIR/ "Andersen et al., 2025/Lat_Lon_CCI_GrIS_RA_dSEC_dh_5km_012011_032025.nc"},
    "Copernicus_Climate_Data_Store" : {"file" : SATELLITE_DATA_DIR/ "Copernicus Climate Data Store/Lat_Lon_C3S_GrIS_RA_SEC_25km_Vers6_199108-202601_2026-04-18.nc"},
    "Zhang et al., 2022" : {"file" : SATELLITE_DATA_DIR/ "Zhang et al., 2022/Time_x_y_Surface_Elevation_Anomaly_Greenland_Monthly_5km_Grid.nc"},
    
}