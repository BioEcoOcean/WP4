This script estimates Annual Net Community Production (ANCP) from BGC-Argo nitrate profiles and produces diagnostics (time series, maps, histograms, and time–depth evolution plots).

Method (summary):

Load & group NetCDF profiles by float ID.
Variable aliases supported:

nitrate: NITRATE_ADJUSTED / NTAW_ADJUSTED

pressure: PRES / PRESSURE

chlorophyll (optional): CHLA_FLUORESCENCE_ADJUSTED / CHLA_ADJUSTED

time: JULD / TIME; position: LATITUDE, LONGITUDE

For each cast, interpolate nitrate (and chlorophyll) to a common 0–150 m grid (1 m). Pressure is made strictly monotonic (unique) prior to interpolation.

Build time series of depth-integrated nitrate and chlorophyll (0–150 m).

Per float and year, select a winter maximum (Dec of previous year + Jan–Mar) and a late-summer minimum (Aug–Nov) in the integrated nitrate.

Compute the drawdown profile (max − min), convert nitrate units:
µmol/kg → mol/m^3 via density (ρ = 1025 kg/m^3), integrate over depth (mol/m^2), then convert to gC/m^2 using molar mass of carbon (12 g/mol) and Redfield C:N (106/16).

Diagnostics:

ANCP vs integrated chlorophyll (linear + log-log),

Map of ANCP at the minimum-nitrate profile locations,

Histogram with example literature ranges annotated,

Time–depth heatmaps for nitrate and chlorophyll.

Inputs:
A folder argo_data_nitrate/ with one Argo profile per NetCDF file.

Outputs:
Arrays of ANCP, locations, years, float IDs, integrated chlorophyll at the minimum, and figures for QA.

Assumptions & caveats:

Seasonal windows are heuristic and may need tuning for region-specific phenology.

Chlorophyll pairing chooses the column that yields the most valid depths within 0–150 m when interpolated.

Map plotting uses geoscatter/geobasemap; basemap downloads may require internet access.

The histogram’s annotated ranges (“Norwegian Sea”, “Siberian Shelf”) are illustrative and can be edited.

Adapting:

Change common_depth for other vertical domains.

Adjust seasonal windows for different ecosystems.

Replace constant density (ρ) with in-situ density if available for more precise units.
