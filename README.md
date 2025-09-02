# nisar-tools

Utilities for reading in simulated NISAR data products from NASA-JPL:
https://science.nasa.gov/mission/nisar/sample-data/

Currently only has rudimentary support for L1 RSLC and L2 GSLC/GCOV products.

## Product Table

> Adapted from Table 2.1 in the L1/L2 ATBD. Bolded/italicized products indicate rudimentary support in this repository.

| Product                                   | Level | Description                                                                                                         |
|-------------------------------------------|-------|---------------------------------------------------------------------------------------------------------------------|
| Radar Raw Science Telemetry (L0A) (RRST) | L0A   | This L0A product is the raw downlinked data delivered to SDS. Communication wrapping has been removed                |
| Radar Raw Signal Data (L0B) (RRSD)        | L0B   | This L0B product is corrected, aligned radar pulse data derived from the RRST products and used for further processing. |
| ***Range Doppler Single Look Complex (RSLC)*** | L1    | Focused SAR Imagery in range-Doppler coordinates.                                                                   |
| Range Doppler Interferogram (RIFG)        | L1    | Multi-looked flattened (ellipsoid) interferogram in range-Doppler coordinates with no removal of topographic fringes. Formed using high-res offsets. |
| Range Doppler Unwrapped Interferogram (RUNW) | L1 | Multi-looked unwrapped differential interferogram in range-Doppler coordinates with topo fringes removed.            |
| Range Doppler pixel Offsets (ROFF)         | L1    | Raw pixel offsets layers in range and azimuth directions derived by speckle tracking with different resolutions (e.g., different chip and search window size) in range-Doppler coordinates. |
| ***Geocoded Single Look Complex (GSLC)***        | L2    | Geocoded SLC product using the MOE state vectors and a DEM.                                                        |
| Geocoded Unwrapped Interferogram (GUNW)    | L2    | Geocoded, multi-looked unwrapped differential Interferogram.                                                      |
| ***Geocoded Polarimetric Covariance (GCOV)***    | L2    | Geocoded, multi-looked polarimetric covariance matrix.                                                             |
| Geocoded pixel Offsets (GOFF)               | L2    | Raw pixel offsets in range and azimuth directions derived by speckle tracking with different resolutions (e.g., different chip and search window size) in geocoded coordinates. |

## Environment setup

```console
conda create -n nisar-tools python=3 h5py rasterio pyproj shapely -c conda-forge
conda activate nisar-tools
```

## Usage examples

### Notebook examples
- [NISAR L2 GSLC](notebooks/gslc.ipynb)
- [NISAR L2 GCOV](notebooks/gcov.ipynb)

### REPL examples

## RSLC
```python
>>> from NISAR import RSLC
# lazily-loaded so no data is immediately read, just metadata
>>> ds = RSLC('NISAR_L1_PR_RSLC_001_030_A_019_2000_SHNA_A_20081012T060910_20081012T060926_D00402_N_F_J_001.h5')
# this is sample data so it makes sense to see the sensor it was derived from
>>> print(ds)
ALOS PALSAR RSLC v0.1.0 2008-10-12T06:09:12.000000000 orbit #0
>>> print(ds.meta.keys())
dict_keys(['attitude', 'calibrationInformation', 'orbit', 'processingInformation', 'geolocation', 'identification'])
# the image footprint is also determined on read
>>> ds.footprint
<POLYGON Z ((-118.229 34.348 1613.109, -118.152 34.361 1357.626, -118.063 34...>
# loading complex data is straightforward (default is HH but can be chosen with polarization arg)
>>> complex_array = ds.load_data()
>>> print(complex_array.keys())
dict_keys(['HH'])
# it is stored as a complex64 ndarray
>>> print(complex_array['HH'].dtype)
dtype('complex64')
# method for converting complex to amplitude
>>> amplitude_array = ds.to_amplitude()
>>> print(f'Amplitude from SLC. Shape: {amp.shape} DType: {amp.dtype}')
Amplitude from SLC. Shape: (30540, 6488) DType: float16
```

## GSLC
```python
>>> from NISAR import GSLC
# lazily-loaded so no data is immediately read, just metadata
>>> ds = GSLC('NISAR_L2_PR_GSLC_001_030_A_019_2800_SHNA_A_20081012T060911_20081012T060925_D00404_N_F_J_001.h5')
# this is sample data so it makes sense to see the sensor it was derived from
>>> print(ds)
NISAR LSAR GSLC v0.1.0 2008-11-27T06:09:59 orbit #0
>>> print(ds.meta.keys())
dict_keys(['attitude', 'calibrationInformation', 'orbit', 'processingInformation', 'radar', 'identification'])
# the image footprint is also determined on read
>>> ds.footprint
<POLYGON Z ((-118.228 34.243 661, -118.144 34.257 661, -118.063 34.271 661, ...>
# as well as the projection
>>> ds.crs
<Projected CRS: EPSG:32611>
Name: WGS 84 / UTM zone 11N
Axis Info [cartesian]:
- [east]: Easting (metre)
- [north]: Northing (metre)
Area of Use:
- undefined
Coordinate Operation:
- name: UTM zone 11N
- method: Transverse Mercator
Datum: World Geodetic System 1984
- Ellipsoid: WGS 84
- Prime Meridian: Greenwich
# loading complex data is straightforward
>>> complex_array = ds.load_data()
>>> print(complex_array.keys())
dict_keys(['HH'])
# it is stored as a complex datatype
>>> print(complex_array['HH'].dtype)
complex64
# method to convert to an amplitude geotiff for downstream applications
>>> ds.to_geotiff('tif_data.tif')
```

## Further Reading

[ASF Tutorials for NISAR Data Manipulation](https://www.earthdata.nasa.gov/learn/tutorials/work-nisar-sample-data)

[ATBD for L1/L2 NISAR Products](https://nisar.asf.earthdatacloud.nasa.gov/NISAR-SAMPLE-DATA/DOCS/NISAR_D-95677_NASA_L1_L2_ATBD_20231112_R3.4_w-sigs.pdf)

[ATBD for L3 NISAR Soil Moisture Products](https://nisar.asf.earthdatacloud.nasa.gov/NISAR-SAMPLE-DATA/DOCS/NISAR_D-107679_L3SM_ATBD_R3.3_20230428_w-sigs.pdf)
