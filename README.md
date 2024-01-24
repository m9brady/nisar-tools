# nisar-tools

Utilities for reading in simulated NISAR data products from NASA-JPL:
https://nisar.jpl.nasa.gov/data/sample-data/

Currently only has rudimentary support for L1 RSLC and L2 GSLC/GCOV products.

### Environment setup

```console
conda create -n nisar-tools python=3 h5py rasterio pyproj shapely -c conda-forge
conda activate nisar-tools
```

### Usage examples
[NISAR L2 GCOV](notebooks/gcov.ipynb)

### RSLC Example
```python
>>> from NISAR import RSLC
# lazily-loaded so no data is immediately read, just metadata
>>> ds = RSLC('NISAR_L1_PR_RSLC_001_005_A_219_2005_DHDH_A_20081127T060959_20081127T061015_P01101_F_N_J_001.h5')
>>> print(ds)
NISAR LSAR RSLC v0.1.0 2008-11-27T06:09:59 orbit #0
>>> print(ds.meta.keys())
dict_keys(['attitude', 'calibrationInformation', 'orbit', 'processingInformation', 'geolocation', 'identification'])
# the image footprint is also determined on read
>>> ds.footprint
<POLYGON Z ((-118.228 34.243 661, -118.144 34.257 661, -118.063 34.271 661, ...>
# loading complex data is straightforward
>>> complex_array = ds.load_data()
>>> print(complex_array.keys())
dict_keys(['HH'])
# it is stored as a structured array instead of complex for some reason
>>> print(complex_array['HH'].dtype)
[('r', '<f2'), ('i', '<f2')]
# method for converting complex to amplitude
>>> amplitude_array = ds.to_amplitude()
>>> print(f'Amplitude from SLC. Shape: {amp.shape} DType: {amp.dtype}')
Amplitude from SLC. Shape: (30540, 6488) DType: float16
```

### GSLC Example
```python
>>> from NISAR import GSLC
# lazily-loaded so no data is immediately read, just metadata
>>> ds = GSLC('NISAR_L2_PR_GSLC_001_005_A_219_2005_DHDH_A_20081127T060959_20081127T061015_P01101_F_N_J_001.h5')
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