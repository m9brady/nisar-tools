import logging
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pyproj
import rasterio as rio
from shapely.wkt import loads

LOGGER = logging.getLogger(__name__)

class NISAR(ABC):
    '''
    Base NISAR product class, extended by various L0/L1/L2 product classes
    '''
    def __init__(self, file_path):
        if not isinstance(file_path, Path):
            file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError('No such file: %s' % file_path.resolve())
        self.file_path = file_path
        self.granule = file_path.name

    @abstractmethod
    def _load_meta(self):
        pass

    @abstractmethod 
    def load_data(self):
        pass
    
class RSLC(NISAR):
    '''
    Subclass for Level-1 Focused SAR product in range/azimuth coordinates (RSLC)
    '''
    def __init__(self, file_path: str | Path):
        super().__init__(file_path)
        self.__meta = self._load_meta(self.file_path)
        self.footprint = loads(self.meta['identification']['boundingPolygon'])
        self.crs = self._get_projection()

    def __str__(self):
        mission = self.meta['identification']['missionId']
        instrument = self.meta['identification']['instrumentName']
        product_type = self.meta['identification']['productType']
        product_version = self.meta['identification']['productVersion']
        abs_orbit_num = self.meta['identification']['absoluteOrbitNumber']
        zero_doppler_start = np.datetime_as_string(self.meta['identification']['zeroDopplerStartTime'])
        return f'{mission} {instrument} {product_type} v{product_version} {zero_doppler_start} orbit #{abs_orbit_num}'
        
    @property
    def meta(self):
        return self.__meta
    
    @meta.setter
    def meta(self, *args, **kwargs):
        LOGGER.warning('Overwriting of metadata property not permitted')
        return
    
    def _load_meta(self, granule):
        bool_map = {"true": True, "false": False}
        with h5py.File(granule, mode='r') as ds:
            identification = {k: v[()] for k,v in ds['science/LSAR/identification'].items()}
            for item, ident_meta in identification.items():
                if isinstance(ident_meta, np.bytes_):
                    ident_meta_str = ident_meta.decode('utf-8')
                    # boolean True/False
                    if item.startswith('is'):
                        identification[item] = bool_map.get(ident_meta_str.lower())
                    # time-aware
                    if 'Time' in item:
                        identification[item] = np.datetime64(ident_meta_str)
                    # otherwise, assume string stored as np.bytes
                    else:
                        identification[item] = ident_meta_str
            attitude = {k: v[()] for k,v in ds['science/LSAR/RSLC/metadata/attitude'].items()}
            for item, attitude_meta in attitude.items():
                if isinstance(attitude_meta, np.bytes_):
                    attitude_meta_str = attitude_meta.decode('utf-8')
                    attitude[item] = attitude_meta_str
            calibration = {}
            for freq_suffix in identification['listOfFrequencies']:
                freq = f'frequency{freq_suffix.decode("utf-8")}'
                for pol in ds[f'science/LSAR/RSLC/swaths/{freq}/listOfPolarizations']:
                    pol_str = pol.decode('utf-8')
                    pol_meta = ds[f'science/LSAR/RSLC/metadata/calibrationInformation/{freq}/{pol_str}']
                    calibration[pol_str] = {k: v[()] for k,v in pol_meta.items()}
                    calibration[pol_str].update({
                        'elevationAntennaPattern':  ds[f'science/LSAR/RSLC/metadata/calibrationInformation/{freq}/elevationAntennaPattern/{pol_str}'][()],
                        'noiseEquivalentBackscatter': ds[f'science/LSAR/RSLC/metadata/calibrationInformation/{freq}/nes0/{pol_str}'][()]
                    })
                calibration[pol_str].update({
                    'commonDelay': ds[f'science/LSAR/RSLC/metadata/calibrationInformation/{freq}/commonDelay'][()],
                    'faradayRotation': ds[f'science/LSAR/RSLC/metadata/calibrationInformation/{freq}/faradayRotation'][()]
                })
            for k, v in ds['science/LSAR/RSLC/metadata/calibrationInformation'].items():
                # frequencyN
                if k.startswith('freq') and isinstance(v, h5py.Group):
                    for pol in v:
                        if isinstance(v[pol], h5py.Dataset):
                            calibration[pol] = v[pol][()]
                        else:
                            calibration[pol] = {kk: vv[()] for kk,vv in v[pol].items()}
                elif isinstance(v, h5py.Dataset):
                    calibration[k] = v[()]
            orbit = {k: v[()] for k,v in ds['science/LSAR/RSLC/metadata/orbit'].items()}
            for item, orbit_meta in orbit.items():
                if isinstance(orbit_meta, np.bytes_):
                    orbit[item] = orbit_meta.decode('utf-8')
            # This is really ugly but I couldn't figure out h5py.Dataset.visit
            processing = {}
            for item, proc_meta in ds['science/LSAR/RSLC/metadata/processingInformation'].items():
                processing[item] = {}
                for k, v in proc_meta.items():
                    if not isinstance(v, h5py.Group):
                        v_val = v[()]
                        if isinstance(v_val, np.bytes_):
                            v_val_str = v_val.decode('utf-8')
                            processing[item][k] = v_val_str
                        else:
                            processing[item][k] = v_val
                    else:
                        processing[item][k] = {}
                        for kk, vv in v.items():
                            if not isinstance(vv, h5py.Group):
                                vv_val = vv[()]
                                if isinstance(vv_val, np.bytes_):
                                    vv_val_str = vv_val.decode('utf-8')
                                    processing[item][k][kk] = vv_val_str
                                else:
                                    processing[item][k][kk] = vv_val
                            else:
                                processing[item][k][kk] = {kkk: vvv[()] for kkk, vvv in vv.items()}
            geolocation = {}
            for item, geo_meta in ds['science/LSAR/RSLC/metadata/geolocationGrid'].items():
                if item == 'projection':
                    geolocation[item] = {k: v for k,v in geo_meta.attrs.items()}
                else:
                    geolocation[item] = geo_meta[()]
        return {
            'attitude': attitude,
            'calibrationInformation': calibration,
            'orbit': orbit,
            'processingInformation': processing,
            'geolocation': geolocation,
            'identification': identification
        }

    def load_data(self, polarisation: str='ALL') -> dict:
        """
        Loads the complex image array(s) from the HDF5 file as numpy array(s)
        """
        data = {}
        with h5py.File(self.file_path, mode='r') as ds:
            # Case A: all polarisations wanted
            if polarisation.upper() == 'ALL':
                for freq in self.meta['identification']['listOfFrequencies']:
                    for pol in ds[f'science/LSAR/RSLC/metadata/calibrationInformation/frequency{freq.decode("utf-8")}']:
                        data[pol] = ds[f'science/LSAR/RSLC/swaths/frequency{freq.decode("utf-8")}/{pol}'][()]
            else:
                # Case B: multiple polarisations wanted
                if ',' in polarisation:
                    for freq in self.meta['identification']['listOfFrequencies']:
                        for pol in polarisation.replace(' ', '').split(','):
                            data[pol] = ds[f'science/LSAR/RSLC/swaths/frequency{freq.decode("utf-8")}/{pol}'][()]
                # Case C: single polarisation wanted
                else:
                    data[polarisation] = None
                    for freq in self.meta['identification']['listOfFrequencies']:
                        try:
                            data[polarisation] = ds[f'science/LSAR/RSLC/swaths/frequency{freq.decode("utf-8")}/{polarisation}'][()]
                        except KeyError:
                            pass
                    if data[polarisation] is None:
                        raise ValueError('Cannot locate polarisation %r in file %s' % (polarisation, self.file_path.resolve()))
        return data
    
    def _get_projection(self) -> pyproj.CRS:
        return pyproj.CRS.from_epsg(self.meta['geolocation']['epsg'])
    
    def to_amplitude(self, polarisation: str='HH') -> np.ndarray:
        """
        Convert RSLC to Amplitude array for fun and enjoyment
        """
        # convert SLC to amplitude
        img = self.load_data(polarisation)[polarisation]
        # return amplitude
        return np.abs(img)

class GSLC(NISAR):
    '''
    Subclass for Level-2 Focused SAR product in geocoded coordinates (GSLC)
    '''
    def __init__(self, file_path: str | Path):
        super().__init__(file_path)
        self.__meta = self._load_meta(self.file_path)
        self.footprint = loads(self.meta['identification']['boundingPolygon'])
        self.crs = self._get_projection()
        self.available_polarisations = list(self.meta['calibrationInformation'].keys())

    def __str__(self):
        mission = self.meta['identification']['missionId']
        instrument = self.meta['identification']['instrumentName']
        product_type = self.meta['identification']['productType']
        product_version = self.meta['identification']['productVersion']
        abs_orbit_num = self.meta['identification']['absoluteOrbitNumber']
        zero_doppler_start = np.datetime_as_string(self.meta['identification']['zeroDopplerStartTime'])
        return f'{mission} {instrument} {product_type} v{product_version} {zero_doppler_start} orbit #{abs_orbit_num}'

    @property
    def meta(self):
        return self.__meta
    
    @meta.setter
    def meta(self, *args, **kwargs):
        LOGGER.warning('Overwriting of metadata property not permitted')
        return

    def _load_meta(self, granule: str | Path) -> dict:
        with h5py.File(granule, mode='r') as ds:
            identification = {k: v[()] for k,v in ds['science/LSAR/identification'].items()}
            attitude = {k: v[()] for k,v in ds['science/LSAR/GSLC/metadata/attitude'].items()}
            calibration = {}
            for freq in ds['science/LSAR/GSLC/metadata/calibrationInformation']:
                for pol, pol_meta in ds[f'science/LSAR/GSLC/metadata/calibrationInformation/{freq}'].items():
                    calibration[pol] = {k: v[()] for k,v in pol_meta.items()}
            orbit = {k: v[()] for k,v in ds['science/LSAR/GSLC/metadata/orbit'].items()}
            processing = {}
            for item, proc_meta in ds['science/LSAR/GSLC/metadata/processingInformation'].items():
                processing[item] = {}
                for k, v in proc_meta.items():
                    if not isinstance(v, h5py.Group):
                        processing[item][k] = v[()]
                    else:
                        processing[item][k] = {kk: vv[()] for kk, vv in v.items()}
            radar = {}
            for item, radar_meta in ds['science/LSAR/GSLC/metadata/radarGrid'].items():
                if item == 'projection':
                    radar[item] = {k: v for k,v in radar_meta.attrs.items()}
                else:
                    radar[item] = radar_meta[()]
        #TODO: consider adding unit information to relevant metadata arrays (e.g. all the timedelta arrays)
        return {
            'attitude': attitude,
            'calibrationInformation': calibration,
            'orbit': orbit,
            'processingInformation': processing,
            'radar': radar,
            'identification': identification
        }
    
    def load_data(self, polarisation: str='ALL') -> dict:
        """
        Loads the complex image array(s) from the HDF5 file as numpy array(s)
        """
        data = {}
        with h5py.File(self.file_path, mode='r') as ds:
            # Case A: all polarisations wanted
            if polarisation.upper() == 'ALL':
                for freq in ds['science/LSAR/GSLC/grids']:
                    for pol in self.available_polarisations:
                        data[pol] = ds[f'science/LSAR/GSLC/grids/{freq}/{pol}'][()]
            else:
                # Case B: multiple polarisations wanted
                if ',' in polarisation:
                    for freq in ds['science/LSAR/GSLC/grids']:
                        for pol in polarisation.replace(' ', '').split(','):
                            data[pol] = ds[f'science/LSAR/GSLC/grids/{freq}/{pol}'][()]
                # Case C: single polarisation wanted
                else:
                    data[polarisation] = None
                    for freq in ds['science/LSAR/GSLC/grids']:
                        try:
                            data[polarisation] = ds[f'science/LSAR/GSLC/grids/{freq}/{polarisation}'][()]
                        except KeyError:
                            pass
                    if data[polarisation] is None:
                        raise ValueError('Cannot locate polarisation %r in file %s' % (polarisation, self.file_path.resolve()))
        return data

    def _get_projection(self) -> pyproj.CRS:
        return pyproj.CRS.from_wkt(
            self.meta['radar']['projection']['spatial_ref'].decode('utf-8')
        )

    def _load_xy_coords(self) -> tuple[np.ndarray, np.ndarray]:
        '''
        This is dumb because it reads in the xy arrays for each freq
        Maybe just check for FrequencyA and be done with it
        '''
        with h5py.File(self.file_path, mode='r') as ds:
            xs = None
            ys = None
            for freq in self.meta['identification']['listOfFrequencies']:
                try:
                    xs = ds[f'science/LSAR/GSLC/grids/frequency{freq.decode("utf-8")}/xCoordinates'][()]
                    ys = ds[f'science/LSAR/GSLC/grids/frequency{freq.decode("utf-8")}/yCoordinates'][()]
                except KeyError:
                    pass
        if xs is None or ys is None:
            raise ValueError('Cannot locate image X/Y coords in file %s' % self.file_path.resolve())
        return xs, ys

    def to_geotiff(self, target_file: str | Path, polarisation: str='HH') -> Path:
        """
        Convert GSLC to Amplitude GeoTIFF for fun and enjoyment
        """
        if not isinstance(target_file, Path):
            target_file = Path(target_file)
        if target_file.exists():
            LOGGER.warning('File exists: %r' % target_file.resolve())
            return target_file
        # convert SLC to amplitude
        img = np.abs(self.load_data(polarisation)[polarisation])
        height, width = img.shape
        xs, ys = self._load_xy_coords()
        xmin, xmax = np.percentile(xs, (0, 100))
        ymin, ymax = np.percentile(ys, (0, 100))
        transform = rio.transform.from_bounds(xmin, ymin, xmax, ymax, width, height)
        tiff_meta = {
            'driver': 'GTiff',
            'width': width,
            'height': height,
            'count': 1,
            'dtype': img.dtype,
            'crs': self.crs,
            'transform': transform,
            'nodata': np.nan,
            'compress': 'LZW',
            'tiled': True,
            'blockxsize': 256,
            'blockysize': 256,
        }
        with rio.open(target_file, 'w', **tiff_meta) as ds:
            ds.write(img, 1)
        return target_file
    
class GCOV(NISAR):
    '''
    Subclass for Level-2 SAR covariance product in geocoded coordinates (GCOV)
    '''
    def __init__(self, file_path: str | Path):
        super().__init__(file_path)
        self.__meta = self._load_meta(self.file_path)
        self.footprint = loads(self.meta['identification']['boundingPolygon'])
        self.crs = self._get_projection()
        self.available_polarisations = list(self.meta['calibrationInformation'].keys())

    def __str__(self):
        mission = self.meta['identification']['missionId']
        instrument = self.meta['identification']['instrumentName']
        product_type = self.meta['identification']['productType']
        product_version = self.meta['identification']['productVersion']
        abs_orbit_num = self.meta['identification']['absoluteOrbitNumber']
        zero_doppler_start = np.datetime_as_string(self.meta['identification']['zeroDopplerStartTime'])
        return f'{mission} {instrument} {product_type} v{product_version} {zero_doppler_start} orbit #{abs_orbit_num}'

    @property
    def meta(self):
        return self.__meta
    
    @meta.setter
    def meta(self, *args, **kwargs):
        LOGGER.warning('Overwriting of metadata property not permitted')
        return

    def _load_meta(self, granule: str | Path) -> dict:
        bool_map = {"true": True, "false": False}
        with h5py.File(granule, mode='r') as ds:
            identification = {k: v[()] for k,v in ds['science/LSAR/identification'].items()}
            for item, ident_meta in identification.items():
                if isinstance(ident_meta, np.bytes_):
                    ident_meta_str = ident_meta.decode('utf-8')
                    # boolean True/False
                    if item.startswith('is'):
                        identification[item] = bool_map.get(ident_meta_str.lower())
                    # time-aware
                    if 'Time' in item:
                        identification[item] = np.datetime64(ident_meta_str)
                    # otherwise, assume string stored as np.bytes
                    else:
                        identification[item] = ident_meta_str
            attitude = {k: v[()] for k,v in ds['science/LSAR/GCOV/metadata/attitude'].items()}
            for item, attitude_meta in attitude.items():
                if isinstance(attitude_meta, np.bytes_):
                    attitude_meta_str = attitude_meta.decode('utf-8')
                    attitude[item] = attitude_meta_str
            calibration = {}
            for freq in ds['science/LSAR/GCOV/grids']:
                for pol in ds[f'science/LSAR/GCOV/grids/{freq}/listOfPolarizations']:
                    pol_str = pol.decode('utf-8')
                    pol_meta = ds[f'science/LSAR/GCOV/metadata/calibrationInformation/{freq}/{pol_str}']
                    calibration[pol_str] = {k: v[()] for k,v in pol_meta.items()}
                    calibration[pol_str].update({
                        'elevationAntennaPattern':  ds[f'science/LSAR/GCOV/metadata/calibrationInformation/{freq}/elevationAntennaPattern/{pol_str}'][()],
                        'noiseEquivalentBackscatter': ds[f'science/LSAR/GCOV/metadata/calibrationInformation/{freq}/noiseEquivalentBackscatter/{pol_str}'][()]
                    })
                calibration[pol_str].update({
                    'commonDelay': ds[f'science/LSAR/GCOV/metadata/calibrationInformation/{freq}/commonDelay'][()],
                    'faradayRotation': ds[f'science/LSAR/GCOV/metadata/calibrationInformation/{freq}/faradayRotation'][()]
                })
            orbit = {k: v[()] for k,v in ds['science/LSAR/GCOV/metadata/orbit'].items()}
            for item, orbit_meta in orbit.items():
                if isinstance(orbit_meta, np.bytes_):
                    orbit[item] = orbit_meta.decode('utf-8')
            # This is really ugly but I couldn't figure out h5py.Dataset.visit
            processing = {}
            for item, proc_meta in ds['science/LSAR/GCOV/metadata/processingInformation'].items():
                processing[item] = {}
                for k, v in proc_meta.items():
                    if not isinstance(v, h5py.Group):
                        v_val = v[()]
                        if isinstance(v_val, np.bytes_):
                            v_val_str = v_val.decode('utf-8')
                            if k.startswith('is') or k.endswith('Applied'):
                                processing[item][k] = bool_map.get(v_val_str.lower())
                            else:
                                processing[item][k] = v_val_str
                        else:
                            processing[item][k] = v_val
                    else:
                        processing[item][k] = {}
                        for kk, vv in v.items():
                            if not isinstance(vv, h5py.Group):
                                vv_val = vv[()]
                                if isinstance(vv_val, np.bytes_):
                                    vv_val_str = vv_val.decode('utf-8')
                                    if kk.startswith('is') or kk.endswith('Applied'):
                                        processing[item][k][kk] = bool_map.get(vv_val_str.lower())
                                    else:
                                        processing[item][k][kk] = vv_val_str
                                else:
                                    processing[item][k][kk] = vv_val
                            else:
                                processing[item][k][kk] = {kkk: vvv[()] for kkk, vvv in vv.items()}
            radar = {}
            for item, radar_meta in ds['science/LSAR/GCOV/metadata/radarGrid'].items():
                # the projection item has several useful CRS-related attributes
                if item == 'projection':
                    projection = {}
                    for k, v in radar_meta.attrs.items():
                        if isinstance(v, np.bytes_):
                            projection[k] = v.decode('utf-8')
                        else:
                            projection[k] = v
                    radar[item] = projection
                else:
                    if isinstance(radar_meta, np.bytes_):
                        radar[item] = radar_meta.decode('utf-8')
                    else:
                        radar[item] = radar_meta[()]
            ceos_ard = {}
            for item, ceos_meta in ds['science/LSAR/GCOV/metadata/ceosAnalysisReadyData'].items():
                if not isinstance(ceos_meta, h5py.Group):
                    value = ceos_meta[()]
                    if isinstance(value, np.bytes_):
                        ceos_ard[item] = value.decode('utf-8')
                    else:
                        ceos_ard[item] = value
                else:
                    for k, v in ceos_meta.items():
                        ceos_ard[item] = {}
                        ceos_ard[item][k] = {kk: vv[()] for kk, vv in v.items()}
        #TODO: consider adding unit information to relevant metadata arrays (e.g. all the timedelta arrays)
        return {
            'attitude': attitude,
            'calibrationInformation': calibration,
            'orbit': orbit,
            'processingInformation': processing,
            'radar': radar,
            'identification': identification,
            'ceosInformation': ceos_ard
        }
    
    def load_data(self, polarisation: str='ALL') -> dict:
        """
        Loads the complex image array(s) from the HDF5 file as numpy array(s)
        """
        data = {}
        with h5py.File(self.file_path, mode='r') as ds:
            # Case A: all polarisations wanted
            if polarisation.upper() == 'ALL':
                for freq in ds['science/LSAR/GCOV/grids']:
                    for pol in self.available_polarisations:
                        #TODO: confirm that the pol is what the sample products show (e.g. HH is actualy "HHHH")
                        data[pol] = ds[f'science/LSAR/GCOV/grids/{freq}/{pol}{pol}'][()]
            else:
                # Case B: multiple polarisations wanted
                if ',' in polarisation:
                    for freq in ds['science/LSAR/GCOV/grids']:
                        for pol in polarisation.replace(' ', '').split(','):
                            data[pol] = ds[f'science/LSAR/GCOV/grids/{freq}/{pol}{pol}'][()]
                # Case C: single polarisation wanted
                else:
                    data[polarisation] = None
                    for freq in ds['science/LSAR/GCOV/grids']:
                        try:
                            data[polarisation] = ds[f'science/LSAR/GCOV/grids/{freq}/{polarisation}{polarisation}'][()]
                        except KeyError:
                            pass
                    if data[polarisation] is None:
                        raise ValueError('Cannot locate polarisation %r in file %s' % (polarisation, self.file_path.resolve()))
        return data

    def _get_projection(self) -> pyproj.CRS:
        return pyproj.CRS.from_wkt(
            self.meta['radar']['projection']['spatial_ref']
        )

    def _load_xy_coords(self) -> tuple[np.ndarray, np.ndarray]:
        '''
        This is dumb because it reads in the xy arrays for each freq
        Maybe just check for FrequencyA and be done with it
        '''
        with h5py.File(self.file_path, mode='r') as ds:
            xs = None
            ys = None
            for freq in self.meta['identification']['listOfFrequencies']:
                try:
                    xs = ds[f'science/LSAR/GCOV/grids/frequency{freq.decode("utf-8")}/xCoordinates'][()]
                    ys = ds[f'science/LSAR/GCOV/grids/frequency{freq.decode("utf-8")}/yCoordinates'][()]
                except KeyError:
                    pass
        if xs is None or ys is None:
            raise ValueError('Cannot locate image X/Y coords in file %s' % self.file_path.resolve())
        return xs, ys

    def to_geotiff(self, target_file: str | Path, polarisation: str='HH') -> Path:
        """
        Convert GCOV to dB GeoTIFF (COG) for fun and enjoyment
        """
        if not isinstance(target_file, Path):
            target_file = Path(target_file)
        if target_file.exists():
            LOGGER.warning('File exists: %r' % target_file.resolve())
            return target_file
        # load covariance data and convert to dB
        img = 10 * np.log10(self.load_data(polarisation)[polarisation])
        height, width = img.shape
        xs, ys = self._load_xy_coords()
        xmin = xs.min()
        ymax = ys.max()
        xres = abs(xs[1] - xs[0])
        yres = abs(ys[1] - ys[0])
        # since the X and Y arrays stored in the H5 data appear to be grid-cell CENTERS,
        # but we need the CORNERS for computing the affine transform, we get a little
        # creative
        transform = rio.transform.from_origin(
            west=xmin - xres * 0.5,
            north=ymax + yres * 0.5,
            xsize=xres,
            ysize=yres
        )
        tiff_meta = {
            'driver': 'GTiff',
            'width': width,
            'height': height,
            'count': 1,
            'dtype': img.dtype,
            'crs': self.crs,
            'transform': transform,
            'nodata': np.nan,
            'compress': 'LZW',
            'tiled': True,
            'blockxsize': 256,
            'blockysize': 256,
        }
        with rio.open(target_file, 'w', **tiff_meta) as ds:
            ds.write(img, 1)
        return target_file    