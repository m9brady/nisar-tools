import logging
from abc import ABC, abstractmethod
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import h5py
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
    def __init__(self, file_path: (str, Path)):
        super().__init__(file_path)
        self.__meta = self._load_meta(self.file_path)
        self.footprint = loads(self.meta['identification']['boundingPolygon'].decode('utf-8'))
        self.crs = self._get_projection()

    def __str__(self):
        mission = self.meta['identification']['missionId'].decode('utf-8')
        instrument = self.meta['identification']['instrumentName'].decode('utf-8')
        product_type = self.meta['identification']['productType'].decode('utf-8')
        product_version = self.meta['identification']['productVersion'].decode('utf-8')
        abs_orbit_num = self.meta['identification']['absoluteOrbitNumber']
        zero_doppler_start = datetime.strptime(
            self.meta['identification']['zeroDopplerStartTime'].decode('utf-8')[:-3],
            '%Y-%m-%dT%H:%M:%S.%f'
        ).strftime('%Y-%m-%dT%H:%M:%S')
        return f'{mission} {instrument} {product_type} v{product_version} {zero_doppler_start} orbit #{abs_orbit_num}'
        
    @property
    def meta(self):
        return self.__meta
    
    @meta.setter
    def meta(self, *args, **kwargs):
        LOGGER.warning('Overwriting of metadata property not permitted')
        return
    
    def _load_meta(self, granule):
        with h5py.File(granule, mode='r') as ds:
            identification = {k: v[()] for k,v in ds['science/LSAR/identification'].items()}
            attitude = {k: v[()] for k,v in ds['science/LSAR/RSLC/metadata/attitude'].items()}
            calibration = {}
            for k, v in ds['science/LSAR/RSLC/metadata/calibrationInformation'].items():
                # frequencyN
                if k.startswith('freq') and isinstance(v, h5py.Group):
                    for pol in v:
                        calibration[pol] = {kk: vv[()] for kk,vv in v[pol].items()}
                elif isinstance(v, h5py.Dataset):
                    calibration[k] = v[()]
            orbit = {k: v[()] for k,v in ds['science/LSAR/RSLC/metadata/orbit'].items()}
            processing = {}
            for item, proc_meta in ds['science/LSAR/RSLC/metadata/processingInformation'].items():
                processing[item] = {}
                for k, v in proc_meta.items():
                    if not isinstance(v, h5py.Group):
                        processing[item][k] = v[()]
                    else:
                        processing[item][k] = {kk: vv[()] for kk, vv in v.items()}
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

class GSLC(NISAR):
    '''
    Subclass for Level-2 Focused SAR product in geocoded coordinates (GSLC)
    '''
    def __init__(self, file_path: (str, Path)):
        super().__init__(file_path)
        self.__meta = self._load_meta(self.file_path)
        self.footprint = loads(self.meta['identification']['boundingPolygon'].decode('utf-8'))
        self.crs = self._get_projection()
        self.available_polarisations = list(self.meta['calibrationInformation'].keys())

    def __str__(self):
        mission = self.meta['identification']['missionId'].decode('utf-8')
        instrument = self.meta['identification']['instrumentName'].decode('utf-8')
        product_type = self.meta['identification']['productType'].decode('utf-8')
        product_version = self.meta['identification']['productVersion'].decode('utf-8')
        abs_orbit_num = self.meta['identification']['absoluteOrbitNumber']
        zero_doppler_start = datetime.strptime(
            self.meta['identification']['zeroDopplerStartTime'].decode('utf-8')[:-3],
            '%Y-%m-%dT%H:%M:%S.%f'
        ).strftime('%Y-%m-%dT%H:%M:%S')
        return f'{mission} {instrument} {product_type} v{product_version} {zero_doppler_start} orbit #{abs_orbit_num}'

    @property
    def meta(self):
        return self.__meta
    
    @meta.setter
    def meta(self, *args, **kwargs):
        LOGGER.warning('Overwriting of metadata property not permitted')
        return

    def _load_meta(self, granule: (str, Path)) -> dict:
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

    def _load_xy_coords(self) -> (np.ndarray, np.ndarray):
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

    def to_geotiff(self, target_file: (str, Path), polarisation: str='HH') -> Path:
        """
        Convert GSLC to Amplitude GeoTIFF for fun and enjoyment

        This might be broken because the sample data I have doesn't seem to align with openstreetmap
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