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
        self.granule = file_path
        
    @abstractmethod 
    def load_data(self):
        pass
    
    @abstractmethod
    def _load_meta(self):
        pass

class GLSC(NISAR):
    '''
    Subclass for Level-2 Focused SAR produce in geocoded coordinates (GSLC)
    '''
    def __init__(self, file_path):
        super().__init__(file_path)
        self.__meta = self._load_meta(self.granule)
        self.footprint = loads(self.meta['identification']['boundingPolygon'].decode('utf-8'))
        self.crs = self._get_projection()
        self.available_polarisations = list(self.meta['calibrationInformation'].keys())

    def __repr__(self):
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
    
    def load_data(self, polarisation='ALL'):
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

    def _get_projection(self):
        return pyproj.CRS.from_wkt(
            self.meta['radar']['projection']['spatial_ref'].decode('utf-8')
        )

    def to_geotiff(self, target_file, polarisation='HH'):
        """
        Convert GLSC to Amplitude GeoTIFF for fun and enjoyment

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
        # this might be the problem since xCoordinates and 
        # yCoordinates aren't the same shape as the image array
        xmin, xmax = np.percentile(self.meta['radar']['xCoordinates'], (0, 100))
        ymin, ymax = np.percentile(self.meta['radar']['yCoordinates'], (0, 100))
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