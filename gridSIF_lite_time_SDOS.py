#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import numpy as np
import h5py
import time
from netCDF4 import Dataset

dict_names = {
    'source': 'OCO-2 L1B',
    'references': '...',
    'Conventions': 'CF-1.6',
    'product_version': 'v1.0',
    'summmary': 'Fraunhofer-line based SIF retrievals',
    'keyword': 'satellite, OCO-2, Solar Induced Fluorescence, SIF',
    'keywords_vocabulary': 'NASA Global Change Master Directory (GCMD)',
    'cdm_data_type': 'grids',
    'comment': 'These data were produced at JPL/Caltech',
    'date_created': 'Created ' + time.ctime(time.time()),
    'creator_name': 'Caltech, Christian Frankenberg',
    'creator_email': 'cfranken@caltech.edu',
    'project': 'OCO-2 NASA/JPL',
    'geospatial_lat_min': '-90.0f; // float',
    'geospatial_lat_max': '90.0f; // float',
    'geospatial_lat_units': 'degrees_north',
    'geospatial_lon_min': '-180.0f; // float',
    'geospatial_lon_max': '180.0f; // float',
    'geospatial_lon_units': 'degrees_east',
    'geospatial_vertical_min': '0.0f; // float',
    'geospatial_vertical_max': '100000.0; // float',
    'standard_name_vocabulary': 'NetCDF Climate and Forecast (CF) Metadata Conventions Version 1.6',
    'platform': 'OCO-2',
    'sensor': 'OCO-2',
    'spatial_resolution': '2km x 1.3km at nadir (along-track x across-track)',
    '                   CoordSysBuilder': 'ucar.nc2.dataset.conv.CF1Convention',
    }

dict_l2 = {
    'sif_757nm': 'SIF_757nm',
    'sif_771nm': 'SIF_771nm',
    'sif_757_1sigma': 'SIF_757nm_uncert',
    'sif_771_1sigma': 'SIF_771nm_uncert',
    'latitude': 'latitude',
    'longitude': 'longitude',
    'sif_757_rel': 'SIF_757nm_relative',
    'sif_771_rel': 'SIF_771nm_relative',
    'cont_radiance_757': 'continuum_radiance_757nm',
    'cont_radiance_771': 'continuum_radiance_771nm',
    'sza': 'solar_zenith_angle',
    'lza': 'sensor_zenith_angle',
    'co2_ratio': 'Cloud/co2_ratio',
    'o2_ratio': 'Cloud/o2_ratio',
    'daily_corr': 'daily_correction_factor',
    'T2m': 'Meteo/2m_temperature',
    'VPD_2m': 'Meteo/vapor_pressure_deficit',
    'Tskin': 'Meteo/skin_temperature',
    'altitude': 'surface_altitude',
    }


def grid(options, args):

    b = np.loadtxt(options.orbits)

    dlat = np.arange(options.latMin, options.latMax + 1e-8,
                     options.dLat)
    dlon = np.arange(options.lonMin, options.lonMax + 1e-8,
                     options.dLon)
    lat = np.arange(options.latMin + options.dLat / 2., (options.latMax
                     + 1e-8) - options.dLat / 2., options.dLat)
    lon = np.arange(options.lonMin + options.dLon / 2., (options.lonMax
                     + 1e-8) - options.dLat / 2., options.dLon)

    print b.shape
    nT = b.shape[0]

    for i in dict_l2:
        cmd = 'vec_' + i + '=np.zeros((nT,len(lat),len(lon)))'
        try:
            exec cmd
        except:
            print 'Error executing ' + cmd
    vec_n = np.zeros((nT, len(lat), len(lon)))

    f = Dataset(options.outFile, 'w', format='NETCDF4')
    for i in dict_names:
        cmd = 'f.' + i + '="' + dict_names[i] + '"'
        try:
            exec cmd
        except:
            print 'Error executing ' + cmd

    time = f.createDimension('time', None)
    lati = f.createDimension('lat', len(lat))
    loni = f.createDimension('lon', len(lon))
    times = f.createVariable('time', 'f8', ('time', ))
    latitudes = f.createVariable('lat', 'f8', 'lat')
    longitudes = f.createVariable('lon', 'f8', 'lon')

    times.units = 'days since 2014-1-1 0:0:0'
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    latitudes.standard_name = 'latitude'
    longitudes.standard_name = 'longitude'
    latitudes.axis = 'Y'
    longitudes.axis = 'X'
    times.long_name = 'time'
    latitudes.long_name = 'latitude'
    longitudes.long_name = 'longitude'
    latitudes[:] = lat
    longitudes[:] = lon

    for i in dict_l2:
        cmd = i + '=f.createVariable("' + i\
             + '","f4",("time","lat","lon",), zlib=True,least_significant_digit=4)'
        try:
            exec cmd
        except:
            print 'Error executing ' + cmd
    # amount of soundings 
    n = f.createVariable('N', 'f4', ('time', 'lat', 'lon'))
    n[:, :, :] = 0
    for i in dict_l2:
        cmd = i + '[:,:,:]=0'
        exec cmd
        cmd = i + '.missing_value=-999.9'
        exec cmd

    counter = 0
    counter_time = 0
    for file in args:
        file_split = os.path.basename(file).split('_')

        day = file_split[-3]
        print file,
        fin = h5py.File(file, 'r')
        lat_in = (fin['latitude'])[:]
        lon_in = (fin['longitude'])[:]
        mode = (fin['measurement_mode'])[:]
        fp = (fin['footprint'])[:]
        biome = (fin['IGBP_index'])[:]
        o2_in = (fin['Cloud/o2_ratio'])[:]
        rad = (fin['continuum_radiance_757nm'])[:]
        s = (fin['SIF_757nm'])[:]
        times_in = (fin['time'])[:]
        orb = (fin['orbit_number'])[:]
        if options.biome == -1:
            ind = np.where((s > -5) & (s < 5) & (lat_in
                            > options.latMin) & (lat_in
                            < options.latMax) & (lon_in
                            < options.lonMax) & (lon_in
                            > options.lonMin) & (o2_in > options.o2Min)
                            & (o2_in < options.o2Max) & ((mode
                            == options.mode) | (options.mode == -1))
                            & ((fp == options.footprint)
                            | (options.footprint == -1)))[0]
        else:
            ind = np.where((s > -5) & (s < 5) & (lat_in
                            > options.latMin) & (lat_in
                            < options.latMax) & (lon_in
                            < options.lonMax) & (lon_ind
                            > options.lonMin) & (biome == options.biome)
                            & ((mode == options.mode) | (options.mode
                            == -1)) & ((fp == options.footprint)
                            | (options.footprint == -1)))[0]
        print ' ', len(ind)
        if len(ind) > 0:

            iLat = ((np.asarray(lat_in[ind]) - options.latMin)
                     / (options.latMax - options.latMin)) * len(lat)
            iLon = ((np.asarray(lon_in[ind]) - options.lonMin)
                     / (options.lonMax - options.lonMin)) * len(lon)

            index_vector = np.asarray(iLon * len(lat) + iLat, dtype=int)
            for i in dict_l2:
                cmd = i + '_in=fin["' + dict_l2[i] + '"][:]'
                exec cmd

            for j in range(len(ind)):
                iT = np.where((b[:, 1] >= orb[ind[j]]) & (b[:, 0]
                               <= orb[ind[j]]))[0]
                for i in dict_l2:
                    cmd = 'vec_' + i + '[iT,iLat[j],iLon[j]]+=' + i\
                         + '_in[ind[j]]'
                    exec cmd
                    #print iT, cmd

                vec_n[iT, iLat[j], iLon[j]] += 1
        fin.close()
    for counter in range(nT):
        for i in dict_l2:
            cmd = i + '[counter,:,:] = -999.9'
            exec cmd
    for counter in range(nT):
        wo = np.asarray(np.where(vec_n[counter, :, :] >= 2.))

        if len(wo[1]) > 0:
            for j in range(len(wo[0])):
                for i in dict_l2:
                    cmd = i + '[counter,wo[0][j],wo[1][j]] = vec_'+i+'[counter,wo[0][j],wo[1][j]]/vec_n[counter,wo[0][j],wo[1][j]]'
                    exec cmd
            print counter
            #print vec_sif_757nm[counter, wo[0], wo[1]] / vec_n[counter,
            #        wo[0], wo[1]]
            #print vec_daily_corr[counter, wo[0], wo[1]]\
            #     // vec_n[counter, wo[0], wo[1]]
            #print sif_757nm[counter, wo[0], wo[1]]
            #print daily_corr[counter, wo[0], wo[1]]
            
        #    print sif_757nm[counter, wo[0], wo[1]]
        else:
            print counter, sif_757nm[counter, 0, 0]
        # Had to change this, indexing somehow didn't work anymore as befores
        
        times[counter] = b[counter, 2.]
        n[counter, :, :] = vec_n[counter, :, :]
   # for i in dict_l2:
   #     cmd = i + ' = vec_' + i + '/vec_n'
   #     print cmd
   #     exec cmd

    f.close()


def standalone_main():
    parser = OptionParser(usage='usage: %prog l2_file2')
    parser.add_option('-o', '--outFile', dest='outFile',
                      default='SIFGridOCO2.nc',
                      help='output filename (default SIFGridOCO2.nc)')
    parser.add_option('--latMin', dest='latMin', type=float,
                      default=-90, help='min latitude region')
    parser.add_option('--dLat', dest='dLat', type=float, default=1,
                      help='latitude resolution (1 degree default)')
    parser.add_option('--dLon', dest='dLon', type=float, default=1,
                      help='longitude resolution (1 degree default)')
    
    parser.add_option('--orbits', dest='orbits',
                      default='/home/cfranken/code/IMAP_DOAS/trunk/tools/orbitRange.dat'
                      , help='file with orbit separation')

    parser.add_option('--latMax', dest='latMax', type=float,
                      default=90, help='max latitude region')
    parser.add_option('--lonMin', dest='lonMin', type=float,
                      default=-180, help='min longitude region')
    parser.add_option('--lonMax', dest='lonMax', type=float,
                      default=180, help='max longitude region')
    parser.add_option('--o2Min', dest='o2Min', type=float, default=0,
                      help='min o2 ratio')
    parser.add_option('--o2Max', dest='o2Max', type=float, default=1.1,
                      help='max o2 ratio')
    parser.add_option('--mode', dest='mode', type=int, default=-1,
                      help='mode (0=ND, 1=GL, 2=TG, etc, -1 for all)')
    parser.add_option('--biome', dest='biome', type=int, default=-1,
                      help='IGBP biome type (-1 default for ALL)')
    parser.add_option('--footprint', dest='footprint', type=int,
                      default=-1,
                      help='Which OCO2 footprint to use (1-8) (-1 default for ALL together)'
                      )

    (options, args) = parser.parse_args()

    grid(options, args)


if __name__ == '__main__':
    standalone_main()

