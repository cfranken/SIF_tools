#!/usr/bin/env python
#
# Copyright 2018 California Institute of Technology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Author: Christian Frankenberg, cfranken@caltech.edu

import numpy as np
import h5py
#from netCDF4 import Dataset
import glob
import traceback


# This dictionary basically determines which data to read in (this can be changed!):
dict_oco2 = {
   #'variable_name': 'standard_name_from_netcdf',
    'lat': 'latitude',
    'lon': 'longitude',
    'sza': 'solar_zenith_angle',
    'vza': 'sensor_zenith_angle',
    'saa': 'solar_azimuth_angle',
    'vaa': 'sensor_azimuth_angle',
    'biome': 'IGBP_index',
    'mode': 'measurement_mode',
    'sif_757': 'SIF_757nm',
    'sif_771': 'SIF_771nm',
    'sif_757_sigma': 'SIF_757nm_uncert',
    'sif_771_sigma': 'SIF_771nm_uncert',
    'dcCorr': 'daily_correction_factor',
    'time': 'time'
    }

dict_tropomi = {
   #'variable_name': 'standard_name_from_netcdf',
    'lat': 'lat',
    'lon': 'lon',
    'sza': 'sza',
    'vza': 'vza',
    'phaseAngle': 'phase_angle',
    'sif': 'sif',
    'sif_sigma': 'sif_err',
    'dcCorr': 'daily_correction_factor',
    'time': 'TIME'
    }

# A more generic class to read in desired HDF5 files (doesn't have to be SIF, can be anything):
# Could add more pre-filters later but lat/lon boxes are as generic as it can be.
class L2:
	def __init__(self, path, dictionary=dict_oco2, latMin=-90, latMax=90,lonMin=-180,lonMax=180):
		files = glob.glob(path)
		# Check whether data has been initialized
		nini = True
		for file in files:
			print(file)

			h = h5py.File(file,'r')
			#h = Dataset("file, "r", format="NETCDF4")
			print('opening ', file)
			try:
				lat = h[dictionary['lat']][:]
				lon = h[dictionary['lon']][:]
				# find right indices:
				wo = (lat>=latMin)&(lat<=latMax)&(lon>=lonMin)&(lon<=lonMax)
				n = len(np.where(wo)[0])
				if nini and n>0:
					# For some reason, this doesn't work yet with my OCO-2 files, something is wrong with the files.
					#try:
					#	self.t_unit = h[dictionary['time']].attrs.get('unit')
					#except:
					#	self.t_unit = h[dictionary['time']].attrs.get('units')
					for k,v in dictionary.items():
						setattr(self, k, h[v][wo])
					nini = False
				elif n>0:
					for k,v in dictionary.items():
						# This could be improved, right now it might resize these variable all the time
						temp = np.hstack((getattr(self, k), h[v][wo]))
						setattr(self, k, temp)
				h.close()
			except:
				print(traceback.format_exc())
				print('Error opening file ', file)
				h.close()

# Compute Phase angle from SZA, SAA, VZA, VAA (via Philipp):
def compPhase(sza, saa, vza, vaa):
	p = 180./np.pi
	phase = np.zeros(len(sza))
	phase[vaa>saa]=-1
	phase[vaa<saa]=1
	# relative azimuth:
	raa = vaa-saa
	raa[raa<-180]+=360
	raa[raa>180]-=360
	raa = np.abs(raa)
	cos_theta = np.cos(vza/p)*np.cos(sza/p) + np.sin(vza/p)*np.sin(sza/p)*np.cos(raa/p);
	return np.arccos(cos_theta)*p*phase
