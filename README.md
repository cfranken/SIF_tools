# SIF_tools
some tools for accessing OCO-2 SIF lite data
Requirements: You will need the following python libraries:
import h5py
from netCDF4 import Dataset

Anaconda has those typically pre-installed.

The main script is now a little slow as I have time-consuming loops (it was faster before, had to change something as netCDF gave me issues with indexing). What you need is the following
1) A file with 3 columns, orbit start, orbit stop and day after 1.1.2014
e.g. if you want to bin ALL data within one lat/lon scheme write something like
0 1e99 100
which will use all orbits (between 0 and 1e99) and just give it a time of Day 100 in 2014. I have an example file in the directory as well

Run the code like this:
gridSIF_lite_time_SDOS.py --mode=1 /scf_archive/oco2/test_runs/Lite/B7101/LtSIF/oco2_LtSIF_150620_B7101_150903004045s.nc4

Look for options like this:
gridSIF_lite_time_SDOS.py --help

NOTE:
This is really a quick&dirty script, so don't expect it to be fool-proof. The good thing is that you can use tools like Panoply to look at the data later (http://www.giss.nasa.gov/tools/panoply/). If you improve the script, let me know and I can make changes available to everyone. 

in the dictionary (dict_l2), you can add variables to be gridded. 
