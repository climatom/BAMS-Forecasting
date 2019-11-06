# Coordinating Python script to run GFS download
# -*- coding: utf-8 -*-
import sys,datetime, os, numpy as np

### Parameters
# Out directory - can be captured as the directory in which this python
# script is stored: 
outdir=os.path.dirname(os.path.realpath(__file__))+"/"
#outdir=sys.argv[1]
# Get number of hour to extract data for. Kept as command-line parameter for now. 
nhrs=sys.argv[1]
###

# Start by purging output directory of GFS.nc files 
fs = [ii for ii in os.listdir(outdir) if "GFS" in ii and ".nc" in ii]
if len(fs)>0:
	for ii in fs: os.remove(outdir+ii)

# Current datetime
date=datetime.datetime.now()

# Convert date to nearest 6 hours just gone
const=date.year*1e6+date.month*1e4+date.day*1e2 # base date 
thresh=date.year*1e6+date.month*1e4+date.day*1e2+date.hour - 5.5 # cant be later than this!
# - the GFS forecst seems to be posted ~ 4-5 hours after initialization time; subtracting 
# 5.5 hours ensures we don't ask for a forecast that's not yet been posted!

# find previous times
poss=np.array([const+ii for ii in range(0,24,6)],dtype=np.int); idx=poss<=thresh 
sel=poss[idx][-1] # get the latest from the previous times

# Pull out the hour
sel_hour=("%02d"%sel)[-2:]

# Create datetime from sel (for conversion to isoformat)
sel="%.0f"%sel
pytime=datetime.datetime(year=np.int(sel[:4]),month=np.int(sel[4:6]),day=np.int(sel[6:8]),hour=np.int(sel[8:]))
isotime=pytime.isoformat()
sel=sel[:-2]+"\%2F"+sel[-2:]

# Call bash script to actually get the data. Note that, in turn, 
# Download.sh also calls Process.py to interpolate to locations. 
# All calls are to programs in outdir (which is the folder holding
# this script!)
fail=os.system("bash %sDownload.sh %s %s %s %s %s" %(outdir,outdir,sel,sel_hour,nhrs,isotime))
count=1
while fail != 0 and count <11:
	# Purge output directory of GFS.nc files 
	fs = [ii for ii in os.listdir(outdir) if "GFS" in ii and ".nc" in ii]
	if len(fs)>0:
		for ii in fs: os.remove(outdir+ii)

	fail=os.system("bash %sDownload.sh %s %s %s %s %s" %(outdir,outdir,sel,sel_hour,nhrs,isotime))
	count +=1 

print fail,count
