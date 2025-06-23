import os
import h5py
import numpy as np
import toml
import glob
import time
import pandas as pd
import datetime as dt
import fileread_module as fr
import ozone_module as om
import argparse

omidir = '/g/data/p66/fd0474/data/OMI/'

# Earthdata user credentials
username = ''
password = ''

stations = [
{'name':'jb1', 'lat':-12.6607,  'lon':132.8931,  'id':853},
{'name':'la3', 'lat':-16.1081 , 'lon':128.7485,  'id':952},
{'name':'bvl', 'lat':-25.89893, 'lon':139.34596, 'id':452},
{'name':'asp1','lat':-38.02694, 'lon':145.1000,  'id':000},
{'name':'bml', 'lat':-35.2713,  'lon':149.1111,  'id':252},
{'name':'fgp', 'lat':-31.08633, 'lon':141.70082, 'id':299},
{'name':'lea', 'lat':-22.2407,  'lon':114.0967,  'id':000},
{'name':'llf', 'lat':-31.255,   'lon':121.7050,  'id':953},
{'name':'luc', 'lat':-18.5198,  'lon':146.3681,  'id':451},
{'name':'pin', 'lat':-30.58361, 'lon':115.155,   'id':000},
{'name':'tum', 'lat':-35.7483,  'lon':147.9499,  'id':000}
]

parser=argparse.ArgumentParser(description="update {site}.o3 file using OMI ozone data")
parser.add_argument('--site' type=str, help='station code e.g. jb1')
parser.add_argument('--date',type=lambda s: dt.datetime.strptime(s, '%Y-%m-%d').date(),
                     help='Specify a date in YYYY-MM-DD format')

with open("/g/data/p66/fd0474/notebooks/AOD/config.toml", "r") as f:
    conf = toml.load(f)
if agrs.site is not None:
    site = args.site
else:    
    site = aod2pconf['site']
rootpath = aod2pconf['rootpath']
if args.date is not None:
    startdate = args.date
    enddate = args.date
else:    
    startdate = aod2pconf['startdate']
    enddate = aod2pconf['enddate']
rootpath = conf['rootpath']
ozonedir = conf['ozonedir']
stationdetail = [s for s in stations if s['name']==site][0]

# check for ozone file
filename = f'{rootpath}{ozonedir}/{site}.o3'
if os.path.isfile(filename):    
    print(f'{site}.o3 exists')
    o3_old = fr.read_ozone(filename)
else:
    print(f'{filename} not present')

# identify missing dates
datelist = pd.date_range(startdate,enddate,freq='d') 
missing = datelist.difference(o3_old.index)
if len(missing)>0:
    print(f'{site}.o3 missing dates:')
    [print(x) for x in missing]
    
    # download OMI data if necessary
    missing_omi = om.check_missing_files(missing, omidir)
    if len(missing_omi)>0:
        print('OMI files missing for dates:')
        [print(x) for x in missing_omi]
        dl = input('download the files (y/n) ?')
        if dl=='y':
            om.get_missing_files(missing_omi, omidir, username, password)
    
    # process OMI data for missing days
    o3data= om.read_omi(omidir,missing,[stationdetail])[0]
    id_date = [f'{stationdetail["id"]}{d.strftime("%Y%m%d")}' for d in missing]
    o3_add = pd.DataFrame({'id_date':id_date, 'ozone': o3data, 'date': [x.date() for x in missing]}).set_index('date') 
    if o3_add['ozone'].isnull().sum()>0:
        print('no OMI overpass in the region for the following dates:')
        [print(x) for x in o3_add[o3_add['ozone'].isnull()].index.to_numpy()]
        o3_add.dropna(inplace=True)
    o3_add['ozone'] = o3_add['ozone'].astype(int)
    o3 = pd.concat([o3_add,o3_old]).sort_index()
    
    # update ozone file
    o3.to_csv(filename, sep=' ', index=False, header=False)
    
    print('ozone file updated')
else:
    print('ozone file up to date')

