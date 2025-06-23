import h5py
import numpy as np
import glob
import pandas as pd
import datetime as dt
import os
import urllib.request
import re
import requests
from pathlib import Path

base_url = 'https://acdisc.gsfc.nasa.gov/data/Aura_OMI_Level3/OMDOAO3e.003/'

# from: https://urs.earthdata.nasa.gov/documentation/for_users/data_access/python

# overriding requests.Session.rebuild_auth to mantain headers when redirected
class SessionWithHeaderRedirection(requests.Session):
    AUTH_HOST = 'urs.earthdata.nasa.gov'
    def __init__(self, username, password):
        super().__init__()
        self.auth = (username, password)

   # Overrides from the library to keep headers when redirected to or from
   # the NASA auth host.
    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url

        if 'Authorization' in headers:
            original_parsed = requests.utils.urlparse(response.request.url)
            redirect_parsed = requests.utils.urlparse(url)

            if (original_parsed.hostname != redirect_parsed.hostname) and \
                redirect_parsed.hostname != self.AUTH_HOST and \
                original_parsed.hostname != self.AUTH_HOST:
                del headers['Authorization']
        return


def get_file(url,downloaddir,username, password):
    # create session with the user credentials that will be used to authenticate access to the data
    session = SessionWithHeaderRedirection(username, password)
    try:    
        # submit the request using the session 
        response = session.get(url, stream=True)    
        #print(response.status_code)

        # raise an exception in case of http errors    
        response.raise_for_status()  
 
        # extract the filename from the url to be used when saving the file
        filename = url[url.rfind('/')+1:]  
        
        # save the file   
        #t1 = time.time() 
        with open(downloaddir+filename, 'wb') as fd:    
            for chunk in response.iter_content(chunk_size=1024*1024):    
                fd.write(chunk)
        #print(tses-t0,time.time()-t1)
    except requests.exceptions.HTTPError as e:    
        # handle any errors here
        print(e)
        
def check_missing_files(date_list, base_directory, file_prefix="OMI-Aura_L3-OMDOAO3e", file_extension=".he5"):
    
    base_path = Path(base_directory)
    missing_dates = []
    
    for date in date_list:        
        # Format date as YYYY-MM-DD for the filename pattern
        date_str = date.strftime("%Ym%m%d")
        year = date.strftime("%Y")
        
        # Construct the directory path (year subdirectory)
        year_dir = base_path / year
        
        # Check if year directory exists
        if not year_dir.exists():
            missing_dates.append(date)
            continue
        
        # Search for files matching the pattern
        # The pattern includes wildcard for the timestamp part
        pattern = f"{file_prefix}_{date_str}_v*{file_extension}"
        matching_files = list(year_dir.glob(pattern))
        
        # If no matching files found, add to missing dates
        if not matching_files:
            missing_dates.append(date)
    
    return pd.DatetimeIndex(missing_dates)  

def get_missing_files(missing_dates, ozonedir, username, password):
    yprev=0
    for d in missing_dates:
        if d.year != yprev:
            response = urllib.request.urlopen(base_url+str(d.year)+'/')
            html = response.read().decode('utf-8')
            yprev=d.year
            # create session with the user credentials that will be used to authenticate access to the data
            session = SessionWithHeaderRedirection(username, password)
            response = session.get(base_url, stream=True)
        
        filename  = f'OMI-Aura_L3-OMDOAO3e_{d.year}m{d.month:02}{d.day:02}_v003*'
        
        if not glob.glob(ozonedir+str(d.year)+'/'+filename):
            pattern = re.compile(rf'{filename}.*?\.he5')
            fn = pattern.findall(html)
            if not fn:
                print(f'no {pattern}')
                continue
            else:
                fn = fn[0]
                
            get_file(base_url+str(d.year)+'/'+fn, ozonedir+str(d.year)+'/', username, password)    
            print(f'downloaded {fn}')
            
def near(x,xarr,d):
    xidx = np.logical_and(xarr>(x-d/2), xarr<(x+d/2))
    return xidx

def read_omi(path, dates,stations):
    omi_o3 = []
    for d in dates:
        fn = glob.glob(f'{path}{d.year}/OMI-Aura_L3-OMDOAO3e_{d.year}m{str(d.month).zfill(2)}{str(d.day).zfill(2)}*')
        if len(fn)==1:
            fn = fn[0]
        else:
            print(f'multiple files for {d}')    
        hf = h5py.File(fn)
        offset = hf['HDFEOS']['GRIDS']['ColumnAmountO3']['Data Fields']['ColumnAmountO3'].attrs['Offset']
        sf = hf['HDFEOS']['GRIDS']['ColumnAmountO3']['Data Fields']['ColumnAmountO3'].attrs['ScaleFactor']
        o3ii = hf['HDFEOS']['GRIDS']['ColumnAmountO3']['Data Fields']['ColumnAmountO3'][:] * sf + offset
        o3ii[np.where(o3ii==hf['HDFEOS']['GRIDS']['ColumnAmountO3']['Data Fields']['ColumnAmountO3'].attrs['_FillValue'])] = np.nan
        omi_o3.append(o3ii)
    omi_o3 = np.stack(omi_o3) 
    
    a = hf['HDFEOS']['GRIDS']['ColumnAmountO3'].attrs
    gspan = [float(x) for x in a['GridSpan'].decode('UTF-8')[1:-1].split(',')]
    gspace = [float(x) for x in a['GridSpacing'].decode('UTF-8')[1:-1].split(',')]
    lon = np.arange(gspan[0],gspan[1],gspace[0])
    lat = np.arange(gspan[2],gspan[3],gspace[1])

    station_o3 = []
    for station in stations:
        x = omi_o3[:,near(station['lat'],lat,5),:][:,:,near(station['lon'],lon,10)]
        x = x.reshape(x.shape[0],-1)
        if np.all(np.isnan(x)):
            station_o3.append(np.nan)
        else:    
            station_o3.append(np.nanmean(x, axis=1))      

    return station_o3
