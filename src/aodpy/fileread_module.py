import pandas as pd
import datetime as dt
import numpy as np
import atmos_module as atm
import os as os
import constants_module as c

def adjust1020(signal,temp):  
    
    gaincoefficient = 0.0025

    if c.mintemp < temp < c.maxtemp:
        tempdep1020 = 1.0 + gaincoefficient * (temp - 25.0)
    else:
        tempdep1020 = 1.0
        
    adjsignal = signal/tempdep1020    
    return adjsignal

def read_all_triple_sun_records(rootpath,site,startdate,enddate,inst,model):
    # load photometer data
    if model == 10:
        ext = '.NSU'
    else:    
        ext = '.sun'    
    numchannels = len(c.wavelength[model-1])
    
    datelist = pd.date_range(startdate-dt.timedelta(days=1),enddate+dt.timedelta(days=1),freq='d') # look in files 1 day either side of desired timeframe because of issues with aligning output files with local days
    sunfiles = [f'{rootpath}agsdat/{site}/#{str(inst).zfill(2)}/{str(di.year)}/{str(di.month).zfill(2)}/{str(di.day).zfill(2)}/{di.strftime("%y%m%d")}{ext}' for di in datelist]
    s_all = []
    if model==10:
        dateformat='%d/%m/%Y'
    else:
        dateformat='%d:%m:%y'
    for sunfile in sunfiles:
        if not os.path.isfile(sunfile):   
            print(f"{sunfile.rsplit('/', 1)[-1]} not present")
        else:
            tempdf = pd.read_csv(sunfile, header=None).drop_duplicates()
            if len(tempdf.columns) != (numchannels*3+3):
                print(f'{filename} has unexpected number of channels,  has:{len(s.columns)} expect:{numchannels*3+3}')
            else:
                tempdf[0] = pd.to_datetime(tempdf[0], format=dateformat)
                tempdf['DateTime'] = pd.to_datetime(tempdf[0], format=dateformat)+pd.to_timedelta(tempdf[1])
                tempdf = tempdf.set_index('DateTime')
                s_all.append(tempdf)
                
    s_all = pd.concat(s_all).drop_duplicates()            
    if not s_all.index.is_monotonic_increasing:
        s_all.sort_index(inplace=True) 
        
    s_all.rename(columns={0:'Date',1:'Time',(3*numchannels+2):'Temperature'}, inplace=True)
    for j in range(numchannels):
        for k in range(3):
            s_all.rename(columns={(k*numchannels+j+2):f'ch{int(c.wavelength[model-1][j])}_{k}'}, inplace=True)
    f1 = s_all
    t_low_end = s_all['Temperature'] > c.mintemp 
    t_high_end = s_all['Temperature'] < c.maxtemp
    low_end = (s_all.loc[:,s_all.columns.str.contains('ch*')]>c.minsignal).all(axis=1)
    high_end = (s_all.loc[:,s_all.columns.str.contains('ch*')]>c.maxsignal).all(axis=1)
    
    # adjust 1020nm channel for temperature dependence
    for k in range(3):
        s_all[f'ch1020_{k}'] = [adjust1020(x,y) for x,y in zip(s_all[f'ch1020_{k}'],s_all['Temperature'])]
        
    s_all
    s = s_all[low_end & t_low_end & t_high_end] # #   & high_end   New inst. has much larger values - so high end filter needs to be adjusted

    return s

def read_adj_file(filename,var):
    df = pd.read_csv(filename, skiprows=1, header=None, delimiter=r'\s+', 
                 names=['year','month','day',var])
    df['datetime'] = pd.to_datetime(df[['year','month','day']])
    df = df.set_index('datetime')
    return df
    
def read_site_configuration(filename, ObsDate):

    Config = pd.read_csv(filename, skiprows=8, header=None, delimiter=r'\s+', usecols=range(0,14))
    Config.rename(columns={0:'StartDate', 1:'EndDate', 2:'TOffSetHour', 3:'TOffSetMinute', 4:'TOffSetSecond',
                           5:'CimelNumber', 6:'CimelModel', 7:'Version', 8:'Card', 9:'Baro', 10:'Wind', 11:'Neph',
                           12:'CimCalFile',13:'SkyCalDate'}, inplace=True)

    Config['TOffset'] = [dt.timedelta(hours=h,minutes=m,seconds=s) for h,m,s in zip(Config.TOffSetHour,Config.TOffSetMinute,Config.TOffSetSecond)]
    Config['StartDate'] = Config['StartDate'].apply(dateint2dateclass)
    Config['EndDate'] = Config['EndDate'].apply(dateint2dateclass)
    #Config['SkyCalDate'] = Config['SkyCalDate'].apply(dateint2dateclass)  # input file contains 'dummy'

    with open(filename,'r') as f:
        next(f)
        next(f)
        Name = f.readline().rstrip()
        BoMNumber = int(f.readline().rstrip())
        Id, lat, lon= f.readline().rstrip().split()
        Height = float(f.readline().rstrip().split()[0])
        TimeZone = float(f.readline().rstrip().split()[0])
    Config.attrs = {'Name':Name, 'BoMNumber':BoMNumber, 'Id':Id, 'Latitude':float(lat), 'Longitude':float(lon), 'Height':Height, 'TimeZone':TimeZone}
    pre_obs = [i for i,x in enumerate(Config['StartDate']) if (x<=ObsDate)]  #
    if not pre_obs : #empty
        print('could not find appropriate date')
    else:
        ConfigOut = Config.iloc[pre_obs[-1]]
        return ConfigOut

def dateint2dateclass(d):
    d=str(d)
    d=dt.date(int(d[0:4]),int(d[4:6]),int(d[6:8]))
    return d


def read_pressure_file(verb, presfile):
    prescolumnnames=['Date','Time','P_UnTCorrect','Tmean','DelP','Pressure']
    p = pd.read_csv(presfile, skiprows=9, header=None, delimiter=r'\s+', names=prescolumnnames, index_col=False)
    p['DateTime'] = pd.to_datetime(p['Date'] + ' ' + p['Time'])
    p = p.set_index('DateTime')
    if np.any(np.diff(p.index.to_list())<dt.timedelta(hours=0)):
        if verb: print(f'time disordered in pressure file: {presfile}')
        p.sort_index(inplace=True)       
    if len(p)==0:
        print(f'problem with pressure file: {presfile}')
    return p  
    
def read_old_pressure_file(presfile):
    prescolumnnames=['y','m','d','H','M','S','P_UnTCorrect','Tmean','DelP','Pressure']
    p = pd.read_csv(presfile, skiprows=9, header=None, delimiter=r'\s+', names=prescolumnnames, index_col=False)
    p['DateTime'] = [dt.datetime(y,m,d,H,M,S) for y,m,d,H,M,S in zip(p['y'],p['m'],p['d'],p['H'],p['M'],p['S'])]
    p = p.set_index('DateTime')
    if np.any(np.diff(p.index.to_list())<dt.timedelta(hours=0)): 
        print('time disordered in pressure file')
        p.sort_index(inplace=True)   
    if len(p)==0:
        print('problem with pressure file')
    return p  
def read_ozone(ozonefile):
    o3 = pd.read_csv(ozonefile, header=None, delimiter=r'\s+')
    o3.rename(columns={0:'id_date',1:'ozone'}, inplace=True)
    o3['date'] = [dateint2dateclass(str(i)[3:11]) for i in o3['id_date']]
    o3 = o3.set_index('date') 
    return o3

def read_cal(calfile):
    columnnames=['wvcal', 'order', 'lnV0coef0', 'lnV0coef1', 'erms']
    cal  = pd.read_csv(calfile, skiprows=10, delimiter=r'\s+', header=None, names=columnnames, usecols=range(5))
    with open(calfile,'r') as f:
        next(f)        
        d1 = [int(x) for x in f.readline().rstrip().split()[1:4]]
        d1 = dt.date(d1[0],d1[1],d1[2])
        d2 = [int(x) for x in f.readline().rstrip().replace('.',' ').split()[1:4]]
        d2 = dt.date(d2[0],d2[1],d2[2])
        istr = int(f.readline().rstrip().split()[1])
        model = int(f.readline().rstrip().split()[1])
        numlangleychannels = int(f.readline().rstrip().split()[1])
        numlangleys = int(f.readline().rstrip().split()[1])
        numgencycles = int(f.readline().rstrip().split()[1])
        d3 = [int(x) for x in f.readline().rstrip().split()[1:4]]
        cal_epoch = dt.date(d3[0],d3[1],d3[2])
    cal.attrs = {'epoch':cal_epoch, 'numlangleys':numlangleys, 'numlangleychannels':numlangleychannels,
                 'numgencycles':numgencycles} 
    return cal

def read_black_record(filename, Model):
    if Model in [1,3]:
        order = [4,3,2,1,7,8,6,5]
    elif Model==2:
        order = [4,3,2,1,5]
    elif Model==4:
        order = [4,3,2,1,6,8,7,5]
    colnames = ['date','time'] + [f'ch{i}' for i in order]
    blkrec = pd.read_csv(filename, header=None, names=colnames, usecols=range(len(colnames)))
    if Model in [1,3,4]:
        blkrec['ch9']=0
        blkrec['ch10']=0
    elif Model==2:
        blkrec['ch6']=-1

    blk = blkrec.iloc[:,2:].mean().to_list()
    return blk      


def checklowsignal(allchan):
    return int(any(ch < c.minsignal for ch in allchan))*10
def checkhighsignal(allchan):
    return int(any(ch > c.maxsignal for ch in allchan))*11 

def adjust1020_trip_old(signal,temp):       
    gaincoefficient = 0.0025

    if c.mintemp < temp < c.maxtemp:
        tempdep1020 = 1.0 + gaincoefficient * (temp - 25.0)
    else:
        tempdep1020 = 1.0
        
    adjsignal = list(np.array(signal)/tempdep1020)    
    return adjsignal

def read_triple_sun_record(filename, Model):
    if Model == 1:
        order = [4,3,2,1,5]
    elif Model == 3:
        order = [4,3,2,1,7,8,6,5]
    elif Model==2:
        order = [4,5,2,1,6,3,8,7]
    elif Model==4:
        order = [4,3,2,1,6,8,7,5]
    elif Model in [5,6,9,10]:
        order = [4,9,3,2,1,7,8,10,6,5]
    elif Model==7:
        order = [4,3,2,1,5,8,6,7]
    elif Model==8:
        order = [4,3,2,1,5,7,9,6,8]    
    nchan = len(order)

    if Model==10:
        dateformat='%d/%m/%Y'
    else:
        dateformat='%d:%m:%y'
        
    try:
        s = pd.read_csv(filename, header=None).drop_duplicates()
        if len(s.columns) != (nchan*3+3):
            print(f'{filename} has unexpected number of channels')
            print(f'has:{len(s.columns)} expect:{nchan*3+3}')
            s = pd.DataFrame([])
    except:
        s = pd.DataFrame([])
        print(f'problem with {filename}')
        
    if not s.empty:   
        s.rename(columns={0:'Date',1:'Time',(3*nchan+2):'Temperature'}, inplace=True)
        s['Date'] = pd.to_datetime(s['Date'], format=dateformat)
        s['DateTime'] = pd.to_datetime(s['Date'], format=dateformat)+pd.to_timedelta(s['Time'])
        s = s.set_index('DateTime')
    
        if not s.index.is_monotonic_increasing:
            #print(f'{filename} has disordered dates')
            s.sort_index(inplace=True)
            
        for ii,ch in enumerate(order):
            name='ch{}'.format(int(atm.wavelength[Model-1][ch-1]))
            s[name] = [[x,y,z] for x,y,z in zip(s[ii+2], s[ii+2+nchan], s[ii+2+2*nchan])] 
    
        s['allchan'] = s.loc[:, 2:(3*nchan+1)].to_numpy().tolist() 
        s['Low'] = s['allchan'].apply(checklowsignal)
        s['High'] = s['allchan'].apply(checkhighsignal)
        # s['wvLow'] = 0; s.loc[s['ch'+str(IWV)] < c.minsignal, 'wvLow'] == 101
        # s['wvHigh'] = 0; s.loc[s['ch'+str(IWV)] > c.maxsignal, 'wvHigh'] == 111 
        s['tempLow'] = 0; s.loc[s['Temperature'] < c.mintemp, 'tempLow'] == 12
        s['tempHigh'] = 0; s.loc[s['Temperature'] > c.maxtemp, 'tempHigh'] == 12
        s['Valid'] = s[['Low','High','tempLow','tempHigh']].sum(axis=1) #,'wvLow','wvHigh'
        
        s['RawSignal1020'] = s['ch1020']
        s['ch1020'] = [adjust1020_trip_old(x,y) for x,y in zip (s['ch1020'],s['Temperature'])]        
        s['allchan'] = s[range(2,(3*nchan+2))].to_numpy().tolist()
    
    return s


def read_single_sun_record(filename, Model):
    if Model == 1:
        order = [4,3,2,1,5]
    elif Model == 3:
        order = [4,3,2,1,7,8,6,5]
    elif Model==2:
        order = [4,5,2,1,6,3,8,7]
    elif Model==4:
        order = [4,3,2,1,6,8,7,5]
    elif Model in [5,6,9,10]:
        order = [4,9,3,2,1,7,8,10,6,5]
    elif Model==7:
        order = [4,3,2,1,5,8,6,7]
    elif Model==8:
        order = [4,3,2,1,5,7,9,6,8]
    nchan = len(order)

    try:
        s = pd.read_csv(filename, header=None).drop_duplicates()
        if len(s.columns) != (nchan+3):
            print(f'{filename} has unexpected number of channels')
            s = pd.DataFrame([])
    except:
        s = pd.DataFrame([])
        print(f'problem with {filename}')
        
    if not s.empty:    
        s.rename(columns={0:'Date',1:'Time',(nchan+2):'Temperature'}, inplace=True)
        s['DateTime'] = pd.to_datetime(s['Date'], format='%d:%m:%y')+pd.to_timedelta(s['Time'])
        s = s.set_index('DateTime')    
        
        for ii,ch in enumerate(order):
            s.rename(columns={(ii+2):('ch{}'.format(int(atm.wavelength[Model-1][ch-1])))}, inplace=True)
            
        s['RawSignal1020'] = s['ch1020']
        s['ch1020'] = [adjust1020_old(x,y) for x,y in zip (s['ch1020'],s['Temperature'])]
        s['allchan'] = s.loc[:, s.columns.str.startswith('ch')].to_numpy().tolist()  
        s['Low'] = s['allchan'].apply(checklowsignal)
        s['High'] = s['allchan'].apply(checkhighsignal)
        s['wvLow'] = 0; s.loc[s['ch936'] < c.minsignal, 'wvLow'] == 101
        s['wvHigh'] = 0; s.loc[s['ch936'] > c.maxsignal, 'wvHigh'] == 111 
        s['tempLow'] = 0; s.loc[s['Temperature'] < c.mintemp, 'tempLow'] == 12
        s['tempHigh'] = 0; s.loc[s['Temperature'] > c.maxtemp, 'tempHigh'] == 12
        s['Valid'] = s[['Low','High','wvLow','wvHigh','tempLow','tempHigh']].sum(axis=1)
    
    return s

