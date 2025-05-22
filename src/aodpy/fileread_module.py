import pandas as pd
import datetime as dt
import numpy as np
import atmos_module as atm

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
        #print(ConfigOut)
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
        p.sort_index(inplace=True)       # jb1 160105 has a dicontinuity, indicate problem?
    if len(p)==0:
        print(f'problem with pressure file: {presfile}')
    return p  
    
def read_old_pressure_file(presfile):
    prescolumnnames=['y','m','d','H','M','S','P_UnTCorrect','Tmean','DelP','Pressure']
    p = pd.read_csv(presfile, skiprows=9, header=None, delimiter=r'\s+', names=prescolumnnames, index_col=False)
    #p['DateTime'] = pd.to_datetime(p['Date'] + ' ' + p['Time'])
    p['DateTime'] = [dt.datetime(y,m,d,H,M,S) for y,m,d,H,M,S in zip(p['y'],p['m'],p['d'],p['H'],p['M'],p['S'])]
    p = p.set_index('DateTime')
    if np.any(np.diff(p.index.to_list())<dt.timedelta(hours=0)): 
        print('time disordered in pressure file')
        p.sort_index(inplace=True)       # jb1 160105 has a dicontinuity, indicate problem?
    if len(p)==0:
        print('problem with pressure file')
    return p  
     
minsignal=50
maxsignal=65000
mintemp=-10.0
maxtemp=60.0 

def checklowsignal(allchan):
    return int(any(ch < minsignal for ch in allchan))*10
def checkhighsignal(allchan):
    return int(any(ch > maxsignal for ch in allchan))*11 

def adjust1020(signal,temp):
    
    gaincoefficient = 0.0025

    if mintemp < temp < maxtemp:
        tempdep1020 = 1.0 + gaincoefficient * (temp - 25.0)
        #TempLastValid = temp
    else:
        tempdep1020 = 1.0
        
    adjsignal = signal/tempdep1020    
    return adjsignal

def adjust1020_trip(signal,temp):   # todo: combine this into a single adjust function
    
    gaincoefficient = 0.0025

    if mintemp < temp < maxtemp:
        tempdep1020 = 1.0 + gaincoefficient * (temp - 25.0)
        #TempLastValid = temp
    else:
        tempdep1020 = temp
        
    adjsignal = list(np.array(signal)/tempdep1020)    
    return adjsignal

def read_triple_sun_record(filename, Model):
    if Model in [1,3]:
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
        if len(s.columns) != (nchan*3+3):
            print(f'{filename} has unexpected number of channels')
            s = pd.DataFrame([])
    except:
        s = pd.DataFrame([])
        print(f'problem with {filename}')
        
    if not s.empty:   
        s.rename(columns={0:'Date',1:'Time',(3*nchan+2):'Temperature'}, inplace=True)
        s['Date'] = pd.to_datetime(s['Date'], format='%d:%m:%y')
        s['DateTime'] = pd.to_datetime(s['Date'], format='%d:%m:%y')+pd.to_timedelta(s['Time'])
        s = s.set_index('DateTime')
    
        if not s.index.is_monotonic_increasing:
            print(f'{filename} has disordered dates')
            s.sort_index(inplace=True)
            
        for ii,ch in enumerate(order):
            name='ch{}'.format(int(atm.wavelength[Model-1][ch-1]))
            s[name] = [[x,y,z] for x,y,z in zip(s[ii+2], s[ii+2+nchan], s[ii+2+2*nchan])] 
    
        s['allchan'] = s.loc[:, 2:(3*nchan+1)].to_numpy().tolist() 
        s['Low'] = s['allchan'].apply(checklowsignal)
        s['High'] = s['allchan'].apply(checkhighsignal)
        # s['wvLow'] = 0; s.loc[s['ch'+str(IWV)] < minsignal, 'wvLow'] == 101
        # s['wvHigh'] = 0; s.loc[s['ch'+str(IWV)] > maxsignal, 'wvHigh'] == 111 
        s['tempLow'] = 0; s.loc[s['Temperature'] < mintemp, 'tempLow'] == 12
        s['tempHigh'] = 0; s.loc[s['Temperature'] > maxtemp, 'tempHigh'] == 12
        s['Valid'] = s[['Low','High','tempLow','tempHigh']].sum(axis=1) #,'wvLow','wvHigh'
        
        s['RawSignal1020'] = s['ch1020']
        s['ch1020'] = [adjust1020_trip(x,y) for x,y in zip (s['ch1020'],s['Temperature'])]        
        s['allchan'] = s[range(2,(3*nchan+2))].to_numpy().tolist()
    
    return s


def read_single_sun_record(filename, Model):
    if Model in [1,3]:
        order = [4,3,2,1,7,8,6,5]
    elif Model==2:
        order = [4,5,2,1,6,3,8,7]
    elif Model==4:
        order = [4,3,2,1,6,8,7,5]
    elif Model in [5,6]:
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
        #s['Date'] = s[0].apply(lambda x: dt.datetime.strptime(x,'%d:%m:%y').date())
        #s['Time'] = s[1].apply(lambda x: dt.datetime.strptime(x,'%H:%M:%S').time())
        
        #s['DateTime'] = [dt.datetime.combine(dt.datetime.strptime(d,'%d:%m:%y').date(),dt.datetime.strptime(t,'%H:%M:%S').time()) for d,t in zip(s['Date'],s['Time'])]
        s['DateTime'] = pd.to_datetime(s['Date'], format='%d:%m:%y')+pd.to_timedelta(s['Time'])
        s = s.set_index('DateTime')    
        
        for ii,ch in enumerate(order):
            s.rename(columns={(ii+2):('ch{}'.format(int(atm.wavelength[Model-1][ch-1])))}, inplace=True)
            
        s['RawSignal1020'] = s['ch1020']
        s['ch1020'] = [adjust1020(x,y) for x,y in zip (s['ch1020'],s['Temperature'])]
        s['allchan'] = s.loc[:, s.columns.str.startswith('ch')].to_numpy().tolist()  
        s['Low'] = s['allchan'].apply(checklowsignal)
        s['High'] = s['allchan'].apply(checkhighsignal)
        s['wvLow'] = 0; s.loc[s['ch936'] < minsignal, 'wvLow'] == 101
        s['wvHigh'] = 0; s.loc[s['ch936'] > maxsignal, 'wvHigh'] == 111 
        s['tempLow'] = 0; s.loc[s['Temperature'] < mintemp, 'tempLow'] == 12
        s['tempHigh'] = 0; s.loc[s['Temperature'] > maxtemp, 'tempHigh'] == 12
        s['Valid'] = s[['Low','High','wvLow','wvHigh','tempLow','tempHigh']].sum(axis=1)
    
    return s



def read_bom_ozone2011(ozonefile):
    o3 = pd.read_csv(ozonefile, header=None, delimiter=r'\s+')
    o3.rename(columns={0:'id_date',1:'ozone'}, inplace=True)
    o3['date'] = [dateint2dateclass(str(i)[3:11]) for i in o3['id_date']]
    

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
