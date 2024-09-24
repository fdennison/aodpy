import datetime as dt
import pandas as pd
import numpy as np
import os as os
import math
import fileread_module as fr
import argparse

parser=argparse.ArgumentParser(description="Process raw environment files (.rnv)")
parser.add_argument('-t', '--test', dest='configtoml', action='store_const', 
                    const='../../tests/testconfig.toml', default='./config.toml',
                    help='run on sample data and check output')
args=parser.parse_args()

# langley config options
with open(args.configtoml, "rb") as f:
    envconf = tomllib.load(f)
print(envconf)    
site = envconf['site']
rootpath = envconf['rootpath']
startdate = envconf['startdate']
enddate = envconf['enddate']

#Logger clock error (+ve for slow) 
ce=dt.timedelta(hours=0,minutes=0,seconds=0)

configfile = rootpath+'config/'+site+'.cfn'
config = fr.Read_Site_Configuration(configfile, startdate)
stationlat = config.attrs['Latitude'] 
stationlon = config.attrs['Longitude'] 
model = config.CimelModel
inst = config.CimelNumber
card = config.Card
loggerversion = 'v' + str(config.Version).zfill(1)

def write_env_header(filename, Station, TZ, Instrument, Barometer, Card, BaroCal, CardforBarometer,
                     BaroCalDate, Anemometer, WCalDate, LoggerVersion, VCalDate):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w") as f:
        f.write('# Station              {}\n'.format(Station))
        f.write('# Time offset          {}\n'.format(TZ))
        f.write('# Cimel                 {}\n'.format(Instrument))
        f.write('# Barometer, Card at site      {0:d}  {1:d}\n'.format(Barometer,Card))
        f.write('# Barometer, Card used in cal  {0:d}  {1:d}   Cal Date {2}\n'.format(int(BaroCal),int(CardforBarometer),BaroCalDate))
        f.write('# Anemometer            {0}            Cal Date    {1}\n'.format(Anemometer,WCalDate))
        f.write('# LoggerVersion         {0}\n'.format(LoggerVersion))
        f.write('# Card                  {0}            Cal Date {1}\n'.format(Card,VCalDate))


caldir = rootpath+'/analog.cal_py/'

if card>0:
    vcal = pd.read_csv(caldir+'vcal.men', skiprows=3, header=None, names=['card','y','m','d','vcoef0','vcoef1'],  delimiter=r'\s+', usecols=range(0,6), index_col=0)
    vcal = vcal.loc[card]
    #vcoef = [vcal.vcoef0, vcal.vcoef1]
    
    tcal = pd.read_csv(caldir+'tcal.men', skiprows=3, header=None, names=['y','m','d','tcoef0','tcoef1'], delimiter=r'\s+')
    tcal = tcal.loc[card]
    #tcoef = [tcal.tcoef0, tcal.tcoef1]

    vcaldate = dt.date(int(vcal.y),int(vcal.m),int(vcal.d))
else:
    vcaldate = []
    
if config.Wind>0:
    wcal = pd.read_csv(caldir+'wcal.men', skiprows=3, header=None, names=['y','m','d','wcoef0','wcoef1'], delimiter=r'\s+')
    wcal = wcal.loc[config.Wind]
    #wcoef = [wcal.wcoef0, wcal.wcoef1]
    wcaldate = dt.date(int(wcal.y),int(wcal.m),int(wcal.d))
else:
    wcaldate = []

if config.Baro>0:
    pcal_list = pd.read_csv(caldir+'pcal.dat', skiprows=4, header=None,  delimiter=r'\s+',
                            names=['y','m','d','baro', 'card', 'pcoef0', 'pcoef1', 'ptcoef0', 'ptcoef1'])
    pcal_list.fillna(0, inplace=True)
    pcal = pcal_list.loc[(pcal_list.card==card) & (pcal_list.baro==config.Baro)]
    if len(pcal)==0:
        # print(f'No P cal that matches both card number [{card}] and baro number [{config.Baro}]')
        # print(pcal_list)
        # pcal_idx = int(input('choose cal by entering the index from the following list'))
        # pcal = pcal_list.loc[pcal_idx]
        pcal = pcal_list[pcal_list.baro==config.Baro].iloc[0]
    else:
        pcal = pcal.iloc[0]
        
    pcaldate = dt.date(int(pcal.y),int(pcal.m),int(pcal.d))
else:
    pcaldate = []

datelist = pd.date_range(startdate,enddate,freq='d')  

outrootpath = rootpath + 'PyOut/' 
inrootpath = rootpath + 'agsdat/' 

for dii in datelist:
    obsdate = pd.to_datetime(dii).date()
    
    subdir = site + '/#' + str(inst).zfill(2) + '/' + obsdate.strftime("%Y") +'/'+ obsdate.strftime("%m") +'/'+ obsdate.strftime("%d") + '/'
    fileroot = obsdate.strftime("%y%m%d")
    

    infile = inrootpath +  subdir + fileroot + '.rnv'
    if os.path.isfile(infile):
        env = pd.read_csv(infile, header=None, names=['Year', 'Month', 'Day', 'Hour', 'minute', 'Second','Pmean_in', 'Pmin_in', 'Pmax_in', 'Psd_in',
                                                      'Tmean_in', 'Tmin_in', 'Tmax_in', 'Tsd_in','WVmean_in', 'WVmin_in', 'WVmax_in', 'WVsd_in',
                                                      'WDmean_in', 'WDmin_in', 'WDmax_in', 'WDsd_in','VBatt_in'], delimiter=r'\s+')
        
        env['DateTime'] = [dt.datetime(y,m,d,h,mi,s) for y,m,d,h,mi,s in zip(env.Year, env.Month, env.Day, env.Hour, env.minute, env.Second)]
        
        env['Tmean'] = tcal.tcoef0 + tcal.tcoef1 * env.Tmean_in
        env['Tmin'] = env.Tmean - tcal.tcoef1 * env.Tmin_in
        env['Tmax'] = env.Tmean + tcal.tcoef1 * env.Tmax_in
        env['Tsd'] = tcal.tcoef1 * env.Tsd_in
        
        if config.Baro>0:
            env['Pmean'] = pcal.pcoef0 + pcal.pcoef1 * env.Pmean_in
            env['Pmin'] = env.Pmean - pcal.pcoef1 * env.Pmin_in
            env['Pmax'] = env.Pmean + pcal.pcoef1 * env.Pmax_in
            env['Psd'] = pcal.pcoef1 * env.Psd_in
        
            env['Pressure_UnTComp'] = env.Pmean
            env['DelPressure'] = pcal.ptcoef1 * (env.Tmean - pcal.ptcoef0)
            env['Pmean'] = env.Pmean - env.DelPressure
            env['Pmin'] = env.Pmin - env.DelPressure
            env['Pmax'] = env.Pmax - env.DelPressure
        else:
            env['POutmean'] = 0
            env['Pmin'] =0
            env['Pmax'] = 0
            env['Psd'] = 0
        
        if config.Wind>0:
            env['WindVmean'] = wcal.wcoef0 + wcal.wcoef1 * env.WVmean_in / 100
            env['WindVmin'] = WindV.mean - wcal.wcoef1 * env.WVmin_in / 100
            env['WindVmax'] = WindV.mean + wcal.wcoef1 * env.WVmax_in / 100
            env['WindVsd'] = wcal.wcoef1 * env.WVsd_in / 100
        else:
            env['WindVmean'] = -1.0
            env['WindVmin'] = -1.0
            env['WindVmax'] = -1.0
            env['WindVsd'] = -1.0
        
        if config.Wind>0:
            env['WindDmean'] = env.WDmean_in / 100
            # Boom points west not north, so rotate axes by -90 (+270) degrees.
            env['WindDmean'] = (env.WindDmean + 270.0) % 360.0
            env['WindDmin'] = env.WDmean - env.WDmin_in / 100
            env['WindDmax'] = env.WDmean + env.WDmax_in / 100
            env['WindDsd'] = env.WDsd_in / 100
        else:
            env['WindDmean'] = 0
            env['WindDmin'] = 0
            env['WindDmax'] = 0
            env['WindDsd'] = 0
        
        env['VBatt'] = vcal.vcoef0 + vcal.vcoef1 * env.VBatt_in
        
        # Add Time zone offset (TZ=-10h for AEST)
        env['DateTimeUTC'] = env.DateTime + config.TOffset
        
        # Add clock error.
        # Negative error: Fast clock
        # Positive error: Slow clock
        # NB This convention differs from that used before October 2000
        env['DateTimeUTC'] = env.DateTimeUTC + ce
    
        
        hpafile = outrootpath + subdir + fileroot + '.hpa'
        os.makedirs(os.path.dirname(hpafile), exist_ok=True)
        write_env_header(hpafile, site, config.TOffset , inst, config.Baro, card, pcal.baro, pcal.card,
                         pcaldate, config.Wind, wcaldate, loggerversion, vcaldate)
        with open(hpafile, 'a') as f:
            f.write('#YYYY MM dd hh mm ss P[raw]  Tmean   DelP P[comp]\n')   
            f.write(env.to_string(header=False, index=False, columns=['DateTimeUTC', 'Pressure_UnTComp', 'Tmean', 'DelPressure', 'Pmean'],
                                                             formatters={'Pressure_UnTComp':'{:6.1f}'.format, 'Tmean':'{:6.1f}'.format, 'DelPressure':'{:6.1f}'.format, 'Pmean':'{:6.1f}'.format}))
            f.write('\n')   
        
        envfile = outrootpath + subdir + fileroot + '.env'
        os.makedirs(os.path.dirname(envfile), exist_ok=True)
        write_env_header(envfile, site, config.TOffset , inst, config.Baro, card, pcal.baro, pcal.card,
                         pcaldate, config.Wind, wcaldate, loggerversion, vcaldate) 
        with open(envfile, 'a') as f:
            f.write('#YYYY MM dd hh mm ss            Pressure               Temperature               Wind speed (m/s)         Wind direction (E of N)   Battery voltage\n')
            f.write('#                     mean    min    max   sdev   mean    min    max   sdev   mean    min    max   sdev   mean    min    max   sdev\n') 
            f.write(env.to_string(header=False, index=False, columns=['DateTimeUTC', 'Pmean', 'Pmin', 'Pmax', 'Psd',
                                                                  'Tmean', 'Tmin', 'Tmax', 'Tsd',
                                                                  'WindVmean', 'WindVmin', 'WindVmax', 'WindVsd',
                                                                  'WindDmean', 'WindDmin', 'WindDmax', 'WindDsd',
                                                                  'VBatt'],
                                                           formatters={'Pmean':'{:6.1f}'.format, 'Pmin':'{:6.1f}'.format, 'Pmax':'{:6.1f}'.format, 'Psd':'{:6.2f}'.format,
                                                                  'Tmean':'{:6.1f}'.format, 'Tmin':'{:6.1f}'.format, 'Tmax':'{:6.1f}'.format, 'Tsd':'{:6.2f}'.format,
                                                                  'WindVmean':'{:6.1f}'.format, 'WindVmin':'{:6.1f}'.format, 'WindVmax':'{:6.1f}'.format, 'WindVsd':'{:6.2f}'.format,
                                                                  'WindDmean':'{:6.1f}'.format, 'WindDmin':'{:6.1f}'.format, 'WindDmax':'{:6.1f}'.format, 'WindDsd':'{:6.2f}'.format,
                                                                  'VBatt':'{:6.2f}'.format}))
            f.write('\n') 

    else:
        print(f'no {infile}')
   

