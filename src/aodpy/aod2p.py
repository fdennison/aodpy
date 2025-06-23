import datetime as dt
import pandas as pd
import numpy as np
import os as os
import math
import fileread_module as fr
import solpos_module as sol
import calfit_module as fit
import atmos_module as atm
import constants_module as c
import tomllib
import argparse

parser=argparse.ArgumentParser(description="Calculate AOD from sun photometer measuments (.sun/.NSU files)")
parser.add_argument('-t', '--test', dest='configtoml', action='store_const', 
                    const='testconfig.toml', default='config.toml',
                    help='run on sample data and check output')
parser.add_argument('-v', '--verbose',action='store_true')
parser.add_argument('--site', type=str, help='station code e.g. jb1')
parser.add_argument('--date',type=lambda s: dt.datetime.strptime(s, '%Y-%m-%d').date(),
                     help='Specify a date in YYYY-MM-DD format')
parser.add_argument('--calfile', type=str, help='full name of cal file')
args=parser.parse_args()

# aod2p config options
with open(args.configtoml, "rb") as f:
    aod2pconf = tomllib.load(f)
if args.site is not None:
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
ozonedir = aod2pconf['ozonedir']
if args.calfile is not None:
    calfile = args.calfile
else:    
    calfile =  aod2pconf['aod2p']['calfile']
clockfix = aod2pconf['aod2p']['clockfix']
cirrusopt = aod2pconf['aod2p']['cirrusopt']
CVmax = aod2pconf['aod2p']['CVmax']
solzenmax = aod2pconf['aod2p']['solzenmax']
minobs = aod2pconf['aod2p']['minobs']
relUaod440thres = aod2pconf['aod2p']['relUaod440thres']
Uangstrom48thres = aod2pconf['aod2p']['Uangstrom48thres']
print(f'Run aod2p for {site} from {startdate} to {enddate}')
verb = args.verbose
if verb: print(aod2pconf) 


# config file
configfile = rootpath+'config/'+site+'.cfn'
config = fr.read_site_configuration(configfile, startdate)
stationlat = config.attrs['Latitude'] 
stationlon = config.attrs['Longitude'] 
model = config.CimelModel
inst = config.CimelNumber
numchannels = len(c.wavelength[model-1])

#cal file
calperiod_filename = str(inst) + str(startdate.year % 100).zfill(2) + str(startdate.month).zfill(2) + str(enddate.month).zfill(2)
if len(calfile)==3:  # if calfile is just the extentsion then construct file name from start/end dates
    calfile = calperiod_filename+'.'+calfile
cal = fr.read_cal(rootpath+'suncals/'+calfile[0:-10].zfill(2)+'/'+calfile)
usignal=0.003
if cal.attrs['numlangleys']>0:
    cal['stdUaod1am'] = np.sqrt(usignal**2 + cal.erms**2/cal.attrs['numlangleys'])
else:
    cal['stdUaod1am'] = usignal

lnV0coef0 = [0]*numchannels
lnV0coef1 = [0]*numchannels
stdUaod1am = [0]*numchannels
for j, wl in enumerate(c.wavelength[model-1]):
    lnV0coef0[j] = cal[cal['wvcal']==wl]['lnV0coef0'].values[0]
    lnV0coef1[j] = cal[cal['wvcal']==wl]['lnV0coef1'].values[0]   
    stdUaod1am[j] = cal[cal['wvcal']==wl]['stdUaod1am'].values[0]
    
# ozone file
ozonefile = rootpath+ozonedir+site+'.o3'
if os.path.isfile(ozonefile):
    o3 = fr.read_ozone(ozonefile)
    ozoneopt=True
else:
    ozoneopt=False
    print('no ozone file, using default ozone: 250 DU')

if clockfix:  clk = fr.read_adj_file(rootpath + 'suncals/' + str(inst).zfill(2) +'/' + calperiod_filename + '.clk', 'timecorr')
    
i440 = c.wl_index(440,c.wavelength[model-1]) 
i670 = c.wl_index(670,c.wavelength[model-1])    
i870 = c.wl_index(870,c.wavelength[model-1]) 
i1020 = c.wl_index(1020,c.wavelength[model-1]) 
iWV = c.wl_index(936,c.wavelength[model-1]) 
langleych = np.delete(np.arange(numchannels),iWV)

s = fr.read_all_triple_sun_records(rootpath,site,startdate,enddate,inst,model)

datelist = pd.date_range(startdate,enddate,freq='d')  
for dii in datelist:
    obsdate = pd.to_datetime(dii).date()
    if verb: print('----'+obsdate.strftime("%Y/%m/%d")+'----')
    filedir = site+'/#'+str(inst).zfill(2)+'/'+str(obsdate.year)+'/'+str(obsdate.month).zfill(2)+'/'+str(obsdate.day).zfill(2)+'/'
    fileroot = obsdate.strftime("%y%m%d")

    # photometer data
    startlocalday = dt.datetime.combine(obsdate,dt.time(4,0,0)) - dt.timedelta(hours=config.attrs['TimeZone'])
    endlocalday = dt.datetime.combine(obsdate,dt.time(21,0,0)) - dt.timedelta(hours=config.attrs['TimeZone'])
    sundata = s.loc[startlocalday:endlocalday]
    nobs = len(sundata)
    
    if nobs>0:     
        # get ozone for the day
        if ozoneopt:
            if sum(o3.index==obsdate): # daily ozone available
                ozonecolumn =  (o3.loc[o3.index==obsdate].ozone / 1000.0).iloc[0]
            else: # use monthly mean
                if verb: print('use monthly mean o3')
                o3m = o3.loc[(pd.to_datetime(o3.index).month==obsdate.month) & 
                             (pd.to_datetime(o3.index).year==obsdate.year)].ozone
                ozonecolumn = sum(o3m)/len(o3m) / 1000.0
        else:
            ozonecolumn = 0.25
    
        # Pressure data
        presfile = rootpath + 'output/' + filedir.replace('#','')  + fileroot + '.hpa'
        if os.path.isfile(presfile):
            p = fr.read_pressure_file(verb, presfile)
        else:
            p = []   
    
        # Black record
        blkfile = rootpath+'agsdat/'+filedir+fileroot+'.blk'
        if os.path.isfile(blkfile):
            if verb: print(f'blackfile : {blkfile}')
            blksun = read_black_record(blkfile,model)
        else:
            blksun = [0] * numchannels

        # Clock Fix
        if clockfix:  
            timecorr = clk.asof(dii).timecorr  
            if np.isnan(timecorr):
                timecorr = clk.timecorr.iloc[0]   # if prior to first entry in clk file, use first entry
      
        volt = np.zeros([nobs,numchannels,3])
        airmass = np.empty([nobs,numchannels,3])
        solzen = np.empty([nobs,3])
        solzenapp = np.empty([nobs,3])
        numsuntriples = [nobs,0,0] # total/morning/afternoon   

        aod = np.zeros([nobs, numchannels, 3])
        stdUaod = np.zeros([nobs, numchannels, 3])
        stdUangstrom48 = np.zeros([nobs, 3])
        Wvap = np.zeros([nobs, 3])
        aod1020c = np.zeros([nobs,3])
        filt = np.zeros(nobs) 
        
        for i in range(nobs): 
            obs = sundata.iloc[i]
            datetime = obs.name
            nday = (obsdate - cal.attrs['epoch']).days
            
            if clockfix: datetime = datetime - dt.timedelta(seconds=timecorr)
                
            if len(p)==0:
                p_i = c.defaultpressure
            else:
                p_i = p.asof(pd.Timestamp(datetime)).Pressure
                if math.isnan(p_i):   # before first pressure measurement use daily mean
                    p_i = p.Pressure.mean()
            
            for k in range(3):
                # Solar Position
                julday, timeday = sol.julian(datetime+dt.timedelta(seconds=(30*k)))            
                sunpos = sol.solar_position_almanac(julday,timeday)   
                DSunEarth=sunpos.DSunEarth
                solzen[i,k], solazi = sol.satellite_position(julday, timeday, stationlat*c.deg2rad, stationlon*c.deg2rad, sol.sun_elements)
                if k==1:
                    if (solazi*c.rad2deg)<180.0:
                        numsuntriples[1] +=1     #morning
                    else:
                        numsuntriples[2] +=1     #afternoon
                solzenapp[i,k] = sol.apparent_zenith(solzen[i,k] * c.rad2deg)
    
                ozoneOD = [x*ozonecolumn for x in c.ozonecoef[model-1]]   
                extraOD = [0]*numchannels 
    
                # Signal 
                for j, wl in enumerate(c.wavelength[model-1]):
                    
                    volt[i,j,k] = obs[f'ch{int(wl)}_{k}'] - blksun[j]
                    
                    rayleighOD = atm.rayleigh(wl, p_i)           
    
                    #Assume aod of 0.03 to calculate airmass. Refine later.
                    airmassODrayleigh = atm.getairmass(1,solzenapp[i,k]) * rayleighOD
                    airmassODaerosol = atm.getairmass(2,solzenapp[i,k]) * 0.03
                    airmassODozone = atm.getairmass(3,solzenapp[i,k]) * ozoneOD[j]
                    airmass[i,j,k] = ( airmassODrayleigh + airmassODaerosol + airmassODozone) \
                                    / (rayleighOD + 0.03 + ozoneOD[j])       
    
                    #print('volt = ',volt)
                    if volt[i,j,k]>0:
                        voltlog = math.log(volt[i,j,k]) 
                    else:
                        voltlog=-9.99   

                        
                    lnV0 = lnV0coef0[j] + lnV0coef1[j] * nday - 2.0 * math.log(DSunEarth)
                    aod[i,j,k] = (lnV0-voltlog) / airmass[i,j,k] - rayleighOD - ozoneOD[j] - extraOD[j]
    
                    stdUaod[i,j,k] =stdUaod1am[j] / airmass[i,j,k]
    
                    if j==iWV: 
                        wlratiolog = math.log(c.wavelength[model-1][i870] / c.wavelength[model-1][i670])
                        if (aod[i,i870,k]>0) & (aod[i,i670,k]>0):
                            angstrom68 = -math.log(aod[i,i870,k]/aod[i,i670,k])/wlratiolog
                        else:
                            angstrom68 = -9.99
                        aod[i,iWV,k] = aod[i,i870,k] * (c.wavelength[model-1][iWV] / c.wavelength[model-1][i870])**(-angstrom68)
                        Y = voltlog + airmass[i,iWV,k] * (rayleighOD + aod[i,iWV,k])
                        WVarg = (lnV0 - Y)/c.acoef
                        if WVarg>0:
                            Wvap[i,k] = WVarg**(1/c.bcoef)/airmass[i,iWV,k]
                            WvapOD = c.h2o_k * Wvap[i,k]**c.h2o_e   # Compute od of H2O at 1020nm.  Based on lowtran 7, probably wrong
                        else:
                            Wvap[i,k] = -9.99
                            WvapOD = -c.h2o_k * 9.99**c.h2o_e
       
                        aod1020c[i,k] = aod[i,i1020,k] - WvapOD
            
        wlratiolog = math.log(c.wavelength[model-1][i870] / c.wavelength[model-1][i440])
        stdUangstrom48 = (stdUaod[:,i440,:]/aod[:,i440,:] + stdUaod[:,i870,:]/aod[:,i870,:]) / wlratiolog
        stdUangstrom48[np.logical_or(aod[:,i440,:]<0, aod[:,i870,:]<0)] = -9.99
        
        if cirrusopt==0:
            cirrusflag=np.ones([nobs,3]) * -1
        elif cirrusopt==1:
            angstrom48 = np.ones(nobs) * -9.99
            angstrom48 = -np.log(aod[:,i870,:]/aod[:,i440,:], where=np.logical_and(aod[:,i870,:]>0, aod[:,i440,:]>0)) / np.log(c.wavelength[model-1][i870] / c.wavelength[model-1][i440])
            cirrusangstromthres = np.ones(nobs) * -9.99
            cirrusangstromthres = 0.862 + 0.556 * np.log10(aod[:,i440,:], where=aod[:,i440,:]>0)
            cirrusflag = angstrom48 < cirrusangstromthres
        elif cirrusopt==2:
            cirrusflag = aod[:,i1020,:] > aod[:,i870,:]
            

        #   1st pass filters
        # ------------------
        filt1 = np.zeros(nobs)  
        stackidx = np.zeros(nobs)
        numOKstacks=0
        for i in range(nobs):                 
            if verb: print(datetime)

            valid = True
            #Exclude if obs obtained before 0300 LST
            LSTcut = dt.datetime.combine(obsdate,dt.time(3,0,0))
            UTcut = LSTcut - dt.timedelta(hours=config.attrs['TimeZone'])
            if datetime < UTcut:
                if verb: print('before LST cutoff')
                continue        
            
            # Covariance of triplet is above specified limit            
            CVvolt = np.std(volt[i,i870,:], ddof=1) / np.mean(volt[i,i870,:])
            if CVvolt > CVmax:
                if verb: print(f'High Covariance {CVvolt:3.5f}')   
                continue
                
            # Sun too low                          
            if solzen[i,1] * c.rad2deg > solzenmax:
                if verb: print(f'Sun too low: zen={solzen[i,1] * c.rad2deg:3.1f} > {solzenmax}')
                valid = False
            if not valid:
                continue

            if verb: print(f'{datetime} GOOD OBS')   
            filt1[i]=1
            numOKstacks += 1
            stackidx[i]=numOKstacks                

        filt2 = filt1.copy()
        if (numOKstacks*3) > minobs:
            # 2nd pass filter
            # ---------------           
            stddev = np.std(aod[filt2==1,i870,:].flatten(), ddof=1) 
            if stddev > c.sd_crit:
                check3sd=True               
                nloop=0
                while check3sd:
                    nloop += 1
                    dailyaod = np.mean(aod[filt2==1,i870,:].flatten())
                    stddev = np.std(aod[filt2==1,i870,:].flatten(), ddof=1) 
                    nobs1 = sum(filt2==1)
                    if verb: print(f'[{nloop}] daily mean AOD 870 = {dailyaod:5.4f}   3*sd = {3*stddev:5.4f}    n = {nobs1}')
            
                    check3sd=False
                    for istack in np.unique(stackidx[stackidx>0]):
                        # stack should be within 3sd of daily mean
                        diff = abs(np.mean(aod[stackidx==istack,i870,:]) - dailyaod)
                        if  diff > 3*stddev:
                            if verb: print(f'{aod[stackidx==istack,i870,:]} outside 3sd')
                            filt2[stackidx==istack]=0
                            stackidx[stackidx==istack]=0
                            check3sd=True  
        
                    if nloop>numOKstacks: 
                        print('too many loops in 2nd pass filter')
                        break

            n1 = sum(filt2)
            filt2[np.all(abs(stdUaod[:,i440,:] / aod[:,i440,:]) > relUaod440thres, axis=1)] = 0
            n2 = sum(filt2)
            filt2[np.all(abs(stdUangstrom48) > Uangstrom48thres, axis=1)] = 0 
            if verb: print(f'{nobs} / {sum(filt1)} / {n1} / {n2} / {sum(filt2)}')
                             
            if sum(filt2)>0:            
#                if cirrusopt>0:
#                    aodfile = rootpath+'output/'+site+'/' + str(inst) + str(obsdate.year % 100).zfill(2) + str(obsdate.month).zfill(2) + str(obsdate.day).zfill(2) +'.aod_'+str(cirrusopt)
#                else:
                aodfile = rootpath+'output/'+site+'/' + str(inst) + str(obsdate.year % 100).zfill(2) + str(obsdate.month).zfill(2) + str(obsdate.day).zfill(2) +'.aod'
                os.makedirs(os.path.dirname(aodfile), exist_ok=True) 
                with open(aodfile, 'w') as f:
                    f.write('date       time    sunzen airmass '+
                            ' '.join([format(f'T{str(int(c.wavelength[model-1][x])).zfill(4)}', "6s") for x in langleych])+
                            ' W0936 T1020w Cirrus\n')
                    for ii in range(len(filt2)):
                        if filt2[ii]==1:
                            for kk in range(3):
                                
                                f.write((sundata.index[ii]+dt.timedelta(seconds=(30*kk))).strftime("%Y-%m-%d %H:%M:%S ")+
                                        format(solzen[ii,kk] * c.rad2deg,"5.2f")+" "+
                                        format(airmass[ii,i870,kk],"5.4f")+" "+
                                        " ".join([format(x, "5.4f") for x in aod[ii,langleych,kk]])+" "+
                                        format(Wvap[ii,kk],"5.4f")+" "+
                                        format(aod1020c[ii,kk],"5.4f")+" "+
                                                    str(int(cirrusflag[ii,kk]))+'\n')
                if verb: print(f'file wirtten to {aodfile}')                
        else:
            print(f'{obsdate} not enough obs passed filters')
    else:
        print(f'{obsdate} no file')
        

