import datetime as dt
import pandas as pd
import numpy as np
import os as os
import math
import fileread_module as fr
import solpos_module as sol
import calfit_module as fit
import atmos_module as atm
import tomllib
import argparse

pi = math.pi
rad2deg = 180.0 / pi
deg2rad = pi /180.0

parser=argparse.ArgumentParser(description="Calculate AOD from sun photometer measuments (.lsu files)")
parser.add_argument('-t', '--test', dest='configtoml', action='store_const', 
                    const='testconfig.toml', default='config.toml',
                    help='run on sample data and check output')
parser.add_argument('-v', '--verbose',action='store_true')
args=parser.parse_args()

# aod2p config options
with open(args.configtoml, "rb") as f:
    aod2pconf = tomllib.load(f)
site = aod2pconf['site']
rootpath = aod2pconf['rootpath']
startdate = aod2pconf['startdate']
enddate = aod2pconf['enddate']
ozonedir = aod2pconf['ozonedir']
calfile =  aod2pconf['aod2p']['calfile']
clockfix = aod2pconf['aod2p']['clockfix']
cirrusopt = aod2pconf['aod2p']['cirrusopt']
window = aod2pconf['aod2p']['window']
CVmax = aod2pconf['aod2p']['CVmax']
solzenmax = aod2pconf['aod2p']['solzenmax']
minobs = aod2pconf['aod2p']['minobs']
sd_crit = aod2pconf['aod2p']['sd_crit']
relUaod440thres = aod2pconf['aod2p']['relUaod440thres']
Uangstrom48thres = aod2pconf['aod2p']['Uangstrom48thres']
print(f'Run aod2p for {site} from {startdate} to {enddate}')
if args.verbose: print(aod2pconf) 

ozonefile = rootpath+ozonedir+site+'.o3'
if os.path.isfile(ozonefile):
    o3 = fr.read_bom_ozone2011(ozonefile)
    ozoneopt = True
else:
    ozoneopt = False
    print('using default ozone: 250 DU')

configfile = rootpath+'config/'+site+'.cfn'
config = fr.read_site_configuration(configfile, startdate)
stationlat = config.attrs['Latitude'] #-12.6607
stationlon = config.attrs['Longitude'] #132.8931
model = config.CimelModel
inst = config.CimelNumber

calperiod_filename = str(inst) + str(startdate.year % 100).zfill(2) + str(startdate.month).zfill(2) + str(enddate.month).zfill(2)
if len(calfile)==3:  # 
    calfile = calperiod_filename+'.'+calfile
cal = fr.read_cal(rootpath+'suncals/'+str(inst).zfill(2)+'/'+calfile)

#See p 87 of 2006/7 notebook
#This kernel has to be divided by airmass to give U(aod)
usignal=0.003
if cal.attrs['numlangleys']>0:
    cal['stdUaod1am'] = np.sqrt(usignal**2 + cal.erms**2/cal.attrs['numlangleys'])
else:
    cal['stdUaod1am'] = usignal

ozonefile = rootpath+ozonedir+site+'.o3'
if os.path.isfile(ozonefile) and ozoneopt:
    o3 = fr.read_bom_ozone2011(ozonefile)
else:
    print('using default ozone: 250 DU')

if clockfix:  clk = fr.read_adj_file(rootpath + 'suncals/' + str(inst).zfill(2) +'/' + calperiod_filename + '.clk', 'timecorr')

if model<=9:   # later models the 2nd channel is 675nm
    lambda2 = 670
else:
    lambda2 = 675

if model == 10:
    ext = '.NSU'
else:    
    ext = '.sun'
    
i440 = atm.wavelength[model-1].index(440)  # 1-1
i670 = atm.wavelength[model-1].index(lambda2)  #2 - 1
i870 = atm.wavelength[model-1].index(870)  #3 - 1
i1020 = atm.wavelength[model-1].index(1020)#4 - 1
iWV = atm.wavelength[model-1].index(936)   #10 - 1
numchannels = len(atm.wavelength[model-1])

stacksize=3
acoef = 0.6548
bcoef = 0.574
h2o_k = 0.9156e-2
h2o_e = 0.5918
defaultpressure = 1013.0 # hPa

datelist = pd.date_range(startdate,enddate,freq='d')  
for dii in datelist:
    obsdate = pd.to_datetime(dii).date()
    if args.verbose: print('----'+obsdate.strftime("%Y/%m/%d")+'----')
    filedir = site+'/#'+str(inst).zfill(2)+'/'+str(obsdate.year)+'/'+str(obsdate.month).zfill(2)+'/'+str(obsdate.day).zfill(2)+'/'
    fileroot = obsdate.strftime("%y%m%d")
    filename = rootpath+'agsdat/'+filedir+fileroot+ext

    if os.path.isfile(filename):
     
        # get ozone for the day
        if ozoneopt:
            if sum(o3['date']==obsdate): # daily ozone available
                ozonecolumn =  (o3.loc[o3['date']==obsdate].ozone / 1000.0).iloc[0]
            else: # use monthly mean
                if args.verbose: print('use monthly mean o3')
                o3m = o3.loc[(pd.to_datetime(o3['date']).dt.month==obsdate.month) & 
                             (pd.to_datetime(o3['date']).dt.year==obsdate.year)].ozone
                ozonecolumn = sum(o3m)/len(o3m) / 1000.0
        else:
            ozonecolumn = 0.25
    
        # Pressure data
        presfile = rootpath + 'PyOut/' + filedir.replace('#','')  + fileroot + '.hpa'
        p = fr.read_pressure_file(args.verbose, presfile)    
    
        # photometer data
        sundata = fr.read_triple_sun_record(filename, model)
        nobs = len(sundata)
    
        # Black record
        blkfile = rootpath+'agsdat/'+filedir+fileroot+'.blk'
        if os.path.isfile(blkfile):
            if args.verbose: print(f'blackfile : {blkfile}')
            blksun = fr.read_black_record(blkfile,model)
        else:
            blksun = [0] * numchannels

        # Clock Fix
        if clockfix:  timecorr = round(clk.asof(dii).timecorr)

       
        volt = np.zeros([nobs,3,numchannels])
        airmass  = np.empty([nobs,3,numchannels])
        airmassODrayleigh = np.empty([nobs,3,numchannels])
        airmassODaerosol  = np.empty([nobs,3,numchannels])
        airmassODozone  = np.empty([nobs,3,numchannels]) 
        solzen  = np.empty([nobs,3])
        solzenapp  = np.empty([nobs,3])
        numsuntriples = [nobs,0,0] # total/morning/afternoon

        aod = np.zeros([nobs, 3, numchannels])
        stdUaod = np.zeros([nobs, 3, numchannels])
        stdUangstrom48 = np.zeros([nobs, 3])
        Wvap = np.zeros([nobs, 3])
        aod1020c = np.zeros([nobs,3])
        filt = np.zeros(nobs) 
        
        for k in range(nobs): 
    
            obs = sundata.iloc[k]
    
            datetime = obs.name
            if clockfix: datetime = datetime - dt.timedelta(seconds=timecorr)
                
            if len(p)==0:
                pr = defaultpressure
            else:
                pr = p.asof(pd.Timestamp(datetime)).Pressure
                if math.isnan(pr):   # before first pressure measurement use daily mean
                    pr = p.Pressure.mean()

            # Cal File
            nday = (obsdate - cal.attrs['epoch']).days
            numlangch = cal.attrs['numlangleychannels']
            lnV0_1AU=[0]*numchannels
            for n in range(numchannels):
                lnV0_1AU[n] = cal.iloc[n].lnV0coef0 +  cal.iloc[n].lnV0coef1 * nday   #
            
            for i in range(3):
            # Solar Position
                julday, timeday = sol.julian(datetime+dt.timedelta(seconds=(30*i)))            
                sunpos = sol.solar_position_almanac(julday,timeday)   
                DSunEarth=sunpos.DSunEarth
                solzen[k,i], solazi = sol.satellite_position(julday, timeday, stationlat*deg2rad, stationlon*deg2rad, sol.sun_elements)
                if i==1:
                    if (solazi*rad2deg)<180.0:
                        numsuntriples[1] +=1     #morning
                    else:
                        numsuntriples[2] +=1     #afternoon
                solzenapp[k,i] = sol.apparent_zenith(solzen[k,i] * rad2deg)
    
                ozoneOD = [x*ozonecolumn for x in atm.ozonecoef[model-1]]   
                extraOD = [0]*numchannels 
    
                # Signal 
                for n, wl in enumerate(atm.wavelength[model-1]):
                    
                    volt[k,i,n] = obs['ch'+str(int(wl))][i] - blksun[n]
                    
                    rayleighOD = atm.rayleigh(wl, pr)           
    
                    #Assume aod of 0.03 to calculate airmass. Refine later.
                    airmassODrayleigh[k,i,n] = atm.getairmass(1,solzenapp[k,i]) * rayleighOD
                    airmassODaerosol[k,i,n] = atm.getairmass(2,solzenapp[k,i]) * 0.03
                    airmassODozone[k,i,n] = atm.getairmass(3,solzenapp[k,i]) * ozoneOD[n]
                    airmass[k,i,n] = ( airmassODrayleigh[k,i,n] + airmassODaerosol[k,i,n] + airmassODozone[k,i,n]) \
                                    / (rayleighOD + 0.03 + ozoneOD[n])       
    
                    #print('volt = ',volt)
                    if volt[k,i,n]>0:
                        voltlog = math.log(volt[k,i,n]) 
                    else:
                        voltlog=-9.99   

                        
                    lnV0 = lnV0_1AU[n] - 2.0 * math.log(DSunEarth)
                    aod[k,i,n] = (lnV0-voltlog) / airmass[k,i,n] - rayleighOD - ozoneOD[n] - extraOD[n]
    
                    stdUaod[k,i,n] = cal.iloc[n].stdUaod1am / airmass[k,i,n]
    
                    if n==iWV: 
                        wlratiolog = math.log(atm.wavelength[model-1][i870] / atm.wavelength[model-1][i670])
                        if (aod[k,i,i870]>0) & (aod[k,i,i670]>0):
                            angstrom68 = -math.log(aod[k,i,i870]/aod[k,i,i670])/wlratiolog
                        else:
                            angstrom68 = -9.99
                        aod[k,i,iWV] = aod[k,i,i870] * (atm.wavelength[model-1][iWV] / atm.wavelength[model-1][i870])**(-angstrom68)
                        Y = voltlog + airmass[k,i,iWV] * (rayleighOD + aod[k,i,iWV])
                        WVarg = (lnV0 - Y)/acoef
                        if WVarg>0:
                            Wvap[k,i] = WVarg**(1/bcoef)/airmass[k,i,iWV]
                            WvapOD = h2o_k * Wvap[k,i]**h2o_e   # Compute od of H2O at 1020nm.  Based on lowtran 7, probably wrong
                        else:
                            Wvap[k,i] = -9.99
                            WvapOD = -h2o_k * 9.99**h2o_e
       
                        aod1020c[k,i] = aod[k,i,i1020] - WvapOD
            
        wlratiolog = math.log(atm.wavelength[model-1][i870] / atm.wavelength[model-1][i440])
        stdUangstrom48 = (stdUaod[:,:,i440]/aod[:,:,i440] + stdUaod[:,:,i870]/aod[:,:,i870]) / wlratiolog
        stdUangstrom48[np.logical_or(aod[:,:,i440]<0, aod[:,:,i870]<0)] = -9.99
        
        if cirrusopt==0:
            cirrusflag=np.ones([nobs,3]) * -1
        elif cirrusopt==1:
            angstrom48 = np.ones(nobs) * -9.99
            angstrom48 = -np.log(aod[:,:,i870]/aod[:,:,i440], where=np.logical_and(aod[:,:,i870]>0, aod[:,:,i440]>0)) / np.log(atm.wavelength[model-1][i870] / atm.wavelength[model-1][i440])
            cirrusangstromthres = np.ones(nobs) * -9.99
            cirrusangstromthres = 0.862 + 0.556 * np.log10(aod[:,:,i440], where=aod[:,:,i440]>0)
            cirrusflag = angstrom48 < cirrusangstromthres
        elif cirrusopt==2:
            cirrusflag = aod[:,:,i1020] > aod[:,:,i870]
            
            
        datetime = np.array(datetime)

        #   1st pass filters
        # ------------------
        filt1 = np.zeros(nobs)  
        stackidx = np.zeros(nobs)
        numOKstacks=0
        for k in range(nobs): 
                
            if args.verbose: print(datetime)

            valid = True
                
            # temperature/signal to high or low - assessed in read function            
            if not sundata.iloc[k].Valid==0:
                if args.verbose: print(f'Stack not valid : {sundata.iloc[k].Valid} ')
                valid = False
            if not valid:
                continue

            #Exclude if obs obtained before 0300 LST
            LSTcut = dt.datetime.combine(obsdate,dt.time(3,0,0))
            UTcut = LSTcut - dt.timedelta(hours=config.attrs['TimeZone'])
            if datetime < UTcut:
                if args.verbose: print('before LST cutoff')
                continue        
            
            # Covariance of stack is above specified limit            
            CVvolt = np.std(volt[k,:,i870], ddof=1) / np.mean(volt[k,:,i870])
            if CVvolt > CVmax:
                if args.verbose: print(f'High Covariance {CVvolt:3.5f}')   
                continue
                
            # Sun too low                          
            if solzen[k,1] * rad2deg > solzenmax:
                if args.verbose: print(f'Sun too low: zen={solzen[k,1] * rad2deg:3.1f} > {solzenmax}')
                valid = False
            if not valid:
                continue


            if args.verbose: 
                print(f'{datetime} GOOD STACK')   
            filt1[k]=1
            numOKstacks += 1
            stackidx[k]=numOKstacks                

        filt2 = filt1.copy()
        if (numOKstacks*stacksize) > minobs:
            # 2nd pass filter
            # ---------------           
            stddev = np.std(aod[filt2==1,:,i870].flatten(), ddof=1) 
            if stddev > sd_crit:
                check3sd=True               
                nloop=0
                while check3sd:
                    nloop += 1
                    dailyaod = np.mean(aod[filt2==1,:,i870].flatten())
                    stddev = np.std(aod[filt2==1,:,i870].flatten(), ddof=1) 
                    nobs1 = sum(filt2==1)
                    if args.verbose: print(f'[{nloop}] daily mean AOD 870 = {dailyaod:5.4f}   3*sd = {3*stddev:5.4f}    n = {nobs1}')
            
                    check3sd=False
                    for istack in np.unique(stackidx[stackidx>0]):
                        # stack should be within 3sd of daily mean
                        diff = abs(np.mean(aod[stackidx==istack,:,i870]) - dailyaod)
                        if  diff > 3*stddev:
                            if args.verbose: print('[{}  {}  {}] outside 3sd'.format(*datetime[stackidx==istack]))
                            filt2[stackidx==istack]=0
                            stackidx[stackidx==istack]=0
                            check3sd=True  
        
                    if nloop>numOKstacks: 
                        print('too many loops in 2nd pass filter')
                        break

            n1 = sum(filt2)
            filt2[np.all(abs(stdUaod[:,:,i440] / aod[:,:,i440]) > relUaod440thres, axis=1)] = 0
            n2 = sum(filt2)
            filt2[np.all(abs(stdUangstrom48) > Uangstrom48thres, axis=1)] = 0 
            if args.verbose: print(f'{nobs} / {sum(filt1)} / {n1} / {n2} / {sum(filt2)}')
                             
            if sum(filt2)>0:            
                if cirrusopt>0:
                    aodfile = rootpath+'PyOut/'+site+'/' + str(inst) + str(obsdate.year % 100).zfill(2) + str(obsdate.month).zfill(2) + str(obsdate.day).zfill(2) +'.aod_'+str(cirrusopt)
                else:
                    aodfile = rootpath+'PyOut/'+site+'/' + str(inst) + str(obsdate.year % 100).zfill(2) + str(obsdate.month).zfill(2) + str(obsdate.day).zfill(2) +'.aod'
                os.makedirs(os.path.dirname(aodfile), exist_ok=True) 
                with open(aodfile, 'w') as f:
                    f.write('date       time    sunzen airmass '+
                            ' '.join([format(f'T{str(int(x)).zfill(4)}', "6s") for x in atm.wavelength[model-1][:numlangch]])+
                            ' W0936 T1020w Cirrus\n')
                    for ii in range(len(filt2)):
                        if filt2[ii]==1:
                            for jj in range(3):
                                
                                f.write((sundata.index[ii]+dt.timedelta(seconds=(30*jj))).strftime("%Y-%m-%d %H:%M:%S ")+
                                        format(solzen[ii,jj] * rad2deg,"5.2f")+" "+
                                        format(airmass[ii,jj,i870],"5.4f")+" "+
                                        " ".join([format(x, "5.4f") for x in aod[ii,jj,:numlangch]])+" "+
                                        format(Wvap[ii,jj],"5.4f")+" "+
                                        format(aod1020c[ii,jj],"5.4f")+" "+
                                        str(int(cirrusflag[ii,jj]))+'\n')
                                    
        else:
            if args.verbose: print(f'{obsdate} not enough obs passed filters')
    else:
        if args.verbose: print(f'{obsdate} no lsu file')
        
print('done')
