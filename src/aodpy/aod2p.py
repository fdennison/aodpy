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
ozoneopt = aod2pconf['ozoneopt']
calfile_in =  aod2pconf['aod2p']['calfile_in']
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

ozonefile = rootpath+'ozone/'+site+'.o3'
if os.path.isfile(ozonefile) and ozoneopt:
    o3 = fr.read_bom_ozone2011(ozonefile)
else:
    print('using default ozone: 250 DU')

configfile = rootpath+'config/'+site+'.cfn'
config = fr.read_site_configuration(configfile, startdate)
stationlat = config.attrs['Latitude'] #-12.6607
stationlon = config.attrs['Longitude'] #132.8931
model = config.CimelModel
inst = config.CimelNumber

#cal = fr.read_cal(rootpath+'suncals/#'+calfile_in[1:-10].zfill(2)+'/'+calfile_in)
cal = fr.read_cal(rootpath+'suncals/'+str(inst).zfill(2)+'/'+calfile_in)
#See p 87 of 2006/7 notebook
#This kernel has to be divided by airmass to give U(aod)
usignal=0.003
if cal.attrs['numlangleys']>0:
    cal['stdUaod1am'] = np.sqrt(usignal**2 + cal.erms**2/cal.attrs['numlangleys'])
else:
    cal['stdUaod1am'] = usignal

if clockfix:
    calperiod_filename = str(inst) + str(startdate.year % 100).zfill(2) + str(startdate.month).zfill(2) + str(enddate.month).zfill(2)
    clk = fr.read_adj_file(rootpath + 'suncals/' + str(inst) +'/' + calperiod_filename + '.clk', 'timecorr')

i440 = 1 - 1
i670 = 2 - 1
i870 = 3 - 1
i1020 = 4 - 1
iWV = 10 - 1

numchannels = 10
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
    filename = rootpath+'agsdat/'+filedir+fileroot+'.lsu'
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
        #presfile = rootpath + 'agsdat/' + filedir + fileroot + '.hpa'  # date formats will be different
        #if os.path.isfile(presfile):
        #    p = fr.read_old_pressure_file(args.verbose, presfile)
        #else:
        #    p=[]
        presfile = rootpath + '/PyOut/' + filedir.replace('#','') + fileroot + '.hpa'
        p = fr.read_pressure_file(args.verbose, presfile)
    
    
        # photometer data
        sundata = fr.read_single_sun_record(filename, model)
        nobs = len(sundata)
    
        # Black record
        blkfile = rootpath+'agsdat/'+filedir+fileroot+'.blk'
        if os.path.isfile(blkfile):
            if args.verbose: print(f'blackfile : {blkfile}')
            blksun = read_black_record(blkfile,model)
        else:
            blksun = [0] * numchannels

        # Clock Fix
        if clockfix:  timecorr = round(clk.asof(dii).timecorr)

        volt = np.zeros([nobs, numchannels])
        datetime = []
        aod = np.zeros([nobs, numchannels])
        airmass = np.zeros([nobs, numchannels])
        stdUaod = np.zeros([nobs, numchannels])
        stdUangstrom48 = np.zeros(nobs)
        solzen = np.zeros(nobs)
        Wvap = np.zeros(nobs)
        aod1020c = np.zeros(nobs)
        filt = np.zeros(nobs)
        
        for iobs in range(nobs): 
            datetime.append(sundata.index[iobs])
            if clockfix: datetime[iobs] = datetime[iobs] - dt.timedelta(seconds=timecorr)
        
            julday, timeday = sol.julian(datetime[iobs])            
            sunpos = sol.solar_position_almanac(julday,timeday)   
            DSunEarth=sunpos.DSunEarth
            solzen[iobs], solazi = sol.satellite_position(julday, timeday, stationlat*deg2rad, stationlon*deg2rad, sol.sun_elements)
            solzenapp = sol.apparent_zenith(solzen[iobs] * rad2deg)             
            
            if len(p)==0:
                pr = defaultpressure
            else:
                pr = p.asof(datetime[iobs]).Pressure
                if math.isnan(pr):   # before first pressure measurement use daily mean
                    pr = p.Pressure.mean()             
    
            # Cal File
            nday = (obsdate - cal.attrs['epoch']).days
            numlangch = cal.attrs['numlangleychannels']
            lnV0_1AU=[0]*numchannels
            for n in range(numlangch):
                lnV0_1AU[n] = cal.iloc[n].lnV0coef0 +  cal.iloc[n].lnV0coef1 * nday   #
            n=iWV
            lnV0_1AU[n] = cal.iloc[n].lnV0coef0 +  cal.iloc[n].lnV0coef1 * nday # point of doing this channel seperately?
       
            # Signal
            ozoneOD = [x*ozonecolumn for x in atm.ozonecoef[model-1]]  
            extraOD = [0]*numchannels 

            for n in range(numchannels):
                rayleighOD = atm.rayleigh(atm.wavelength[model-1][n], pr)  
                #Assume aod of 0.03 to calculate airmass. Refine later.
                airmass[iobs,n] = (atm.getairmass(1,solzenapp) * rayleighOD + \
                                atm.getairmass(2,solzenapp) * 0.03 + \
                                atm.getairmass(3,solzenapp) * ozoneOD[n] ) \
                                / (rayleighOD + 0.03 + ozoneOD[n]) 

                volt[iobs,n] = sundata.iloc[iobs]['Ch'+str(n+1)] - blksun[n]
                if volt[iobs,n]>0:
                    voltlog = math.log(volt[iobs,n]) 
                else:
                    voltlog=-9.99        

                lnV0 = lnV0_1AU[n] - 2.0 * math.log(DSunEarth)

                aod[iobs,n] = (lnV0-voltlog) / airmass[iobs,n] - rayleighOD - ozoneOD[n] - extraOD[n]

                stdUaod[iobs,n] = cal.iloc[n].stdUaod1am / airmass[iobs,n]

                if n==iWV: 
                    wlratiolog = math.log(atm.wavelength[model-1][i870] / atm.wavelength[model-1][i670])
                    if (aod[iobs,i870]>0) & (aod[iobs,i670]>0):
                        angstrom68 = -math.log(aod[iobs,i870]/aod[iobs,i670])/wlratiolog
                    else:
                        angstrom68 = -9.99
                    aod[iobs,iWV] = aod[iobs,i870] * (atm.wavelength[model-1][iWV] / atm.wavelength[model-1][i870])**(-angstrom68)
                    Y = voltlog + airmass[iobs,iWV] * (rayleighOD + aod[iobs,iWV])
                    WVarg = (lnV0 - Y)/acoef
                    if WVarg>0:
                        Wvap[iobs] = WVarg**(1/bcoef)/airmass[iobs,iWV]
                        WvapOD = h2o_k * Wvap[iobs]**h2o_e   # Compute od of H2O at 1020nm.  Based on lowtran 7, probably wrong
                    else:
                        Wvap[iobs] = -9.99
                        WvapOD = -h2o_k * 9.99**h2o_e
   
                    aod1020c[iobs] = aod[iobs,i1020] - WvapOD
            
        wlratiolog = math.log(atm.wavelength[model-1][i870] / atm.wavelength[model-1][i440])
        stdUangstrom48 = (stdUaod[:,i440]/aod[:,i440] + stdUaod[:,i870]/aod[:,i870]) / wlratiolog
        stdUangstrom48[np.logical_or(aod[:,i440]<0, aod[:,i870]<0)] = -9.99
        
        if cirrusopt==0:
            cirrusflag=np.ones(nobs) * -1
        elif cirrusopt==1:
            angstrom48 = np.ones(nobs) * -9.99
            angstrom48 = -np.log(aod[:,i870]/aod[:,i440], where=np.logical_and(aod[:,i870]>0, aod[:,i440]>0)) / np.log(atm.wavelength[model-1][i870] / atm.wavelength[model-1][i440])
            cirrusangstromthres = np.ones(nobs) * -9.99
            cirrusangstromthres = 0.862 + 0.556 * np.log10(aod[:,i440], where=aod[:,i440]>0)
            cirrusflag = angstrom48 < cirrusangstromthres
        elif cirrusopt==2:
            cirrusflag = aod[:,i1020] > aod[:,i870]
            
            
        datetime = np.array(datetime)

        #   1st pass filters
        # ------------------
        filt1 = np.zeros(nobs)  
        stackidx = np.zeros(nobs)
        numOKstacks=0
        iobs=0
        nloop=0
        while iobs<=(nobs-stacksize):
            nloop += 1
            if nloop>1000:
                break
                
            if args.verbose: print(datetime[iobs])

            valid = True

            # stack of obs does not fit in specified time window
            if ((datetime[iobs+stacksize-1] - datetime[iobs]).total_seconds())/60.0 > window:         
                if args.verbose: print('stack of obs does not fit in the {} minute window : {} - {}'.format(
                               window,datetime[iobs], datetime[iobs+stacksize-1]))
                iobs += 1    
                continue
                
            # temperature/signal to high or low - assessed in read function
            for i in range(stacksize):
                if not sundata.iloc[iobs+i].Valid==0:
                    if args.verbose: print(f'Stack not valid : {sundata.iloc[iobs+i].Valid} [{iobs}]+[{i}]')
                    iobs += 1
                    valid = False
                    break
            if not valid:
                continue

            #Exclude if obs obtained before 0300 LST
            LSTcut = dt.datetime.combine(obsdate,dt.time(3,0,0))
            UTcut = LSTcut - dt.timedelta(hours=config.attrs['TimeZone'])
            if datetime[iobs] < UTcut:
                if args.verbose: print('before LST cutoff')
                iobs += 1
                continue        
            
            # Covariance of stack is above specified limit            
            CVvolt = np.std(volt[iobs:(iobs+stacksize),i870], ddof=1) / np.mean(volt[iobs:(iobs+stacksize),i870])
            if CVvolt > CVmax:
                if args.verbose: 
                    [print(d) for d in datetime[(iobs+1):(iobs+stacksize)]]
                    print(f'High Covariance {CVvolt:3.5f}')  
                iobs += 3   
                continue
                
            # Sun too low                          
            for i in range(stacksize):
                if solzen[iobs+i] * rad2deg > solzenmax:
                    if args.verbose: print(f'Sun too low: zen={solzen[iobs+i] * rad2deg:3.1f} > {solzenmax}')
                    iobs += 1
                    valid = False
                    break
            if not valid:
                continue


            if args.verbose: 
                [print(d) for d in datetime[(iobs+1):(iobs+stacksize)]]
                print('GOOD STACK')   
            filt1[iobs:(iobs+stacksize)]=1
            numOKstacks += 1
            stackidx[iobs:(iobs+stacksize)]=numOKstacks                
            iobs += stacksize

        filt2 = filt1.copy()
        if (numOKstacks*stacksize) > minobs:
            # 2nd pass filter
            # ---------------           
            stddev = np.std(aod[filt2==1,i870], ddof=1) 
            if stddev > sd_crit:
                check3sd=True               
                nloop=0
                while check3sd:
                    nloop += 1
                    dailyaod = np.mean(aod[filt2==1,i870])
                    stddev = np.std(aod[filt2==1,i870], ddof=1) 
                    nobs1 = sum(filt2==1)
                    if args.verbose: print(f'[{nloop}] daily mean AOD 870 = {dailyaod:5.4f}   3*sd = {3*stddev:5.4f}    n = {nobs1}')
            
                    check3sd=False
                    for istack in np.unique(stackidx[stackidx>0]):
                        # stack should be within 3sd of daily mean
                        diff = abs(np.mean(aod[stackidx==istack,i870]) - dailyaod)
                        if args.verbose:print('[{}  {}  {}]  diff = {:5.4f}'.format(*datetime[stackidx==istack], diff))
                        if  diff > 3*stddev:
                            if args.verbose: print('[{}  {}  {}] outside 3sd'.format(*datetime[stackidx==istack]))
                            filt2[stackidx==istack]=0
                            stackidx[stackidx==istack]=0
                            check3sd=True  
        
                    if nloop>numOKstacks: 
                        print('too many loops in 2nd pass filter')
                        break

            n1 = sum(filt2)
            filt2[abs(stdUaod[:,i440] / aod[:,i440]) > relUaod440thres] = 0
            n2 = sum(filt2)
            filt2[abs(stdUangstrom48) > Uangstrom48thres] = 0 
            if args.verbose: print(f'{nobs} / {sum(filt1)} / {n1} / {n2} / {sum(filt2)}')
                             
            if sum(filt2)>0:            
                out = " ".join([format(x, "6.4f") for x in np.mean(aod[filt2==1,:numlangch], axis=0)])
                if args.verbose: print(f'{obsdate}, '+out)
                if cirrusopt>0:
                    aodfile = rootpath+'PyOut/'+site+'/' + str(inst) + str(obsdate.year % 100).zfill(2) + str(obsdate.month).zfill(2) + str(obsdate.day).zfill(2) +'.aod_'+str(cirrusopt)
                else:
                    aodfile = rootpath+'PyOut/'+site+'/' + str(inst) + str(obsdate.year % 100).zfill(2) + str(obsdate.month).zfill(2) + str(obsdate.day).zfill(2) +'.aod'
                os.makedirs(os.path.dirname(aodfile), exist_ok=True) 
                with open(aodfile, 'w') as f:
                    f.write('date       time    sunzen airmass '+
                            ' '.join([format(f'T{str(int(x)).zfill(4)}', "6s") for x in atm.wavelength[model-1][:numlangch]])+
                            ' W0936 T1020w Cirrus\n')
                    for ii in range(len(solzen)):
                        if filt2[ii]==1:
                            f.write(sundata.index[ii].strftime("%Y-%m-%d %H:%M:%S ")+
                                    format(solzen[ii] * rad2deg,"5.2f")+" "+
                                    format(airmass[ii,i870],"5.4f")+" "+
                                    " ".join([format(x, "5.4f") for x in aod[ii,:numlangch]])+" "+
                                    format(Wvap[ii],"5.4f")+" "+
                                    format(aod1020c[ii],"5.4f")+" "+
                                    str(int(cirrusflag[ii]))+'\n')
                                    
        else:
            if args.verbose: print(f'{obsdate} not enough obs passed filters')
    else:
        print(f'{obsdate} no lsu file')
        

