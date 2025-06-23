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

parser=argparse.ArgumentParser(description="Apply general method to sun photometer triplet measuments (.sun files)")
parser.add_argument('-t', '--test', dest='configtoml', action='store_const', 
                    const='testconfig.toml', default='config.toml',
                    help='run on sample data and check output')
parser.add_argument('-v', '--verbose',action='store_true')
args=parser.parse_args()

# langley config options
with open(args.configtoml, "rb") as f:
    genconf = tomllib.load(f)
site = genconf['site']
rootpath = genconf['rootpath']
startdate = genconf['startdate']
enddate = genconf['enddate']
ozonedir = genconf['ozonedir']
fitV0 = genconf['cal']['fitV0']
clockfix = genconf['cal']['clockfix']
MaxSdevFit = genconf['cal']['MaxSdevFit']
UseRefChannel4QA = genconf['cal']['UseRefChannel4QA']
calfile_in =  genconf['cal']['calfile_in']
refchan = genconf['cal']['refchannel']
print(f'Run gencal for {site} from {startdate} to {enddate}, ref chan = {refchan}')
verb = args.verbose
if verb: print(genconf)

# site config
configfile = rootpath+'config/'+site+'.cfn'
config = fr.read_site_configuration(configfile, startdate)
stationlat = config.attrs['Latitude'] 
stationlon = config.attrs['Longitude'] 
model = config.CimelModel
inst = config.CimelNumber
numchannels = len(c.wavelength[model-1])

caloutext = '.'+str(refchan)
iref = c.wl_index(refchan,c.wavelength[model-1]) 

calperiod_filename = str(inst) + str(startdate.year % 100).zfill(2) + str(startdate.month).zfill(2) + str(enddate.month).zfill(2)

# cut file
cutfile = rootpath + 'suncals/' + str(inst).zfill(2) +'/' + calperiod_filename + '.cut'
if os.path.isfile(cutfile):
    cut = fr.read_adj_file(cutfile, 'AmPm')
    usecut = True
else:
    usecut = False
    print('no cut file')

if clockfix: clk = fr.read_adj_file(rootpath + 'clockdrift/' + site + '.clk', 'timecorr')

if len(calfile_in)==3:  # if only lcl or wavelength then use default naming convention, otherwise must be whole filename
    calfile_in = calperiod_filename +'.'+ calfile_in
if verb: print(calfile_in)    
cal = fr.read_cal(rootpath+'suncals/'+calfile_in[0:-10].zfill(2)+'/'+calfile_in)
numlangleychannels = cal.attrs['numlangleychannels']
calepoch = cal.attrs['epoch']
numgencyc = cal.attrs['numgencycles'] + 1

# ozone file
ozonefile = rootpath+ozonedir+site+'.o3'
if os.path.isfile(ozonefile):
    o3 = fr.read_ozone(ozonefile)
    ozoneopt=True
else:
    ozoneopt=False
    print('using default ozone: 250 DU')
   
i870 = c.wl_index(870,c.wavelength[model-1])    
iWV = c.wl_index(936,c.wavelength[model-1])   

s = fr.read_all_triple_sun_records(rootpath,site,startdate,enddate,inst,model)

numlangleys=0
dayssinceepoch=[]
lnV0glob=[]
ampm = ['Am','Pm']
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
        blksun = fr.read_black_record(blkfile,model)
    else:
        blksun = [0] * numchannels

    # Cal File In
    nday = (obsdate - cal.attrs['epoch']).days
    lnV0in_1AU = cal.iloc[iref].lnV0coef0 +  cal.iloc[iref].lnV0coef1 * nday
    
    # Clock Fix
    if clockfix:  
        timecorr = clk.asof(dii).timecorr  
        if np.isnan(timecorr):
            timecorr = clk.timecorr.iloc[0]   # if prior to first entry in clk file, use first entry
  
    volt = np.zeros([nobs,numchannels,3])
    voltlog = np.empty([nobs,numchannels,3])
    airmass  = np.empty([nobs,numchannels,3])
    airmassODrayleigh = np.empty([nobs,numchannels,3])
    airmassODozone = np.empty([nobs,numchannels,3])
    solzenapp = np.empty([nobs,3])
    numsuntriples = [nobs,0,0] # total/morning/afternoon   
    pr=[]
    for i in range(nobs): 
        obs = sundata.iloc[i]
        datetime = obs.name
        
        if clockfix: datetime = datetime - dt.timedelta(seconds=timecorr)
        
        if len(p)==0:
            p_i = c.defaultpressure
        else:
            p_i = p.asof(pd.Timestamp(datetime)).Pressure
            if math.isnan(p_i):   # before first pressure measurement use daily mean
                p_i = p.Pressure.mean()
        pr.append(p_i)

        ozoneOD = [x*ozonecolumn for x in c.ozonecoef[model-1]]
        
        for k in range(3):
            # Solar Position
            julday, timeday = sol.julian(datetime+dt.timedelta(seconds=(30*k)))            
            sunpos = sol.solar_position_almanac(julday,timeday) 
            if (i==0) & (k==0): DSunEarth=sunpos.DSunEarth
            solzen, solazi = sol.satellite_position(julday, timeday, stationlat*c.deg2rad, stationlon*c.deg2rad, sol.sun_elements)
            if k==1:
                if (solazi*c.rad2deg)<180.0:
                    numsuntriples[1] +=1     #morning
                else:
                    numsuntriples[2] +=1     #afternoon
            solzenapp[i,k] = sol.apparent_zenith(solzen * c.rad2deg)

            # Signal 
            for j, wl in enumerate(c.wavelength[model-1]):
                
                volt[i,j,k] = obs[f'ch{int(c.wavelength[model-1][j])}_{k}'] - blksun[j]
                
                rayleighOD = atm.rayleigh(wl, p_i)           

                # [F90: Assume aod of 0.03 to calculate airmass. Refine later.]
                airmassODrayleigh[i,j,k] = atm.getairmass(1,solzenapp[i,k]) * rayleighOD
                airmassODaerosol = atm.getairmass(2,solzenapp[i,k]) * 0.03
                airmassODozone[i,j,k] = atm.getairmass(3,solzenapp[i,k]) * ozoneOD[j]
                airmass[i,j,k] = ( airmassODrayleigh[i,j,k] + airmassODaerosol + airmassODozone[i,j,k]) \
                                / (rayleighOD + 0.03 + ozoneOD[j])       

                if volt[i,j,k]>0:
                    voltlog[i,j,k] = math.log(volt[i,j,k]) 
                else:
                    voltlog[i,j,k]=-9.99

    tripletmean = np.mean(volt, axis=2)
    tripletCv = np.divide(np.std(volt, axis=2, ddof=1), tripletmean, out=np.ones_like(tripletmean)*99.99, where=tripletmean != 0)  * 100.0
  
    if verb: print('General method, processing...')
    nstart = [0,0]
    nend = [0,0]
    for iap in [1,2]:
        nstart[iap-1] = numsuntriples[1] * (iap - 1)
        nend[iap-1] = nstart[iap-1] + numsuntriples[iap]
    if verb: print(f'total:{numsuntriples[0]} triples, {numsuntriples[1]} am [{nstart[0]}:{nend[0]}], {numsuntriples[2]} pm [{nstart[1]}:{nend[1]}]')
    
    lnV0in=lnV0in_1AU-2.0*math.log(DSunEarth)
    MinSunTriples = 10
    for iap in [0,1]: #  morning / afternoon
        include=True
        # cutList
        if usecut:
            if dii in cut.index:
                if ampm[iap] in cut.loc[dii].AmPm:
                    include=False
                    print(f'cut {obsdate} {ampm[iap]}')
        
        
        if include:
            if numsuntriples[iap]>=MinSunTriples :
                SpreadFlag, numOk, indOk = fit.CheckTripletCv(verb, numsuntriples[iap+1], nstart[iap],airmass[:,1,i870],tripletCv[:,i870])
              
                if (numOk>=MinSunTriples) & SpreadFlag :
                    if verb: print(' Passed triplet cv langley filter')

                    V = np.empty([numOk*3])
                    W = np.empty([numOk*3])
                    X = np.empty([numOk*3,numchannels])
                    Y = np.empty([numOk*3,numchannels])
                    Z = np.empty([numOk*3,numchannels])
                    ii=0 # Singlet index       
                    for i in range(numOk):              # Ok Triplet index
                        for k in range(3):                     # Triplet components
                            X[ii,:]=airmass[indOk[i],:,k]
                            Y[ii,:]=voltlog[indOk[i],:,k]
                            V[ii]=atm.getairmass(2,solzenapp[indOk[i],k])
                            #W=m*aod at reference wavelength
                            W[ii]=lnV0in - voltlog[indOk[i],iref,k] - airmassODrayleigh[indOk[i],iref,k] - airmassODozone[indOk[i],iref,k]
                            #Z=Dependent variable for regression
                            Z[ii,:]=voltlog[indOk[i],:,k] + airmassODrayleigh[indOk[i],:,k] + airmassODozone[indOk[i],:,k]
                            ii=ii+1
                    npts_ap=ii
                     
                    if  UseRefChannel4QA:
                        QARef=iref
                    else: 
                        QARef=i870
    
                    SpreadFlag,FitFlag,indOk,npts_ap = fit.CheckFitQuality(verb,npts_ap,X,Y,QARef,MaxSdevFit)
                    
                    V=V[indOk]
                    W=W[indOk]
                    Z=Z[indOk,:]    
                    if SpreadFlag & FitFlag: 
                        if verb: print(' Passed fit quality filter')
                        weight = [1.0] * npts_ap
                        intercept = np.zeros(numchannels)  
                        slope = np.zeros(numchannels)
                        for j in range(numchannels):
                            intercept[j], slope[j], res, erms = fit.elfit(npts_ap,weight,W,Z[:,j])                           
                            intercept[j] = intercept[j] + 2.0*math.log(DSunEarth)  # correct intercept to 1 AU
                        lnV0glob.append(intercept)
                        dayssinceepoch.append((obs['Date'] - pd.Timestamp(calepoch)).days)

                        numlangleys += 1
                        if verb: print(f'{obsdate} {ampm[iap]}')                       
                    else:
                        if verb: print(' Failed fit quality filter')
                else:
                    if verb: print(' Failed triplet cv filter') 
            else:
                if verb: print(' Insufficient triplets')

lnV0glob = np.vstack(lnV0glob)
lnV0coef0 = np.zeros(numchannels)
lnV0coef1 = np.zeros(numchannels)
erms = np.zeros(numchannels)

if cal.order[i870]==1 :
    weight = [1] * numlangleys
    for j in range(numchannels):
        lnV0coef0[j], lnV0coef1[j], res, erms[j] = fit.elfit(numlangleys,weight,dayssinceepoch,lnV0glob[:,j])
    lnV0coef0[iWV] = cal.lnV0coef0[iWV]
    lnV0coef1[iWV] = cal.lnV0coef1[iWV]
    lnV0coef = np.vstack([lnV0coef0, lnV0coef1]) 
    erms[iWV] = cal.erms[iWV]
else:
    lnV0conf = np.mean(lnV0glob, axis=0)
    lnV0coef.append(cal.lnV0coef0[iWV])
    
genfile = rootpath+'suncals/' + str(inst).zfill(2)+'/'+ calperiod_filename + caloutext
calfile = open(genfile, 'w')
calfile.write(f'# Calibration file generated from General Method between:\n'
       f'#{startdate.year:5d}{startdate.month:3d}{startdate.day:3d} and\n'
       f'#{enddate.year:5d}{enddate.month:3d}{enddate.day:3d}\n'
       f'#{inst:5d}        -- Instrument number\n'
       f'#{model:5d}        -- Model number     \n'
       f'#{numlangleychannels:5d}        -- Number of Langley wavelengths\n'
       f'#{numlangleys:5d}        -- Number of Langley intervals in period\n'
       f'#{numgencyc:5d}        -- Number of General Method cycles applied\n'
       f'#{calepoch.year:5d}{calepoch.month:3d}{calepoch.day:3d}  -- Calibration epoch\n')


calfile.write(f'# Wavel(nm) Order      {" ".join(f"LnV0({i})   " for i in range(cal.order[0] + 1))}     Erms\n')

for j in range(numchannels):
   calfile.write(f'{c.wavelength[model-1][j]:10.1f} {cal.order[j]:6d} ' +
           ' '.join(f'{lnV0coef[i, j]:13.5e}' for i in range(cal.order[j] + 1)) +
           f' {erms[j]:13.5e}\n')

calfile.close()
print(f'General calibration file written to: {genfile}')         

