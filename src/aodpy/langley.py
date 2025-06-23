import datetime as dt
import numpy as np
import pandas as pd
import os
import math
import fileread_module as fr
import solpos_module as sol
import calfit_module as fit
import atmos_module as atm
import constants_module as c 
import tomllib
import argparse
import matplotlib.pyplot as plt        

parser=argparse.ArgumentParser(description="Apply langley method to sun photometer triplet measuments (.sun files)")
parser.add_argument('-t', '--test', dest='configtoml', action='store_const', 
                    const='testconfig.toml', default='config.toml',
                    help='run on sample data and check output')
parser.add_argument('-v', '--verbose',action='store_true')
args=parser.parse_args()

# langley config options
with open(args.configtoml, "rb") as f:
    langleyconf = tomllib.load(f)
site = langleyconf['site']
rootpath = langleyconf['rootpath']
startdate = langleyconf['startdate']
enddate = langleyconf['enddate']
ozonedir = langleyconf['ozonedir']
langleytable = langleyconf['cal']['maketable']
makeplot = langleyconf['cal']['makeplot']
calepoch = langleyconf['cal']['calepoch']
fitV0 = langleyconf['cal']['fitV0']
clockfix = langleyconf['cal']['clockfix']
MaxSdevFit = langleyconf['cal']['MaxSdevFit']
print(f'Run langley for {site} from {startdate} to {enddate}')
verb = args.verbose
if verb: print(langleyconf)    

# site config
configfile = rootpath+'config/'+site+'.cfn'
config = fr.read_site_configuration(configfile, startdate)
stationlat = config.attrs['Latitude'] 
stationlon = config.attrs['Longitude'] 
model = config.CimelModel
inst = config.CimelNumber

numchannels = len(c.wavelength[model-1])

calperiod_filename = str(inst) + str(startdate.year % 100).zfill(2) + str(startdate.month).zfill(2) + str(enddate.month).zfill(2)
if calepoch=="":  calepoch = dt.date(startdate.year,1,1)

# cut file
cutfile = rootpath + 'suncals/' + str(inst).zfill(2) +'/' + calperiod_filename + '.cut'
if os.path.isfile(cutfile):
    cut = fr.read_adj_file(cutfile, 'AmPm')
    usecut = True
else:
    usecut = False
    print('no cut file')

# clock fix file
if clockfix:
    clk = fr.read_adj_file(rootpath + 'clockdrift/' + site + '.clk', 'timecorr')

# ozone file
ozonefile = rootpath+ozonedir+site+'.o3'
if os.path.isfile(ozonefile):
    o3 = fr.read_ozone(ozonefile)
    ozoneopt=True
else:
    ozoneopt=False
    print('using default ozone: 250 DU')

if langleytable:
    ltbfile = rootpath + 'suncals/' + str(inst).zfill(2) +'/' + calperiod_filename +'.ltb'
    os.makedirs(os.path.dirname(ltbfile), exist_ok=True)
    with open(ltbfile, 'w') as f:
        f.write('datetime CalD T P(hPa)  N '+' '.join([format(f'LV0{str(int(x)).zfill(4)}', "6s") for x in c.wavelength[model-1]])+'\n')

i500 = c.wl_index(500,c.wavelength[model-1])   
i670 = c.wl_index(670,c.wavelength[model-1])   
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

    # Clock Fix
    if clockfix:  
        timecorr = clk.asof(dii).timecorr  
        if np.isnan(timecorr):
            timecorr = clk.timecorr.iloc[0]   # if prior to first entry in clk file, use first entry
  
    volt = np.zeros([nobs,numchannels,3])
    voltlog = np.empty([nobs,numchannels,3])
    airmass  = np.empty([nobs,numchannels,3])
    solzenapp  = np.empty([nobs,3])
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
                airmassODrayleigh = atm.getairmass(1,solzenapp[i,k]) * rayleighOD
                airmassODaerosol = atm.getairmass(2,solzenapp[i,k]) * 0.03
                airmassODozone = atm.getairmass(3,solzenapp[i,k]) * ozoneOD[j]
                airmass[i,j,k] = ( airmassODrayleigh + airmassODaerosol + airmassODozone) \
                                / (rayleighOD + 0.03 + ozoneOD[j])       

                if volt[i,j,k]>0:
                    voltlog[i,j,k] = math.log(volt[i,j,k]) 
                else:
                    voltlog[i,j,k]=-9.99

    tripletmean = np.mean(volt, axis=2)
    tripletCv = np.divide(np.std(volt, axis=2, ddof=1), tripletmean, out=np.ones_like(tripletmean)*99.99, where=tripletmean != 0)  * 100.0
    
    if verb: print('Langley processing...')
    if verb: print(f'N Sun Triples: {numsuntriples}')  
    nstart = [0,0]    
    nend = [0,0]
    for iap in [1,2]:
        nstart[iap-1] = numsuntriples[1] * (iap - 1)
        nend[iap-1] = nstart[iap-1] + numsuntriples[iap]

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
                SpreadFlag, numOk, indOk = fit.CheckTripletCv(verb, numsuntriples[iap+1], nstart[iap],airmass[:,i870,1],tripletCv[:,i870])                
                
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
                            ii=ii+1
                    npts_ap=ii

                    SpreadFlag,FitFlag, fitindex, npts_ap = fit.CheckFitQuality(verb,npts_ap,X,Y,i870,MaxSdevFit)    
                    if SpreadFlag & FitFlag: 
                        
                        weight = [1.0] * npts_ap
                        intercept = np.zeros(numchannels)  
                        slope = np.zeros(numchannels) 
                        for j in range(numchannels):
                            intercept[j], slope[j], Residual, Erms, Delintercept, DelSlope = fit.boxfit(npts_ap,X[:,j],Y[:,j])
                            intercept[j] = intercept[j] + 2.0*math.log(DSunEarth)  # correct intercept to 1 AU
                        lnV0glob.append(intercept)

                        meanpressure = np.mean(pr[nstart[iap]:nend[iap]])
                        rayleighOD = [atm.rayleigh(wl, meanpressure) for wl in c.wavelength[model-1]]
                        aerosolOD = -slope - rayleighOD - ozoneOD
                            
                        angstrom = np.zeros([numchannels,numchannels])
                        for n in range(numchannels):
                            for m in range(numchannels):
                                if aerosolOD[n] >= 0 and aerosolOD[m] >= 0 and (n!=m):
                                    angstrom[n, m] = -math.log(aerosolOD[n] / aerosolOD[m]) / \
                                            math.log(c.wavelength[model-1][n] / c.wavelength[model-1][m])
                                else:
                                    angstrom[n, m] = -9.999

                        aerosolOD[iWV] = aerosolOD[i870] * (c.wavelength[model-1][iWV]/c.wavelength[model-1][i870])**(-angstrom[i670,i870])
                        Bcoef=0.574
                        U = X[:,iWV]**Bcoef
                        V  = Y[:,iWV] + X[:,iWV]*(rayleighOD[iWV]+aerosolOD[iWV])
                        weight = [1] * npts_ap
                        intercept[iWV], slope[iWV], res, erms = fit.elfit(npts_ap,weight,U,V)
                        intercept[iWV] = intercept[iWV] + 2.0*math.log(DSunEarth)
                        if slope[iWV]<0:
                            lnV0glob[numlangleys][iWV] = intercept[iWV]
                        else:
                            if verb: print('Water vapour below detection threshold')
                        cal_d = (obs['Date'] - pd.Timestamp(calepoch)).days        
                        dayssinceepoch.append(cal_d)    
                        numlangleys += 1

                        if langleytable:
                            with open(ltbfile, 'a') as f:
                                f.write( f'{obsdate} {cal_d} {ampm[iap]} {meanpressure:6.1f} {npts_ap} ' + " ".join([format(x, "6.4f") for x in intercept]) +'\n')    
                    else: # SpreadFlag & FitFlag:  
                        if verb: print(' Failed fit quality filter')                        
                else:  # (numOk>=MinSunTriples) & SpreadFlag :
                    if verb: print('Failed triplet cv langley filter')                    
            else: # numsuntriples[iap]>=MinSunTriples
                if verb: print(f'insufficient triplets [{numsuntriples[iap]}]')

print(f'{numlangleys} langleys in this time period')

if numlangleys>0:

    lnV0glob = np.vstack(lnV0glob)
    lnV0coef0 = np.zeros(numchannels)
    lnV0coef1 = np.zeros(numchannels)
    erms = np.zeros(numchannels)

    if numlangleys>=2 & fitV0 :
        Iorder = 1
        weight = [1] * numlangleys
        for n in range(numchannels):
            lnV0coef0[n], lnV0coef1[n], res, erms[n] = fit.elfit(numlangleys,weight,dayssinceepoch,lnV0glob[:,n])
        lnV0coef = np.vstack([lnV0coef0, lnV0coef1]) 
    else:
        Iorder = 0
        lnV0conf = np.mean(lnV0glob, axis=0)
    
    langleyfile = rootpath+'suncals/'+str(inst).zfill(2)+'/'+calperiod_filename+'.lcl'
    calfile = open(langleyfile, 'w')
    calfile.write(f'# Calibration file generated from Langleys between:\n'
           f'#{startdate.year:5d}{startdate.month:3d}{startdate.day:3d} and\n'
           f'#{enddate.year:5d}{enddate.month:3d}{enddate.day:3d}\n'
           f'#{inst:5d}        -- Instrument number\n'
           f'#{model:5d}        -- Model number     \n'
           f'#{(numchannels-1):5d}        -- Number of Langley wavelengths\n'
           f'#{numlangleys:5d}        -- Number of Langley intervals in period\n'
           '#    0            -- Number of General Method cycles applied\n'
           f'#{calepoch.year:5d}{calepoch.month:3d}{calepoch.day:3d}  -- Calibration epoch\n')
      
    calfile.write(f'# Wavel(nm) Order      {" ".join(f"LnV0({I})   " for I in range(Iorder + 1))}     Erms\n')
    
    for N in range(numchannels):
       calfile.write(f'{c.wavelength[model-1][N]:10.1f} {Iorder:6d} ' +
               ' '.join(f'{lnV0coef[I, N]:13.5e}' for I in range(Iorder + 1)) +
               f' {erms[N]:13.5e}\n')        
    calfile.close()

    print(f'Langley calibration file written to: {langleyfile}')    

    if makeplot:
        fig, ax = plt.subplots()
        for n in [i500, i670, i870]: #range(numchannels):
            p = ax.plot(dayssinceepoch,lnV0glob[:,n],'.')
            if fitV0:
                linfit = np.polyval([lnV0coef1[n], lnV0coef0[n]], dayssinceepoch)
                col = p[-1].get_color()
                ax.plot(dayssinceepoch, linfit, color=col, label=c.wavelength[model-1][n])
        ax.set_ylabel(r'lnV$_0$')
        ax.set_xlabel('days since epoch')
        ax.set_title(calperiod_filename+'.lcl')
        plt.legend()
        figfile = rootpath+'suncals/'+str(inst).zfill(2)+'/'+calperiod_filename+'.lcl.ps'
        plt.savefig(figfile,bbox_inches='tight')
        print(f'Figure written to: {figfile}, check for outliers')

    if langleytable:
        print(f'Langley table written to: {ltbfile}')
