import datetime as dt
import numpy as np
import pandas as pd
import os
import math
import fileread_module as fr
import solpos_module as sol
import calfit_module as fit
import atmos_module as atm
import tomllib
import argparse
#import matplotlib.pyplot as plt        #problems getting matplotlib to work in my env

pi = math.pi
rad2deg = 180.0 / pi
deg2rad = pi /180.0

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
if args.verbose: print(langleyconf)    

# site config
configfile = rootpath+'config/'+site+'.cfn'
config = fr.read_site_configuration(configfile, startdate)
stationlat = config.attrs['Latitude'] 
stationlon = config.attrs['Longitude'] 
model = config.CimelModel
inst = config.CimelNumber
numchannels = len(atm.wavelength[model-1])
numlangleychannels = numchannels - 1

calperiod_filename = str(inst) + str(startdate.year % 100).zfill(2) + str(startdate.month).zfill(2) + str(enddate.month).zfill(2)

if calepoch=="":  calepoch = dt.date(startdate.year,1,1)

cutfile = rootpath + 'suncals/' + str(inst).zfill(2) +'/' + calperiod_filename + '.cut'
if os.path.isfile(cutfile):
    cut = fr.read_adj_file(cutfile, 'AmPm')
    usecut = True
else:
    usecut = False
    print('no cut file')

if clockfix:
    #clk = fr.read_adj_file(rootpath + 'suncals/' + str(inst).zfill(2) +'/' + calperiod_filename + '.clk', 'timecorr')
    clk = fr.read_adj_file(rootpath + 'clockdrift/' + site + '.clk', 'timecorr')
    
ozonefile = rootpath+ozonedir+site+'.o3'
if os.path.isfile(ozonefile):
    o3 = fr.read_bom_ozone2011(ozonefile)
    ozoneopt = True
else:
    ozoneopt = False
    print('using default ozone: 250 DU')
    
if langleytable:
    ltbfile = rootpath + 'suncals/' + str(inst).zfill(2) +'/' + calperiod_filename +'.ltb'
    os.makedirs(os.path.dirname(ltbfile), exist_ok=True)
    with open(ltbfile, 'w') as f:
        f.write('datetime CalD T P(hPa)  N '+' '.join([format(f'LV0{str(int(x)).zfill(4)}', "6s") for x in atm.wavelength[model-1]])+'\n')

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

numlangleys=0
dayssinceepoch=[]
lnV0glob=[]
ampm = ['Am','Pm']
datelist = pd.date_range(startdate,enddate,freq='d')  

for dii in datelist:
    obsdate = pd.to_datetime(dii).date()
    if args.verbose: print('----'+obsdate.strftime("%Y/%m/%d")+'----')
    filedir = site+'/#'+str(inst).zfill(2)+'/'+str(obsdate.year)+'/'+str(obsdate.month).zfill(2)+'/'+str(obsdate.day).zfill(2)+'/'
    fileroot = obsdate.strftime("%y%m%d")
    
    # photometer data
    sunfile = rootpath+'agsdat/'+filedir+fileroot+'.sun'
    if os.path.isfile(sunfile):
        sundata = fr.read_triple_sun_record(sunfile, model)
        if sundata.empty or (sum(sundata.Valid==0)==0):  # Fortran does not have valid data check here - without it it seems to produce spurious data on days where all data is invalid
            if args.verbose: print(f'skipping {obsdate}')
            continue
        startlocalday = dt.datetime.combine(obsdate,dt.time(4,0,0)) - dt.timedelta(hours=config.attrs['TimeZone'])
        endlocalday = dt.datetime.combine(obsdate,dt.time(21,0,0)) - dt.timedelta(hours=config.attrs['TimeZone'])
        n1=len(sundata)
        sundata = sundata.loc[startlocalday:endlocalday]
        nobs = len(sundata)
    else:
        if args.verbose: print(f'no obs on {obsdate}')
        continue
    
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
    presfile = rootpath + 'PyOut/' + filedir.replace('#','')  + fileroot + '.hpa'
    p = fr.read_pressure_file(args.verbose, presfile)

    # Black record
    blkfile = rootpath+'agsdat/'+filedir+fileroot+'.blk'
    if os.path.isfile(blkfile):
        if args.verbose: print(f'blackfile : {blkfile}')
        blksun = fr.read_black_record(blkfile,model)
    else:
        blksun = [0] * numchannels

    # Clock Fix
    if clockfix:  
        timecorr = clk.asof(dii).timecorr  # this is an integer in fortran version, hence use of round function
                                                            # there is no good reason to round this, so perhaps remove this in the future
                                                            # however it will substansially change output
        if np.isnan(timecorr):
            timecorr = round(clk.timecorr.iloc[0])   # prior to first entry in clk file
        else:
            timecorr = round(timecorr)
  
    volt = np.zeros([nobs,3,numchannels])
    voltlog = np.empty([nobs,3,numchannels])
    airmass  = np.empty([nobs,3,numchannels])
    airmassODrayleigh = np.empty([nobs,3,numchannels])
    airmassODaerosol  = np.empty([nobs,3,numchannels])
    airmassODozone  = np.empty([nobs,3,numchannels]) 
    solzenapp  = np.empty([nobs,3])
    numsuntriples = [nobs,0,0] # total/morning/afternoon
    pr=[]
    for k in range(nobs): 
        #print(f'[{j}]')
        rayleighOD = np.empty([numchannels])
        #lnV0 = np.empty([numchannels])
        
        obs = sundata.iloc[k]
        datetime = obs.name
        #datetime = dt.datetime.combine(obs['Date'] , obs['Time'])
        
        if clockfix: datetime = datetime - dt.timedelta(seconds=timecorr)
        #if args.verbose: print(datetime)
        
        if len(p)==0:
            pk = defaultpressure
        else:
            pk = p.asof(pd.Timestamp(datetime)).Pressure
            if math.isnan(pk):   # before first pressure measurement use daily mean
                pk = p.Pressure.mean()
        pr.append(pk)
        
        for i in range(3):
        # Solar Position
            julday, timeday = sol.julian(datetime+dt.timedelta(seconds=(30*i)))            
            sunpos = sol.solar_position_almanac(julday,timeday) 
            if (i==0) & (k==0): DSunEarth=sunpos.DSunEarth
            solzen, solazi = sol.satellite_position(julday, timeday, stationlat*deg2rad, stationlon*deg2rad, sol.sun_elements)
            if i==1:
                if (solazi*rad2deg)<180.0:
                    numsuntriples[1] +=1     #morning
                else:
                    numsuntriples[2] +=1     #afternoon
            solzenapp[k,i] = sol.apparent_zenith(solzen * rad2deg)

            ozoneOD = [x*ozonecolumn for x in atm.ozonecoef[model-1]]            

            # Signal 
            for n, wl in enumerate(atm.wavelength[model-1]):
                
                volt[k,i,n] = obs['ch'+str(int(wl))][i] - blksun[n]
                
                rayleighOD[n] = atm.rayleigh(wl, pk)           

                #Assume aod of 0.03 to calculate airmass. Refine later.
                airmassODrayleigh[k,i,n] = atm.getairmass(1,solzenapp[k,i]) * rayleighOD[n]
                airmassODaerosol[k,i,n] = atm.getairmass(2,solzenapp[k,i]) * 0.03
                airmassODozone[k,i,n] = atm.getairmass(3,solzenapp[k,i]) * ozoneOD[n]
                airmass[k,i,n] = ( airmassODrayleigh[k,i,n] + airmassODaerosol[k,i,n] + airmassODozone[k,i,n]) \
                                / (rayleighOD[n] + 0.03 + ozoneOD[n])       

                #print('volt = ',volt)
                if volt[k,i,n]>0:
                    voltlog[k,i,n] = math.log(volt[k,i,n]) 
                else:
                    voltlog[k,i,n]=-9.99

    tripletmean = np.mean(volt, axis=1)
    tripletCv = np.divide(np.std(volt, axis=1, ddof=1), tripletmean, out=np.ones_like(tripletmean)*99.99, where=tripletmean != 0)  * 100.0
    
    if args.verbose: print('Langley processing...')
    if args.verbose: print(f'N Sun Triples: {numsuntriples}')  
    nstart = [0,0]    
    nend = [0,0]
    for iap in [1,2]:
        nstart[iap-1] = numsuntriples[1] * (iap - 1)
        nend[iap-1] = nstart[iap-1] + numsuntriples[iap]
    #print(f'total:{numsuntriples[0]} triples, {numsuntriples[1]} am [{nstart[0]}:{nend[0]}], {numsuntriples[2]} pm [{nstart[1]}:{nend[1]}]')

    
    MinSunTriples = 10
    for iap in [0,1]: #  morning / afternoon        
        include=True
        # cutList
        if usecut:
            if dii in cut.index:
                if ampm[iap] in cut.loc[dii].AmPm:
                    include=False
                    if args.verbose: print(f'cut {obsdate} {ampm[iap]}')
        
        if include:
            if numsuntriples[iap]>=MinSunTriples :
                SpreadFlag, numOk, indOk = fit.CheckTripletCv(args.verbose, numsuntriples[iap+1], nstart[iap],airmass[:,1,i870],tripletCv[:,i870])                
                
                if (numOk>=MinSunTriples) & SpreadFlag :
                    if args.verbose: print(' Passed triplet cv langley filter')

                    n=numchannels
                    V = np.empty([numOk*3])
                    W = np.empty([numOk*3])
                    X = np.empty([numOk*3,n])
                    Y = np.empty([numOk*3,n])
                    Z = np.empty([numOk*3,n])
                    i=0 # Singlet index        
                    for k in range(numOk):              # Ok Triplet index
                        for j in range(3):                     # Triplet components
                            X[i,:] = airmass[indOk[k],j,:]
                            Y[i,:] = voltlog[indOk[k],j,:]
                            i = i+1
                    npts_ap = i

                    SpreadFlag,FitFlag, fitindex, npts_ap = fit.CheckFitQuality(args.verbose,npts_ap,X,Y,i870,MaxSdevFit)    
                    if SpreadFlag & FitFlag: 
                        
                        weight = [1.0] * npts_ap
                        intercept = np.zeros(numchannels)  
                        slope = np.zeros(numchannels) 
                        for n in range(numlangleychannels):
                            intercept[n], slope[n], Residual, Erms, Delintercept, DelSlope = fit.boxfit(npts_ap,X[:,n],Y[:,n])
                            # correct intercept to 1 AU
                            intercept[n]=intercept[n] + 2.0*math.log(DSunEarth)
                        lnV0glob.append(intercept)

                        meanpressure = np.mean(pr[nstart[iap]:nend[iap]])
                        rayleighOD = [atm.rayleigh(wl, meanpressure) for wl in atm.wavelength[model-1]]
                        aerosolOD = -slope - rayleighOD - ozoneOD
                            
                        angstrom = np.zeros([numlangleychannels,numlangleychannels])
                        for N in range(numlangleychannels-1):
                            for M in range(N+1, numlangleychannels):
                                if aerosolOD[N] >= 0 and aerosolOD[M] >= 0:
                                    angstrom[N, M] = -math.log(aerosolOD[N] / aerosolOD[M]) / \
                                            math.log(atm.wavelength[model-1][N] / atm.wavelength[model-1][M])
                                else:
                                    angstrom[N, M] = -9.999

                        aerosolOD[iWV] = aerosolOD[i870] * (atm.wavelength[model-1][iWV]/atm.wavelength[model-1][i870])**(-angstrom[i670,i870])
                        Bcoef=0.574
                        U = X[:,iWV]**Bcoef
                        V = Y[:,iWV] + X[:,iWV]*(rayleighOD[iWV]+aerosolOD[iWV])
                        weight = [1] * npts_ap
                        intercept[iWV], slope[iWV], res, erms = fit.elfit(npts_ap,weight,U,V)
                        intercept[iWV] = intercept[iWV] + 2*math.log(DSunEarth)
                        if slope[iWV]<0:
                            lnV0glob[numlangleys][iWV] = intercept[iWV]
                        else:
                            if args.verbose: print('Water vapour below detection threshold')
                        cald = (obs['Date'] - pd.Timestamp(calepoch)).days
                        dayssinceepoch.append(cald)    

                        numlangleys += 1
                        if args.verbose: print(f'{obsdate} {ampm[iap]}  {npts_ap} {meanpressure:6.1f}')
                        if langleytable:
                            with open(ltbfile, 'a') as f:
                                f.write( f'{obsdate} {cald} {ampm[iap]} {meanpressure:6.1f} {npts_ap} ' + " ".join([format(x, "6.4f") for x in intercept]) +'\n')
                    else: # SpreadFlag & FitFlag:  
                        if args.verbose: print(' Failed fit quality filter')
                else:  # (numOk>=MinSunTriples) & SpreadFlag :
                    if args.verbose: print('Failed triplet cv langley filter')
            else: # numsuntriples[iap]>=MinSunTriples
                if args.verbose: print(f'insufficient triplets [{numsuntriples[iap]}]')

print(f'{numlangleys} langleys in this time period')


lnV0glob = np.vstack(lnV0glob)

lnV0coef0 = np.zeros(numchannels)
lnV0coef1 = np.zeros(numchannels)
erms = np.zeros(numchannels)

if numlangleys>0:
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
           f'#{numlangleychannels:5d}        -- Number of Langley wavelengths\n'
           f'#{numlangleys:5d}        -- Number of Langley intervals in period\n'
           '#    0        -- Number of General Method cycles applied\n'
           f'#{calepoch.year:5d}{calepoch.month:3d}{calepoch.day:3d}  -- Calibration epoch\n')
    
    
    calfile.write(f'# Wavel(nm) Order      {" ".join(f"LnV0({I})   " for I in range(Iorder + 1))}     Erms\n')
    
    for N in range(numchannels):
       calfile.write(f'{atm.wavelength[model-1][N]:10.1f} {Iorder:6d} ' +
               ' '.join(f'{lnV0coef[I, N]:13.5e}' for I in range(Iorder + 1)) +
               f' {erms[N]:13.5e}\n')
        
    calfile.close()

if makeplot:
    fig, ax = plt.subplots(figsize=(12,9))
    for n in range(numchannels):
        p = ax.plot(dayssinceepoch,lnV0glob[:,n],'.')
        linfit = np.polyval([lnV0coef1[n], lnV0coef0[n]], dayssinceepoch)
        col = p[-1].get_color()
        ax.plot(dayssinceepoch, linfit, color=col, label=atm.wavelength[model-1][n])
    ax.set_ylabel(r'lnV$_0$')
    ax.set_xlabel('days since epoch')
    ax.set_title(calperiod_filename+'.lcl')
    plt.legend()
    figfile = rootpath+'suncals/'+str(inst).zfill(2)+'/'+calperiod_filename+'.lcl.ps'
    plt.savefig(figfile,bbox_inches='tight')
    print(f'Figure written to: {figfile}, check for outliers')

if langleytable:
    print(f'Langley table written to: {ltbfile}')
print(f'Langley calibration file written to: {langleyfile}')    
