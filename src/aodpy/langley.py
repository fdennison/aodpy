import datetime as dt
import pandas as pd
import numpy as np
import os as os
import math
import fileread_module as fr
from solpos_module import *
import calfit_module as fit
import atmos_module as atm

pi = math.pi
rad2deg = 180.0 / pi
deg2rad = pi /180.0

site = 'jb1'
rootpath = '/home/599/fd0474/AODcode/SampleData/'
startdate = dt.date(2015,6,2)  #
enddate = dt.date(2016,6,13)   #
verb=False

calepoch = dt.date(2015,1,1)
fitV0 = True

configfile = rootpath+'config/'+site+'.cfn'
config = fr.Read_Site_Configuration(configfile, startdate)
stationlat = config.attrs['Latitude'] #-12.6607
stationlon = config.attrs['Longitude'] #132.8931
model = config.CimelModel
inst = config.CimelNumber
numchannels = 10
numlangleychannels = 9

# Cut File
calperiod_filename = str(inst) + str(startdate.year % 100).zfill(2) + str(startdate.month).zfill(2) + str(enddate.month).zfill(2)
cut = pd.read_csv(rootpath + 'suncals/#' + str(inst) +'/#' + calperiod_filename + '.cut', skiprows=1, header=None, delimiter=r'\s+', 
                 names=['year','month','day','AmPm'])
cut['datetime'] = pd.to_datetime(cut[['year','month','day']])
cut = cut.set_index('datetime')

## Use reference channel for regression QA?(F:870nm)
#UseRefChannel4QA = True

#Max sd of fit for regression (normally 0.005)
MaxSdevFit = 0.005

o3 = fr.read_bom_ozone2011(rootpath+'ozone/'+site+'.o3')

i440 = 1 - 1
i670 = 2 - 1
i870 = 3 - 1
i1020 = 4 - 1
iWV = 10 - 1

datelist = pd.date_range(startdate,enddate,freq='d')  

numlangleys=0
dayssinceepoch=[]
lnV0glob=[]
ampm = ['Am','Pm']
for dii in datelist:
    obsdate = pd.to_datetime(dii).date()
    if verb: print('----'+obsdate.strftime("%Y/%m/%d")+'----')
    filedir = site+'/#'+str(inst)+'/'+str(obsdate.year)+'/'+str(obsdate.month).zfill(2)+'/'+str(obsdate.day).zfill(2)+'/'
    fileroot = obsdate.strftime("%y%m%d")
    
    # photometer data
    sunfile = rootpath+'agsdat/'+filedir+fileroot+'.sun'
    if os.path.isfile(sunfile):
        sundata = fr.read_triple_sun_record(sunfile, model)
        startlocalday = dt.datetime.combine(obsdate,dt.time(4,0,0)) - dt.timedelta(hours=config.attrs['TimeZone'])
        endlocalday = dt.datetime.combine(obsdate,dt.time(21,0,0)) - dt.timedelta(hours=config.attrs['TimeZone'])
        n1=len(sundata)
        sundata = sundata.loc[startlocalday:endlocalday]
        nobs = len(sundata)
    else:
        print(f'no obs on {obsdate}')
        continue
    
    # get ozone for the day
    if sum(o3['date']==obsdate): # daily ozone available
        ozonecolumn =  (o3.loc[o3['date']==obsdate].ozone / 1000.0).iloc[0]
    else: # use monthly mean
        if verb: print('use monthly mean o3')
        o3m = o3.loc[(pd.to_datetime(o3['date']).dt.month==obsdate.month) & 
                     (pd.to_datetime(o3['date']).dt.year==obsdate.year)].ozone
        ozonecolumn = sum(o3m)/len(o3m) / 1000.0

    # Pressure data
    #presfile = rootpath + 'agsdat/' + filedir + fileroot + '.hpa'  # date formats will be different
    presfile = rootpath + 'PyOut/' + filedir  + fileroot + '.hpa'
    p = fr.read_pressure_file(presfile)


    #print(f'{n1} to {nobs} in sundata')

    # Black record
    blkfile = rootpath+'agsdat/'+filedir+fileroot+'.blk'
    if os.path.isfile(blkfile):
        print(f'blackfile : {blkfile}')
        blksun = read_black_record(blkfile,model)
    else:
        blksun = [0] * numchannels

    
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
        
        obs = sundata.iloc[k]

        datetime = dt.datetime.combine(obs['Date'] , obs['Time'])
        #if verb: print(datetime)
        
        if len(p)==0:
            pk = defaultpressure
        else:
            pk = p.asof(pd.Timestamp(datetime)).Pressure
        pr.append(pk)
        
        for i in range(3):
        # Solar Position
            julday, timeday = julian(datetime+dt.timedelta(seconds=(30*i)))            
            sunpos = solar_position_almanac(julday,timeday)   
            DSunEarth=sunpos.DSunEarth
            solzen, solazi = satellite_position(julday, timeday, stationlat*deg2rad, stationlon*deg2rad, sun_elements)
            if i==1:
                if (solazi*rad2deg)<180.0:
                    numsuntriples[1] +=1     #morning
                else:
                    numsuntriples[2] +=1     #afternoon
            solzenapp[k,i] = apparent_zenith(solzen * rad2deg)

            ozoneOD = [x*ozonecolumn for x in atm.ozonecoef[model-1]]            

            # Signal 
            for n in range(numchannels):
                
                volt[k,i,n] = obs['Ch'+str(n+1)][i] - blksun[n]
                
                rayleighOD[n] = atm.rayleigh(atm.wavelength[model-1][n], pk)           

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
    
    if verb: print('Langley processing...')
    if verb: print(f'N Sun Triples: {numsuntriples}')    
    nstart = [0,0]
    nend = [0,0]
    for iap in [1,2]:
        nstart[iap-1] = numsuntriples[1] * (iap - 1)
        nend[iap-1] = nstart[iap-1] + numsuntriples[iap]
    #print(f'total:{numsuntriples[0]} triples, {numsuntriples[1]} am [{nstart[0]}:{nend[0]}], {numsuntriples[2]} pm [{nstart[1]}:{nend[1]}]')

    
    MinSunTriples = 10
    for iap in [0,1]: #  morning / afternoon
        meanpressure = np.mean(pr[nstart[iap]:nend[iap]])
        rayleighOD = [atm.rayleigh(wl, pk) for wl in atm.wavelength[model-1]]
        
        # cutList
        if dii in cut.index:
            if ampm[iap] in cut.loc[dii].AmPm:
                include=False
                print(f'cut {obsdate} {ampm[iap]}')
            else:
                include=True
        else:        
            include=True
        
        if include:
            if numsuntriples[iap]>=MinSunTriples :
                SpreadFlag, numOk, indOk = fit.CheckTripletCv(verb, numsuntriples[iap+1], nstart[iap],airmass[:,1,i870],tripletCv[:,i870])                
                
                if (numOk>=MinSunTriples) & SpreadFlag :
                    if verb: print(' Passed triplet cv langley filter')

                    n=numchannels
                    V = np.empty([numOk*3])
                    W = np.empty([numOk*3])
                    X = np.empty([numOk*3,n])
                    Y = np.empty([numOk*3,n])
                    Z = np.empty([numOk*3,n])
                    i=0 # Singlet index        
                    for k in range(numOk):              # Ok Triplet index
                        for j in range(3):                     # Triplet components
                            X[i,:]=airmass[indOk[k],j,:]
                            Y[i,:]=voltlog[indOk[k],j,:]
                            i=i+1
                    npts_ap=i

                    SpreadFlag,FitFlag, fitindex = fit.CheckFitQuality(verb,npts_ap,X,Y,i870,MaxSdevFit)    
                    if SpreadFlag & FitFlag: 
                        
                        weight = [1.0] * npts_ap
                        intercept = np.zeros(numchannels)  
                        slope = np.zeros(numchannels) 
                        for n in range(numlangleychannels):
                            intercept[n], slope[n], Residual, Erms, Delintercept, DelSlope = fit.boxfit(npts_ap,X[:,n],Y[:,n])
                            # correct intercept to 1 AU
                            intercept[n]=intercept[n] + 2.0*math.log(DSunEarth)
                        lnV0glob.append(intercept)
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
                        V  = Y[:,iWV] + X[:,iWV]*(rayleighOD[iWV]+aerosolOD[iWV])
                        weight = [1] * npts_ap
                        intercept[iWV], slope[iWV], res, erms = fit.elfit(npts_ap,weight,U,V)
                        intercept[iWV] = intercept[iWV] + 2*math.log(DSunEarth)
                        if slope[iWV]<0:
                            lnV0glob[numlangleys][iWV] = intercept[iWV]
                        else:
                            if verb: print('Water vapour below detection threshold')
                        dayssinceepoch.append((obs['Date'] - pd.Timestamp(calepoch)).days)    

                        numlangleys += 1
                        
                        print(f'{obsdate} {ampm[iap]}')
                            
                    else: # SpreadFlag & FitFlag:  
                        if verb: print(' Failed fit quality filter')
                else:  # (numOk>=MinSunTriples) & SpreadFlag :
                    if verb: print('Failed triplet cv langley filter')
            else: # numsuntriples[iap]>=MinSunTriples
                if verb: print(f'insufficient triplets [{numsuntriples[iap]}]')

print(f'{numlangleys} langleys')
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
    
    calfile = open(rootpath+'PyOut/' + site + '/' + calperiod_filename +'.lcl', 'w')
    calfile.write(f'# Calibration file generated from Langleys between:\n'
           f'#{startdate.year:5d}{startdate.month:3d}{startdate.day:3d} and\n'
           f'#{enddate.year:5d}{enddate.month:3d}{enddate.day:3d}\n'
           f'#{inst:5d}        -- Instrument number\n'
           f'#{model:5d}        -- Model number     \n'
           f'#{numlangleychannels:5d}        -- Number of Langley wavelengths\n'
           f'#{numlangleys:5d}        -- Number of Langley intervals in period\n'
           #f'#{numgeneralcycles:5d}        -- Number of General Method cycles applied\n'
           f'#{calepoch}  -- Calibration epoch\n')
    
    calfile.write(f'# Wavel(nm) Order      {" ".join(f"LnV0({I})   " for I in range(Iorder + 1))}     Erms\n')
    
    for N in range(numchannels):
       calfile.write(f'{atm.wavelength[model-1][N]:10.1f} {Iorder:6d} ' +
               ' '.join(f'{lnV0coef[I, N]:13.5e}' for I in range(Iorder + 1)) +
               f' {erms[N]:13.5e}\n')
    
    
    calfile.close()
