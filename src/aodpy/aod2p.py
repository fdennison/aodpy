import datetime as dt
import pandas as pd
import numpy as np
import os as os
import math
import fileread_module as fr
from solpos_module import *
from atmos_module import *

pi = math.pi
rad2deg = 180.0 / pi
deg2rad = pi /180.0
avogadro = 6.022529e+23  # mol^{-1}
M_dryair = 2.8964e-02  # kg/mol

site = 'jb1'
rootpath = '/home/599/fd0474/AODcode/SampleData/'
startdate = dt.date(2016,1,1)
enddate = dt.date(2016,1,30)

calfile_in = '#12150606.500'


cal = fr.read_cal(rootpath+'suncals/'+calfile_in[0:3]+'/'+calfile_in)

o3 = fr.read_bom_ozone2011(rootpath+'ozone/'+site+'.o3')

configfile = rootpath+'config/'+site+'.cfn'
config = fr.Read_Site_Configuration(configfile, startdate)

I440 = 1 - 1
I670 = 2 - 1
I870 = 3 - 1
I1020 = 4 - 1
IWV = 10 - 1

stationlat = config.attrs['Latitude'] #-12.6607
stationlon = config.attrs['Longitude'] #132.8931
model = config.CimelModel
inst = config.CimelNumber
numchannels = 10
stacksize=3

window = 2 # minutes
CVmax = 0.01
defaultpressure = 1013.0 # hPa

datelist = pd.date_range(startdate,enddate,freq='d')  

verb=False

for dii in datelist:
    obsdate = pd.to_datetime(dii).date()
    print('----'+obsdate.strftime("%Y/%m/%d")+'----')
    filedir = site+'/#'+str(inst)+'/'+str(obsdate.year)+'/'+str(obsdate.month).zfill(2)+'/'+str(obsdate.day).zfill(2)+'/'
    fileroot = obsdate.strftime("%y%m%d")

    # get ozone for the day
    if sum(o3['date']==obsdate): # daily ozone available
        ozonecolumn =  (o3.loc[o3['date']==obsdate].ozone / 1000.0).iloc[0]
    else: # use monthly mean
        print('use monthly mean o3')
        o3m = o3.loc[(pd.to_datetime(o3['date']).dt.month==obsdate.month) & 
                     (pd.to_datetime(o3['date']).dt.year==obsdate.year)].ozone
        ozonecolumn = sum(o3m)/len(o3m) / 1000.0

    # Pressure data
    #presfile = rootpath + 'agsdat/' + filedir + fileroot + '.hpa'  # date formats will be different
    presfile = rootpath + '/PyOut/' + config.attrs['Id'] + '/' + fileroot + '.hpa'
    p = fr.read_pressure_file(presfile)
    # prescolumnnames=['Date','Time','P_UnTCorrect','Tmean','DelP','Pressure']
    # p = pd.read_csv(presfile, skiprows=9, header=None, delimiter=r'\s+', names=prescolumnnames, index_col=False)
    # p['DateTime'] = pd.to_datetime(p['Date'] + ' ' + p['Time'])
    # p = p.set_index('DateTime')
    # if np.any(np.diff(p.index.to_list())<dt.timedelta(hours=0)): 
    #     print('time disordered in pressure file')
    #     p.sort_index(inplace=True)       # jb1 160105 has a dicontinuity, indicate problem?
    # if len(p)==0:
    #     print('problem with pressure file')


    # photometer data
    sundata = fr.read_single_sun_record(rootpath+'agsdat/'+filedir+fileroot+'.lsu', model)
    nobs = len(sundata)
    stackloop = nobs - stacksize + 1

    # Black record
    blkfile = rootpath+'agsdat/'+filedir+fileroot+'.blk'
    if os.path.isfile(blkfile):
        print(f'blackfile : {blkfile}')
        blksun = read_black_record(blkfile,model)
    else:
        blksun = [[0]] * numchannels

    filt1 = set()
    aod = np.zeros([nobs, numchannels])
    for j in range(stackloop): 
        #print(f'[{j}]')
        rayleighOD = np.empty([numchannels])
        #lnV0 = np.empty([numchannels])
        volt = np.empty([stacksize, numchannels])
        voltlog = np.empty([stacksize, numchannels])
        solzen = np.empty([stacksize])
        solzenapp = np.empty([stacksize])
        airmass  = np.empty([stacksize, numchannels])


        stack = sundata.iloc[j:(j+stacksize)]
        #print(stack)

        datetime = [dt.datetime.combine(d,t) for d,t in zip(stack['Date'] , stack['Time'])]
        if verb: print(datetime[0])


        #   1st pass Filters
        #--------------------------
        # stack of obs does not fit in specified time window
        if ((datetime[stacksize-1] - datetime[0]).total_seconds())/60.0 > window:         
            if verb: print('stack of obs does not fit in the {} minute window : {} {} {}'.format(
                           window,*datetime))
            continue

        Valid=True  
        for i in range(stacksize):
            # temperature/signal to high or low - assessed in read function
            if not stack.iloc[i].Valid==0:
                if verb: print(f'Stack not valid : {stack.iloc[i].Valid}')
                Valid=False
                break

            #Exclude if obs obtained before 0300 LST
            julday, timeday = julian(datetime[i])
            timeUT = timeday - 0.5   # why 0.5?
            timeAEST = timeUT*24 + 10
            LSTcut = dt.datetime.combine(obsdate,dt.time(3,0,0))
            UTcut = LSTcut - dt.timedelta(hours=config.attrs['TimeZone'])
            if datetime[i] < UTcut:
                if verb: print('before LST cutoff')
                Valid=False
                break
        if not Valid:
            continue

        # Covariance of stack is above specified limit   
        for n in range(numchannels):
            for i in range(stacksize):
                # Calculate signal
                volt[i,n] = stack.iloc[i]['Ch'+str(n+1)] - blksun[n]        

        if verb: print(volt[0,:])
        CVvolt = np.std(volt[:,I870]) / np.mean(volt[:,I870])
        if CVvolt > CVmax:
            if verb: print(f'High Covariance {CVvolt:3.5f}')  
            continue



        if Valid: 
            if verb: print('add : {} {} {}'.format(*list(range(j,(j+stacksize)))))
            [filt1.add(x) for x in range(j,(j+stacksize))]
        #    [print(d) for d in datetime]
            if verb: print('GOOD STACK')
        #----------------------------


        if len(p)==0:
            pr = defaultpressure
        else:
            pr = p.asof(pd.Timestamp(datetime[0])).Pressure
        #print(pr)

        # Cal File
        nday = (obsdate - cal.attrs['epoch']).days
        numlangch = cal.attrs['numlangleychannels']
        lnV0_1AU=[0]*numchannels
        for n in range(numlangch):
            lnV0_1AU[n] = cal.iloc[n].lnV0coef0 +  cal.iloc[n].lnV0coef1 * nday   #
        n=IWV
        lnV0_1AU[n] = cal.iloc[n].lnV0coef0 +  cal.iloc[n].lnV0coef1 * nday # point of doing this channel seperately?

        # Solar Position
        for i in range(stacksize):
            dti = datetime[i]
            julday, timeday = julian(dti)            
            sunpos = solar_position_almanac(julday,timeday)   
            DSunEarth=sunpos.DSunEarth
            solzen[i], solazi = satellite_position(julday, timeday, stationlat*deg2rad, stationlon*deg2rad, sun_elements)
            solzenapp[i] = apparent_zenith(solzen[i] * rad2deg)

        # Signal
        ozoneOD = [x*ozonecolumn for x in ozonecoef[model-1]]  
        extraOD = [0]*numchannels       
        for n in range(numchannels):
            rayleighOD[n] = rayleigh(wavelength[model-1][n], pr)
            for i in range(stacksize):            

                #print('rayleighOD = ',rayleighOD)
                #Assume aod of 0.03 to calculate airmass. Refine later.
                airmass[i,n] = (getairmass(1,solzenapp[i]) * rayleighOD[n] + \
                                getairmass(2,solzenapp[i]) * 0.03 + \
                                getairmass(3,solzenapp[i]) * ozoneOD[n] ) \
                                / (rayleighOD[n] + 0.03 + ozoneOD[n])       

                #print('volt = ',volt)
                if volt[i,n]>0:
                    voltlog[i,n] = math.log(volt[i,n]) 
                else:
                    voltlog[i,n]=-9.99        

                lnV0 = lnV0_1AU[n] - 2.0 * math.log(DSunEarth)
                #print('lnV0 = ',lnV0)
                #print('voltlog = ',voltlog)
                #tran
# 22:24:35
# for i=0,n=2             9.6515  9.1456        2.7108        0.0149          0
                aod[j+i,n] = (lnV0-voltlog[i,n]) / airmass[i,n] - rayleighOD[n] - ozoneOD[n] - extraOD[n]
                #
        if verb: 
            print('SolZen = ',solzen[0] * rad2deg)
            print('rayl = ',rayleighOD[I870])
            print('airmass = ',airmass[0,I870])
            print('voltlog = ', voltlog[0,I870])
            print('AOD 870 = ',aod[j,I870])

    idx1 = sorted(list(filt1))
    if len(idx1)==0:
        print('no valid obs on this day')
    else:    
        if verb: [print(f'{d} {t} {a:5.4f}') for d,t,a in zip(sundata.iloc[idx1]['Date'],
                                                sundata.iloc[idx1]['Time'],
                                                aod[idx1,I870])]

        dailyaod = np.mean(aod[idx1,I870])
        stddev = np.std(aod[idx1,I870])
        nobs1 = len(idx1)
        print(f'daily mean AOD 870 = {dailyaod:5.4f}    n = {nobs1}')

        # 2nd pass filter
        for k, j in enumerate(idx1[:-(stacksize-1)]):
            # stack should be within 3sd of daily mean
            stackidx = idx1[k:(k+stacksize)]
            if abs(np.mean(aod[stackidx,I870]) - dailyaod) > 3*stddev:
                print('outside 3sd')
                [print(f'{d} {t} {a:5.4f}') for d,t,a in zip(sundata.iloc[stackidx]['Date'],
                                                sundata.iloc[stackidx]['Time'],
                                                aod[stackidx,I870])]


