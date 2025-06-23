import math
import datetime as dt
from cimel_module import *
from constants_module import *
import os
import time_module as tm
import navigation_module as nav
import solpos_module as sol
import statistics_module as st
import numpy as np

MaxLangleys = 2048
MaxPolyOrder = 8

Dbug = [False] * 4
LangFlag = [False] * 2
NewStation = False
MakeCal = False
MakeTab = False
ApplyClockFix = False
ClockFixFileOK = False
FitLineV0 = False
Include = False

DateObs = [None] * 3
CalEpoch = [None] * 3

CurrentTimeValues = [None] * 8
RunDateTime = [None] * 8

FormatString = [None] * 10

PData = [None] * MaxPressurePoints

AlmData = [None] * NumAerosolChannels
AlnData = [None] * NumAerosolChannels
PPData = [None] * NumAerosolChannels
CutList = [None] * MaxCutPeriods
JulDay = [None] * 3
YMDHMS = [[None] * 3 for _ in range(7)]
IWavelength = [None] * NumChannels


SolarZenithDeg = np.full((MaxTriplets, 3), None)
SolarZenithApp = np.full((MaxTriplets, 3), None)
SolarAzimuthDeg = np.full((MaxTriplets, 3), None)
Temperature = np.full(MaxTriplets, None)
Pressure = np.full(MaxTriplets, None)
DSunEarth = np.full(MaxTriplets, None)
MeanDetectorTemp = np.full(2, None)
PresMean = np.full(2, None)
AirMassODRayleigh = np.full((MaxTriplets, 3, NumChannels), None)
AirMassODAerosol = np.full((MaxTriplets, 3, NumChannels), None)
AirMassODOzone = np.full((MaxTriplets, 3, NumChannels), None)
AirMass = np.full((MaxTriplets, 3, NumChannels), None)
Signal = np.full((MaxTriplets, 3, NumChannels), None)
LnSignal = np.full((MaxTriplets, 3, NumChannels), None)
Triplet = np.full(3, None)
TripletMean = np.full((MaxTriplets, NumChannels), None)
TripletSdev = np.full((MaxTriplets, NumChannels), None)
TripletCv = np.full((MaxTriplets, NumChannels), None)
TripletLog = np.full((MaxTriplets, NumChannels), None)
BlackSun = np.full(NumChannels, None)
BlackSkyA = np.full(NumSkyA, None)
BlackSkyK = np.full(NumSkyK, None)
RayleighOD = np.full(NumChannels, None)
OzoneOD = np.full(NumChannels, None)
AerosolOD = np.full((2, NumChannels), None)
LnV0Glob = np.full((MaxLangleys, NumChannels), None)
AodGlob = np.full((MaxLangleys, NumChannels), None)
LnV0Mean = np.full(NumChannels, None)
LnV0Sdev = np.full(NumChannels, None)
AodMean = np.full(NumChannels, None)
AodSdev = np.full(NumChannels, None)
DaysSinceEpoch = np.full(MaxLangleys, None)
Weight = np.full(MaxLangleys, None)
Residual = np.full(NumChannels, None)
Erms = np.full(NumChannels, None)
LnV0Coef = np.full((MaxPolyOrder + 1, NumChannels), None)
ElnV0Coef = np.full((MaxPolyOrder + 1, NumChannels), None)
Angstrom = np.full((2, NumChannels, NumChannels), None)

FileOK = False
There = False
SpreadFlag = False
FitFlag = False
Iap = None
MaxRecLength = None
NumLangleys = None
NumWVLangleys = None
NumGeneralCycles = None
MinSunTriples = 10

Nstart = [None] * 2
Nend = [None] * 2
NumOK = [None] * 2
NumPoints = [None] * 2
IndOK = [[None] * 2 for _ in range(MaxTriplets)]
NumSunTriples = [None] * 3
NumSunSingles = [None] * 3

MaxSdevFit = None  # 0.005 canonical
X = [[None] * NumChannels for _ in range(MaxSinglets)]
Y = [[None] * NumChannels for _ in range(MaxSinglets)]
Intercept = [[None] * NumChannels for _ in range(2)]
DelIntercept = [[None] * NumChannels for _ in range(2)]
Slope = [[None] * NumChannels for _ in range(2)]
DelSlope = [[None] * NumChannels for _ in range(2)]
AmPm = ['Am', 'Pm']

H2O_k = 0.9156e-2
H2O_e = 0.5918

U = [None] * MaxSinglets
V = [None] * MaxSinglets
OD_WvapMean = [[None] * NumChannels for _ in range(2)]
WvapMean = [None] * 2

CurrentMonth = None
Ozone = [None] * 32

Slash = '/'
nav.nav_setup()
CurrentDateTime = dt.datetime.now(dt.timezone.utc).astimezone()
#CurrentDateTime.utcoffset() #gives timedelta in seconds
CurrentDate, CurrentTime, TimeZone = CurrentDateTime.strftime('%Y%m%d,%H%M%S,%:z').split(',')

#Menu
Dbug=[False,False,False,False]
MakeTab=True
MakeCal=True
FitLineV0=True
ApplyClockFix=False
MaxSdevFit=0.005
Extension='sun'
DataRoot='/home/599/fd0474/AODcode/SampleData'
DataSuffix='/agsdat'
ResultPath='/suncals'
ConfigRoot='/config'
SiteCode='jb1'
YearMin=2016
MonthMin=1
DayMin=1
YearMax=2016
MonthMax=1 #6
DayMax=2 #13
PressureOption=1
DefaultSurfacePressure=1013.0
BoMPresFile='mba-19980101-19991231.sp'
OzoneOption=0 #1
DefaultOzone=0.25
OzoneFile='jb1.o3'
CalEpoch=dt.date(2015, 1, 1)

ConfigPath = DataRoot+ConfigRoot
print('config path=', ConfigPath)
ConfigFile = ConfigPath+'/'+SiteCode+'.cfn'
print(' Config File=', ConfigFile)
ObsDate = dt.date(YearMin, MonthMin, DayMin)
Config = Read_Site_Configuration(Dbug[0], ConfigFile, ObsDate)
#print(Config)
Instrument = Config.CimelNumber
Model = Config.CimelModel
StationLatDeg = Config.attrs['Latitude']
StationLonDeg = Config.attrs['Longitude']
StationLatRad = StationLatDeg * RADIANS_PER_DEGREE
StationLonRad = StationLonDeg * RADIANS_PER_DEGREE
#TZ_hms = Hour2HMS(Config.TimeZone)
#TimeZone_DT = Date_Time(0, 0, 0, -TZ_hms.Hour, -TZ_hms.Minute, -TZ_hms.Second)
TimeZone_DT = tm.timedelta(hours=-Config.attrs['TimeZone'])
OzoneFile = DataRoot + Slash + 'ozone' + Slash + OzoneFile
BoMPresFile = DataRoot + Slash + 'pressure' + Slash + BoMPresFile

ResultPath = DataRoot + ResultPath + Slash + '#' + str(Instrument).zfill(2) + Slash
SkyRoot = DataRoot + Slash + 'agsout' + Slash + SiteCode + Slash + 'sky' + Slash

DataType = Extension
print(Model)
IWavelength = [int(w) for w in Wavelength[Model-1]]
print(IWavelength)
#if PressureOption == 3:
#    print(f'PressureFile={BoMPresFile}')
#    FileOK = os.path.isfile(BoMPresFile)
#    if FileOK:
#        open(BoMPresUnit, BoMPresFile, 'r')
#    else:
#        print(f'Error opening Pressure data file {BoMPresFile}')
#        halt()
#
print(f'OzoneFile={OzoneFile}')
FileOK = os.path.isfile(OzoneFile)

FileRoot = '#' + str(Instrument) + str(YearMin % 100).zfill(2) + str(MonthMin).zfill(2) + str(MonthMax).zfill(2)

ParFile = ResultPath + FileRoot, '.lpar'
print(f'par File Name={ParFile}')
#open(ParUnit, ParFile, 'w')
print(ParUnit, f'*** Langley run on {RunDate} at {RunTime}***')
#
FileRoot = '#' + str(Instrument) + str(YearMin % 100).zfill(2) + str(MonthMin).zfill(2) + str(MonthMax).zfill(2)
LtbFile = ResultPath + FileRoot, '.ltb'
print(f'Langley tab File Name={LtbFile}')
CalFile = ResultPath + FileRoot + '.lcl'
print(f'Cal File Name={CalFile}')
#
CutFile = ResultPath + FileRoot + '.cut'
print(f'Cut File Name={CutFile}')

There = os.path.isfile(CutFile)
if There:
    open(CutUnit, CutFile, 'r')
    NumCut, CutList = Read_Cut_Periods(True, CutUnit, NumCut, CutList)
    print(f'Num periods cut={NumCut}')
else:
    print('No valid .cut file found')
    NumCut = 0
NumLangleys = 0   # Count of Langley periods to enable statistcal analysis of entire period
NumWVLangleys = 0 # Count of Water Vapour Langley periods to enable statistcal analysis of entire period
CurrentMonth = 0
for Year in range(YearMin, YearMax+1):
    Year_Code = str(Year)
    MMin = 1
    MMax = 12
    if Year == YearMin:
        MMin = MonthMin
    if Year == YearMax:
        MMax = MonthMax
    for Month in range(MMin, MMax+1):
        Month_Code = str(Month).zfill(2)
        Dmin = 1
        Dmax = (dt.date(Year, Month, 1) - dt.date(Year, Month, 1)).days
        if Month == MonthMin and Year == YearMin:
            Dmin = DayMin
        if Month == MonthMax and Year == YearMax:
            Dmax = DayMax
        for Day in range(Dmin, Dmax+1):
            Day_Code = str(Day).zfill(2)
            print('*************************************')
            print(' Instrument #', Instrument, Year_Code + Month_Code + Day_Code)
            DataPath = DataRoot + DataSuffix + Slash + Config.attrs['Id']+ Slash + '#' + str(Instrument).zfill(2) + Slash + Year_Code + Slash + Month_Code + Slash + Day_Code + Slash
            DateObs = dt.date(Year, Month, Day)
            CalDay = (CalEpoch - DateObs).days
            ObsDate = dt.date(Year, Month, Day)
            if Month != CurrentMonth:
                #Read_BoM_Ozone_2011(Dbug(0), OzoneUnit, ObsDate, Ozone, Eof)
                CurrentMonth = Month
            FileRoot = Year_Code[2:4] + Month_Code + Day_Code
            InFile = DataPath + FileRoot + '.' + Extension
            print('InFile=', InFile)
            #INQUIRE(FILE=InFile, EXIST=FileOK)
            FileOK=True
            if FileOK:
                if PressureOption==1:
                    PresPath = DataRoot + '/PyOut/' + Config.attrs['Id'] + '/'
                    #if Instrument == 2 and Year == 1998 and Month >= 4 and Month <= 7:
                    #    I = INDEX(DataPath, 'tta')
                    #    PresPath = DataPath[1:I-1] + 'ttb/#03' + DataPath[I+7:LEN_TRIM(DataPath)]
                    PresFile = PresPath + FileRoot + '.hpa'
                    #PresFile = PresPath + FileRoot + '.hpa'
                    print('Pressure File Name=', PresFile)
                    #INQUIRE(FILE=PresFile, EXIST=FileOK)
                    if FileOK:
                        #Read_Pressure_File(Dbug(0), PresFile, NumPressurePoints, NumValidPresPoints, PData, DailyMeanPressure)
                        PresColumnNames=['Date','Time','P_UnTCorrect','Tmean','DelP','Pressure']
                        p = pd.read_csv(PresFile, skiprows=9, header=None, delimiter=r'\s+', names=PresColumnNames, index_col=False, parse_dates=[[0,1]]) #depreciated
                        p.rename(columns={'Date_Time':'DateTime'}, inplace=True)
                        print(p)
                        NumPressurePoints = len(p.index)
                        NumValidPresPoints = sum(p.Pressure>0)
                        DailyMeanPressure = sum(p.Pressure[p.Pressure>0])/NumValidPresPoints
                        if NumValidPresPoints <= MinValidPresPoints:
                            print('Insufficient valid pressure data in file')
                    else:
                        NumPressurePoints = 0
                        print('Pressure file not found (.hpa)-has envcal been run?')
                        
                if NumPressurePoints <= 0:
                    print('Assuming Pdefault =', DefaultSurfacePressure)
#                CASE(3):
#                    Read_BoM_Pressure(Dbug(0), BoMPresUnit, ObsDate, PData, NumPressurePoints, NumValidPresPoints, DailyMeanPressure)
            if OzoneOption == 0:
                OzoneColumn = DefaultOzone
#            else:
#                SELECT CASE(OzoneOption)
#                    CASE(1):
#                        if Ozone(Day) > 0:
#                            OzoneColumn = Ozone(Day) / 1000.0
#                        else:
#                            print('Ozone file: Missing day, using monthly mean=', Ozone(0))
#                            OzoneColumn = Ozone(0) / 1000.0
#                    CASE(2):
#                        OzoneColumn = Ozone(0) / 1000.0
            OzoneOD[:NumChannels] = [x*OzoneColumn for x in OzoneCoef[Model-1][:NumChannels]]
            print('Ozone(atm-cm), od(670):', OzoneColumn, OzoneOD[I670])
#            if ApplyClockFix:
#                Get_Time_Correction(Dbug(0), ClockUnit, ObsDate, TimeCorr, Eof)
            BlkFile = DataPath + FileRoot + '.blk'
            print('BlackFile=', BlkFile)
            BlackSun = [0 for _ in BlackSun]
 

            Irec = 0
            Eof = 0
            NumSunTriples = [0, 0, 0]  # Counter for sun triple
            K = 0
            NumSunSingles = [0, 0, 0]  # Counter for sun singles
            L = 0
            print('Processing ' + DataType + ' records...')
            if DataType == 'sun':
                #Read_Triple_Sun_Record(Dbug(1), InUnit, Model, DateTime, Sun3Data, ValidData, Eof)
                Sun3Data = read_triple_sun_record(InFile,Model)
            #elif DataType == 'lsu':
            #    Read_Single_Sun_Record(Dbug(1), InUnit, Model, DateTime, SunData, ValidData, Eof)
            #    SunData = SunData
            #    ValidData = ValidData
            #elif DataType == 'alm':
            #    Read_Almucantar_Record(Dbug(1), InUnit, Model, AlmData, Eof)
            #elif DataType == 'aln':
            #    Read_Almucantar_Split_Record(Dbug(0), InUnit, Model, AlnData, Eof)
            #elif DataType == 'pp1':
            #    Read_PP_Record(Dbug(1), InUnit, Model, PPData, Eof)
            #elif DataType == 'ppp':
            #    Read_PPP_Record(Dbug(1), InUnit, PPPData, Eof)
            #else:
            #    print(' Data type ', DataType, ' not supported yet')
            #while Eof == 0:

            StartLocalDay = dt.datetime(Sun3Data.iloc[0].Date.year, Sun3Data.iloc[0].Date.month, Sun3Data.iloc[0].Date.day, 4, 0, 0) + TimeZone_DT
            EndLocalDay = dt.datetime(Sun3Data.iloc[0].Date.year, Sun3Data.iloc[0].Date.month, Sun3Data.iloc[0].Date.day, 21, 0, 0) + TimeZone_DT

            
            
            SolarZenithDeg = np.full((MaxTriplets, 3), None)
            SolarZenithApp = np.full((MaxTriplets, 3), None)
            SolarAzimuthDeg = np.full((MaxTriplets, 3), None)
            Temperature = np.full(MaxTriplets, None)
            Pressure = np.full(MaxTriplets, None)
            DSunEarth = np.full(MaxTriplets, None)
            AirMassODRayleigh = np.full((MaxTriplets, 3, NumChannels), None)
            AirMassODAerosol = np.full((MaxTriplets, 3, NumChannels), None)
            AirMassODOzone = np.full((MaxTriplets, 3, NumChannels), None)
            AirMass = np.full((MaxTriplets, 3, NumChannels), None)
            Signal = np.full((MaxTriplets, 3, NumChannels), None)
            LnSignal = np.full((MaxTriplets, 3, NumChannels), None)
            TripletMean = np.full((MaxTriplets, NumChannels), None)
            TripletSdev = np.full((MaxTriplets, NumChannels), None)
            TripletCv = np.full((MaxTriplets, NumChannels), None)
            TripletLog = np.full((MaxTriplets, NumChannels), None)
            
            
            for rowindex,row in Sun3Data.iterrows():
                DateTime = dt.datetime.combine(row.Date,row.Time)
                #print(DateTime)
                if (DateTime > StartLocalDay) or (DateTime < EndLocalDay):
            #    while (DateTime < StartLocalDay or DateTime > EndLocalDay) and Eof == 0:
            #        print('Time outside local day, skipping')
            #        print('Startday:', StartLocalDay)
            #        print('Endday  :', EndLocalDay)
            #        print('ObsDT   :', DateTime)
            #        if DataType == 'sun':
            #            Read_Triple_Sun_Record(Dbug(1), InUnit, Model, DateTime, Sun3Data, ValidData, Eof)
            #        elif DataType == 'lsu':
            #            Read_Single_Sun_Record(Dbug(1), InUnit, Model, DateTime, SunData, ValidData, Eof)
            #            SunData = SunData
            #            ValidData = ValidData
            #        elif DataType == 'alm':
            #            Read_Almucantar_Record(Dbug(1), InUnit, Model, AlmData, Eof)
            #        elif DataType == 'aln':
            #            Read_Almucantar_Split_Record(Dbug(0), InUnit, Model, AlnData, Eof)
            #        elif DataType == 'pp1':
            #            Read_PP_Record(Dbug(1), InUnit, Model, PPData, Eof)
            #        elif DataType == 'ppp':
            #            Read_PPP_Record(Dbug(1), InUnit, PPPData, Eof)
            #        else:
            #            print(' Data type ', DataType, ' not supported yet')
            #    if Eof != 0:
            #        break
                    if DataType == 'sun':
                         #YMDHMS[:, 1] = [DateTime.year, DateTime.month, DateTime.day, DateTime.hour, DateTime.minute, DateTime.second, 0]
                #        if ApplyClockFix:
                #            Add_Time(YMDHMS[:, 1], -TimeCorr)
                #            WRITE(FixTimeCheckUnit, 26) DateTime, TimeCorr, YMDHMS[1:6, 1]
                         
                         NumSunTriples[0] = K+1
                         #print(DateTime)
                         DateTimeI=DateTime
                         JulDay, TimeDay = nav.julian(DateTimeI)
                         SunPos = sol.solar_position_almanac(JulDay, TimeDay)
                         #print('JD = {0}, TimeDay = {1}'.format(JulDay,TimeDay))
                         #print('geossun = {}, dsunearth = {}, rA = {}, Dec = {}, EqnTime = {}'.format(SunPos.GeoSun,SunPos.DSunEarth,SunPos.RightAscension,
                         #    SunPos.Declination, SunPos.EqnTime))
                         DSunEarth = SunPos.DSunEarth
            #        if Dbug(1):
            #            print('Year             :', YMDHMS[1, 1])
            #            print('Month            :', YMDHMS[2, 1])
            #            print('Day              :', YMDHMS[3, 1])
            #            print('Hour             :', YMDHMS[4, 1])
            #            print('Minute           :', YMDHMS[5, 1])
            #            print('Second           :', YMDHMS[6, 1])
            #            print('Julian Day       :', JulDay[1])
            #            print('TimeDay          :', TimeDay[K, 1])
            #            print('DSunEarth        :', DSunEarth[K])
            #            HALT()
                         I = 0
                         NewStation = True
                         SolarZenith, SolarAzimuth, NewStation = nav.satellite_position(JulDay, TimeDay, NewStation, StationLatRad, StationLonRad, nav.sun_elements)
                         SolarZenithDeg[K,I] = SolarZenith * DegreesPerRadian
                         SolarZenithApp[K,I] = sol.apparent_zenith(SolarZenithDeg[K,I])
                         SolarAzimuthDeg[K,I] = SolarAzimuth * DegreesPerRadian
                    #if Dbug[1]:
                    #     print('Latitude      :', StationLatDeg)
                    #     print('Longitude     :', StationLonDeg)
                    #     print('Sol Zen (true):', SolarZenithDeg[K,I])
                         #print('Sol Zen (app) :', SolarZenithApp[K,I])
                         #print('Solar Azimuth :', SolarAzimuthDeg[K,I])
                         
            #            HALT()
                    for I in range(1, 3):
                        #YMDHMS[:, I] = YMDHMS[:, I - 1]
                        #Add_Time(YMDHMS[1, I], 30)
                        DateTimeI = DateTimeI + dt.timedelta(seconds=30)
                        JulDay, TimeDay = nav.julian(DateTimeI)
                        SunPos = sol.solar_position_almanac(JulDay, TimeDay)
                                    #            if Dbug(1):
            #                print('Day                 :', DateTime.Day)
            #                print('Month               :', DateTime.Month)
            #                print('Year                :', DateTime.Year)
            #                print('Julian Day number   :', JulDay[I])
            #                print('Fractional Day      :', TimeDay[K,I])
            #                HALT()
                        NewStation = True
                        SolarZenith, SolarAzimuth, NewStation = nav.satellite_position(JulDay, TimeDay, NewStation, StationLatRad, StationLonRad, nav.sun_elements)
                        SolarZenithDeg[K,I] = SolarZenith * DegreesPerRadian
                        SolarZenithApp[K,I] = sol.apparent_zenith(SolarZenithDeg[K,I])
                        SolarAzimuthDeg[K,I] = SolarAzimuth * DegreesPerRadian
                        #Satellite_Position(JulDay[I], TimeDay[K,I], NewStation, StationLatRad, StationLonRad, SunElements, SolarZenith, SolarAzimuth)
                        #SolarZenithDeg[K,I] = SolarZenith * DegreesPerRadian
                        #SolarZenithApp[K,I] = Apparent_Zenith(SolarZenithDeg[K,I])
                        #SolarAzimuthDeg[K,I] = SolarAzimuth * DegreesPerRadian
            #if Dbug(1):
            #    print('Latitude      :', StationLatDeg)
            #    print('Longitude     :', StationLonDeg)
            #    print('Sol Zen (true):', SolarZenithDeg[K,I])
                        #print('Sol Zen (app) :', SolarZenithApp[K,I])
                        #print('Solar Azimuth :', SolarAzimuthDeg[K,I])
            #    HALT()
                    if SolarAzimuthDeg[K][2] < 180.0:
                        NumSunTriples[1] += 1  # Morning
                        NumSunSingles[1] += 3  # Morning
                    else:
                        NumSunTriples[2] += 1  # Afternoon
                        NumSunSingles[2] += 3  # Afternoon
                    
                    Pressure[K] = GetPressure(DateTime, p, PressureOption, NumPressurePoints, NumValidPresPoints, DailyMeanPressure, DefaultSurfacePressure)
                    
                    for N in range(NumChannels):
                        RayleighOD[N] = Rayleigh(Dbug[2], Wavelength[Model-1][N], Pressure[K])
                        #if Dbug[2]:
                        #print('Wavelength:', Wavelength[Model-1][N])
                        #print('Pressure  :', Pressure[K])
                        #print('RayleighOD:', RayleighOD[N])
                    for N in range(NumChannels):
                        for I in range(3):
                            AirMassODRayleigh[K,I,N] = GetAirMass(1, SolarZenithApp[K,I]) * RayleighOD[N]
                            AirMassODAerosol[K,I,N] = GetAirMass(2, SolarZenithApp[K,I]) * 0.03
                            AirMassODOzone[K,I,N] = GetAirMass(3, SolarZenithApp[K,I]) * OzoneOD[N]
                            AirMass[K,I,N] = (AirMassODRayleigh[K,I,N] + AirMassODAerosol[K,I,N] + AirMassODOzone[K,I,N]) / (RayleighOD[N] + 0.03 + OzoneOD[N])
                            if Dbug[2]:
                                #print('Time          :', YMDHMS[1:6, I])
                                print('Sol Zen (true):', SolarZenithDeg[K,I])
                                print('Sol Zen (app) :', SolarZenithApp[K,I])
                                print('Rayleigh am   :', GetAirMass(1, SolarZenithApp[K,I]))
                                print('Aerosol  am   :', GetAirMass(2, SolarZenithApp[K,I]))
                                print('Ozone    am   :', GetAirMass(3, SolarZenithApp[K,I]))
                                print('Weighted am   :', AirMass[K,I,N])
                            
                    for N in range(NumChannels):
                        #print(Sun3Data.iloc[K]['Ch'+str(N+1)][0])
                        Triplet = [(Sun3Data.iloc[K]['Ch'+str(N+1)][0] - BlackSun[N]), (Sun3Data.iloc[K]['Ch'+str(N+1)][1] - BlackSun[N]), (Sun3Data.iloc[K]['Ch'+str(N+1)][2] - BlackSun[N])]
                        TripletMean[K,N], TripletSdev[K,N], tripmin, tripmax = st.stat(Triplet)
                        if TripletSdev[K,N] >= 1000.0:
                            TripletSdev[K,N] = 999.99
                        if TripletMean[K,N] > 0.0:
                            TripletCv[K,N] = 100 * TripletSdev[K,N] / TripletMean[K,N]
                        else:
                            TripletCv[K,N] = 99.99
                        if TripletCv[K,N] >= 99.99:
                            TripletCv[K,N] = 99.99
                        if TripletMean[K,N] > 0.0:
                            TripletLog[K,N] = math.log(TripletMean[K,N])
                        else:
                            TripletLog[K,N] = -9.99
                        #if Dbug[1]:
                        #print('In Channel loop, Cha #=', N)
                        #print('Triplet(1)=', Triplet[1])
                        #print('TripletMean=', TripletMean[K,N])
                        #print('TripletSdev=', TripletSdev[K,N])
                    for I in range(3):
                        L += 1
                        if L > MaxTriplets:
                            print('Array bounds exceeded')
                            print('Increase MaxTriplets!')
                        NumSunSingles[0] = L
                        for N in range(NumChannels):
                            Signal[K,I,N] = Sun3Data.iloc[K]['Ch'+str(N+1)][I] - BlackSun[N]
                            if Signal[K,I,N] > 0.0:
                                LnSignal[K,I,N] = math.log(Signal[K,I,N])
                            else:
                                LnSignal[K,I,N] = -9.99
                    Temperature[K] = Sun3Data.iloc[K].Temperature
                    K += 1

            # End of K loop
            Temperature = Temperature[0:K]
            AirMass = AirMass[0:K,:,:]
            Pressure = Pressure[0:K]

            #if MakeTab:
            #    WRITE(OutUnit, 6) (YMDHMS[J, 2], J=1, 6), Pressure[K], SolarZenithApp[K, 2], AirMass[K, 2, I870], (TripletMean[K, J], TripletSdev[K, J], TripletCv[K, J], TripletLog[K, J], J=1, NumLangleyChannels(Model)), (TripletMean[K, J], TripletSdev[K, J], TripletCv[K, J], TripletLog[K, J], J=IWV, IWV)
            #    if Dbug(1):
            #        WRITE(*          , 6) (YMDHMS[J, 2], J=1, 6), SolarZenithApp[K, 2], AirMass[K, 2, I870], (TripletMean[K, J], TripletSdev[K, J], TripletCv[K, J], TripletLog[K, J], J=1, 1)
            #    for I in range(1, 4):
            #        WRITE(LangleyUnit, 8) (YMDHMS[J, I], J=1, 6), Pressure[K], SolarZenithApp[K,I], AirMass[K, I, I870], (Signal[K, I, J], J=1, NumLangleyChannels(Model)), (Signal[K, I, J], J=IWV, IWV)

            #elif DataType == 'alm':
            #    K += 1  # Index of this almucantar measurement
            #    DateTime = AlmData[4].DateTime  # Time at start of first scan
            #    SkyPath = TRIM(SkyRoot) + Year_Code + Slash + Month_Code + Slash + Day_Code + Slash
            #    OutFile = TRIM(SkyPath) + TRIM(Integer_to_String(DateTime.Day, 2)) + TRIM(Integer_to_String(DateTime.Hour, 2)) + TRIM(Integer_to_String(DateTime.Minute, 2)) + TRIM(Integer_to_String(DateTime.Second, 2)) + '.alm'
            #print('OutFile=', OutFile)
            #HALT()
            #Next_Unit_Number(OutUnit)
            #OPEN(UNIT=OutUnit, FILE=OutFile)
            #YMDHMS[:, 1] = [DateTime.Year, DateTime.Month, DateTime.Day, DateTime.Hour, DateTime.Minute, DateTime.Second, 0]
            #if ApplyClockFix:
            #    Add_Time(YMDHMS[:, 1], -TimeCorr)
            #    WRITE(FixTimeCheckUnit, 26) DateTime, TimeCorr, YMDHMS[1:6, 1]
            #Julian(YMDHMS[1:7, 1], JulDay[1], TimeDay[K, 1])
            #Solar_Position_Almanac(JulDay[1], TimeDay[K, 1], SunPos)
            #DSunEarth[K] = SunPos.DsunEarth
            #if Dbug(1):
            #    print('Year             :', YMDHMS[1, 1])
            #    print('Month            :', YMDHMS[2, 1])
            #    print('Day              :', YMDHMS[3, 1])
            #    print('Hour             :', YMDHMS[4, 1])
            #    print('Minute           :', YMDHMS[5, 1])
            #    print('Second           :', YMDHMS[6, 1])
            #    print('Julian Day       :', JulDay[1])
            #    print('TimeDay          :', TimeDay[K, 1])
            #    print('DSunEarth        :', DSunEarth[K])
            #    HALT()
            #I = 1
            #NewStation = True
            #Satellite_Position(JulDay[I], TimeDay[K,I], NewStation, StationLatRad, StationLonRad, SunElements, SolarZenith, SolarAzimuth)
            #SolarZenithDeg[K,I] = SolarZenith * DegreesPerRadian
            #SolarZenithApp[K,I] = Apparent_Zenith(SolarZenithDeg[K,I])
            #SolarAzimuthDeg[K,I] = SolarAzimuth * DegreesPerRadian
            #if Dbug(1):
            #    print('Latitude      :', StationLatDeg)
            #    print('Longitude     :', StationLonDeg)
            #    print('Sol Zen (true):', SolarZenithDeg[K,I])
            #    print('Sol Zen (app) :', SolarZenithApp[K,I])
            #    print('Solar Azimuth :', SolarAzimuthDeg[K,I])
            #    HALT()
            #WRITE(OutUnit, 9) InFile, DateTime.Year, DateTime.Month, DateTime.Day, DateTime.Hour, DateTime.Minute, DateTime.Second, SolarZenithApp[K, 1], (Wavelength_Aerosol[I], I=1, NumAerosolChannels)
            #for N in range(1, NumAzimuths+1):
            #    WRITE(OutUnit, 10) N, '', AlmAzimuth[N], CollGainAlm[N], (AlmData[J].Signal[N], J=1, NumAerosolChannels)
            #CLOSE(UNIT=OutUnit)
            #elif DataType == 'aln':
            #    K += 1  # Index of this almucantar measurement
            #    DateTime = AlnData[4].DateTime  # Time at start of first scan
            #    SkyPath = TRIM(SkyRoot) + Year_Code + Slash + Month_Code + Slash + Day_Code + Slash
            #    OutFile = TRIM(SkyPath) + TRIM(Integer_to_String(DateTime.Day, 2)) + TRIM(Integer_to_String(DateTime.Hour, 2)) + TRIM(Integer_to_String(DateTime.Minute, 2)) + TRIM(Integer_to_String(DateTime.Second, 2)) + '.alm'
            #    print('OutFile=', OutFile)
            #    HALT()
            #    Next_Unit_Number(OutUnit)
            #    OPEN(UNIT=OutUnit, FILE=OutFile)
            #    YMDHMS[:, 1] = [DateTime.Year, DateTime.Month, DateTime.Day, DateTime.Hour, DateTime.Minute, DateTime.Second, 0]
            #    if ApplyClockFix:
            #        Add_Time(YMDHMS[:, 1], -TimeCorr)
            #        WRITE(FixTimeCheckUnit, 26) DateTime, TimeCorr, YMDHMS[1:6, 1]
            #    Julian(YMDHMS[1:7, 1], JulDay[1], TimeDay[K, 1])
            #    Solar_Position_Almanac(JulDay[1], TimeDay[K, 1], SunPos)
            #    DSunEarth[K] = SunPos.DsunEarth
            #    if Dbug(1):
            #        print('Year             :', YMDHMS[1, 1])
            #        print('Month            :', YMDHMS[2, 1])
            #        print('Day              :', YMDHMS[3, 1])
            #        print('Hour             :', YMDHMS[4, 1])
            #        print('Minute           :', YMDHMS[5, 1])
            #        print('Second           :', YMDHMS[6, 1])
            #        print('Julian Day       :', JulDay[1])
            #        print('TimeDay          :', TimeDay[K, 1])
            #        print('DSunEarth        :', DSunEarth[K])
            #        HALT()
            #    I = 1
            #    NewStation = True
            #    Satellite_Position(JulDay[I], TimeDay[K,I], NewStation, StationLatRad, StationLonRad, SunElements, SolarZenith, SolarAzimuth)
            #    SolarZenithDeg[K,I] = SolarZenith * DegreesPerRadian
            #    SolarZenithApp[K,I] = Apparent_Zenith(SolarZenithDeg[K,I])
            #    SolarAzimuthDeg[K,I] = SolarAzimuth * DegreesPerRadian
            #    if Dbug(1):
            #        print('Latitude      :', StationLatDeg)
            #        print('Longitude     :', StationLonDeg)
            #        print('Sol Zen (true):', SolarZenithDeg[K,I])
            #        print('Sol Zen (app) :', SolarZenithApp[K,I])
            #        print('Solar Azimuth :', SolarAzimuthDeg[K,I])
            #        HALT()
            #    WRITE(OutUnit, 9) InFile, DateTime.Year, DateTime.Month, DateTime.Day, DateTime.Hour, DateTime.Minute, DateTime.Second, SolarZenithApp[K, 1], (Wavelength_Aerosol[I], I=1, NumAerosolChannels)
            #    # Right hemisphere, 0-180 azimuth
            #    for N in range(1, NumNAzimuths+1):
            #        WRITE(OutUnit, 10) N, 'R', AlmRightAzimuth[N], CollGainAlmN[N], (AlnData[J].SignalRight[N], J=1, NumAerosolChannels)
            #    # Left hemisphere, 180-360 azimuth
            #    for N in range(1, NumNAzimuths+1):
            #        WRITE(OutUnit, 10) N, 'L', AlmLeftAzimuth[N], CollGainAlmN[N], (AlnData[J].SignalLeft[N], J=1, NumAerosolChannels)
            #    CLOSE(UNIT=OutUnit)
            #elif DataType == 'pp1':
            #    # 27/4/2006       Read all 4 wavelengths
            #    K += 1  # Index of this measurement
            #    DateTime = PPData[4].DateTime  # Time at start of this scan
            #    SkyPath = TRIM(SkyRoot) + Year_Code + Slash + Month_Code + Slash + Day_Code + Slash
            #    SkyFile = TRIM(SkyPath) + TRIM(Integer_to_String(DateTime.Day, 2)) + TRIM(Integer_to_String(DateTime.Hour, 2)) + TRIM(Integer_to_String(DateTime.Minute, 2)) + TRIM(Integer_to_String(DateTime.Second, 2)) + '.pp'
            #print('PP SkyFile=', SkyFile)
            #HALT()
            #Next_Unit_Number(OutUnit)
            #OPEN(UNIT=OutUnit, FILE=SkyFile)
            #YMDHMS[:, 1] = [DateTime.Year, DateTime.Month, DateTime.Day, DateTime.Hour, DateTime.Minute, DateTime.Second, 0]
            #if ApplyClockFix:
            #    Add_Time(YMDHMS[:, 1], -TimeCorr)
            #    WRITE(FixTimeCheckUnit, 26) DateTime, TimeCorr, YMDHMS[1:6, 1]
            #Julian(YMDHMS[1:7, 1], JulDay[1], TimeDay[K, 1])
            #Solar_Position_Almanac(JulDay[1], TimeDay[K, 1], SunPos)
            #DSunEarth[K] = SunPos.DsunEarth
            #if Dbug(1):
            #    print('Year             :', YMDHMS[1, 1])
            #    print('Month            :', YMDHMS[2, 1])
            #    print('Day              :', YMDHMS[3, 1])
            #    print('Hour             :', YMDHMS[4, 1])
            #    print('Minute           :', YMDHMS[5, 1])
            #    print('Second           :', YMDHMS[6, 1])
            #    print('Julian Day       :', JulDay[1])
            #    print('TimeDay          :', TimeDay[K, 1])
            #    print('DSunEarth        :', DSunEarth[K])
            #    HALT()
            #I = 1
            #NewStation = True
            #Satellite_Position(JulDay[I], TimeDay[K,I], NewStation, StationLatRad, StationLonRad, SunElements, SolarZenith, SolarAzimuth)
            #SolarZenithDeg[K,I] = SolarZenith * DegreesPerRadian
            #SolarZenithApp[K,I] = Apparent_Zenith(SolarZenithDeg[K,I])
            #SolarAzimuthDeg[K,I] = SolarAzimuth * DegreesPerRadian
            #if Dbug(1):
            #    print('Latitude      :', StationLatDeg)
            #    print('Longitude     :', StationLonDeg)
            #    print('Sol Zen (true):', SolarZenithDeg[K,I])
            #    print('Sol Zen (app) :', SolarZenithApp[K,I])
            #    print('Solar Azimuth :', SolarAzimuthDeg[K,I])
            #    HALT()
            #WRITE(OutUnit, 91) InFile, DateTime.Year, DateTime.Month, DateTime.Day, DateTime.Hour, DateTime.Minute, DateTime.Second, SolarZenithApp[K, 1], (Wavelength_Aerosol[I], I=1, NumAerosolChannels)
            #for N in range(1, NumPPZeniths+1):
            #    WRITE(OutUnit, 10) N, '', PPZenith[N], CollGainPP[N], (PPData[J].Signal[N], J=1, NumAerosolChannels)
            #CLOSE(UNIT=OutUnit)
            #elif DataType == 'ppp':
            #    K += 1  # Index of this measurement
            #    DateTime = PPPData.DateTime  # Time at start of first scan
            #    OutFile = TRIM(ResultPath) + TRIM(Integer_to_String(DateTime.Day, 2)) + TRIM(Integer_to_String(DateTime.Hour, 2)) + TRIM(Integer_to_String(DateTime.Minute, 2)) + TRIM(Integer_to_String(DateTime.Second, 2)) + '.ppp'
            #    YMDHMS[:, 1] = [DateTime.Year, DateTime.Month, DateTime.Day, DateTime.Hour, DateTime.Minute, DateTime.Second, 0]
            #    if ApplyClockFix:
            #        Add_Time(YMDHMS[:, 1], -TimeCorr)
            #        WRITE(FixTimeCheckUnit, 26) DateTime, TimeCorr, YMDHMS[1:6, 1]
            #    Julian(YMDHMS[1:7, 1], JulDay[1], TimeDay[K, 1])
            #    Solar_Position_Almanac(JulDay[1], TimeDay[K, 1], SunPos)
            #    DSunEarth[K] = SunPos.DsunEarth
            #    if Dbug(1):
            #        print('Year             :', YMDHMS[1, 1])
            #        print('Month            :', YMDHMS[2, 1])
            #        print('Day              :', YMDHMS[3, 1])
            #        print('Hour             :', YMDHMS[4, 1])
            #        print('Minute           :', YMDHMS[5, 1])
            #        print('Second           :', YMDHMS[6, 1])
            #        print('Julian Day       :', JulDay[1])
            #        print('TimeDay          :', TimeDay[K, 1])
            #        print('DSunEarth        :', DSunEarth[K])
            #        HALT()
            #    I = 1
            #    NewStation = True
            #    Satellite_Position(JulDay[I], TimeDay[K,I], NewStation, StationLatRad, StationLonRad, SunElements, SolarZenith, SolarAzimuth)
            #    SolarZenithDeg[K,I] = SolarZenith * DegreesPerRadian
            #    SolarZenithApp[K,I] = Apparent_Zenith(SolarZenithDeg[K,I])
            #    SolarAzimuthDeg[K,I] = SolarAzimuth * DegreesPerRadian
            #    if Dbug(1):
            #        print('Latitude      :', StationLatDeg)
            #        print('Longitude     :', StationLonDeg)
            #        print('Sol Zen (true):', SolarZenithDeg[K,I])
            #        print('Sol Zen (app) :', SolarZenithApp[K,I])
            #        print('Solar Azimuth :', SolarAzimuthDeg[K,I])
            #        HALT()
            #    WRITE(OutUnit, 92) InFile, Instrument, DateTime.Year, DateTime.Month, DateTime.Day, DateTime.Hour, DateTime.Minute, DateTime.Second, SolarZenithApp[K, 1]
            #    for N in range(1, NumPPPZeniths+1):
            #        WRITE(OutUnit, 10) N, PPPZenith[N], CollGainPPP[N], (PPPData.Signal[J, N], J=1, 3)
            #    CLOSE(UNIT=OutUnit)
            #else:
            #    print(' Data type ', DataType, ' not supported yet')
            #    HALT()
            #if DataType == 'sun':
            #    Read_Triple_Sun_Record(Dbug(1), InUnit, Model, DateTime, Sun3Data, ValidData, Eof)
            #elif DataType == 'lsu':
            #    Read_Single_Sun_Record(Dbug(1), InUnit, Model, DateTime, SunData, ValidData, Eof)
            #    SunData = SunData
            #    ValidData = ValidData
            #elif DataType == 'alm':
            #    Read_Almucantar_Record(Dbug(1), InUnit, Model, AlmData, Eof)
            #elif DataType == 'aln':
            #    Read_Almucantar_Split_Record(Dbug(0), InUnit, Model, AlnData, Eof)
            #elif DataType == 'pp1':
            #    Read_PP_Record(Dbug(1), InUnit, Model, PPData, Eof)
            #elif DataType == 'ppp':
            #    Read_PPP_Record(Dbug(1), InUnit, PPPData, Eof)
            #else:
            #    print(' Data type ', DataType, ' not supported yet')
            #
            #
            #    if DataType == 'sun':
            #        InUnit.close()
            #        if MakeTab:
            #            OutUnit.close()
            #            LangleyUnit.close()
            #    else:
            #        InUnit.close()
            #
            LangFlag[0:2] = [False, False]
            if DataType == 'sun' or DataType == 'lsu':
                print('Langley processing...')
                for I in [1,2]:
                    Nstart[I-1] = NumSunTriples[1] * (I - 1)
                    Nend[I-1] = Nstart[I-1] + NumSunTriples[I] 
                    #print(NumSunTriples,Nstart,Nend)
            #        if Dbug[2]:
            #            print('Number of sun triples(am)  :', NumSunTriples[1])
            #            print('Number of sun triples(pm)  :', NumSunTriples[2])
            #            print('Number of sun triples(tot) :', NumSunTriples[0])
            #            print('Number of sun singles(am)  :', NumSunSingles[1])
            #            print('Number of sun singles(pm)  :', NumSunSingles[2])
            #            print('Number of sun singles(tot) :', NumSunSingles[0])
            #            print('Nstart (am)                :', Nstart[1])
            #            print('Nend   (am)                :', Nend[1])
            #            print('Nstart (pm)                :', Nstart[2])
            #            print('Nend   (pm)                :', Nend[2])
            #            HALT()
            #
                FitFile = ResultPath + 'Langley.dat'
        #        open(FitFile, 'w')
        #
                for Iap in [0,1]:  # Morning/Afternoon
                    ObsDate = dt.date(Year, Month, Day)
                    Include = True
                    for Ic in range(NumCut):
                        if ObsDate == CutList[Ic].Date and AmPm[Iap] == CutList[Ic].AP:
                            Include = False
        
                    #if Include:
                    #    print(f"#{Year:5d}{Month:5d}{Day:5d} {AmPm[Iap]}") # write to fitfile
        
                    PresMean[Iap] = 0
                    N = 0
                    for K in range(Nstart[Iap], Nend[Iap]):
                        if Pressure[K] > 0:
                            N += 1
                            PresMean[Iap] += Pressure[K]
                    PresMean[Iap] /= N
            
                    for N in range(NumChannels):
                        RayleighOD[N] = Rayleigh(Dbug[2], Wavelength[Model-1][N], PresMean[Iap])
                        if Dbug[2]:
                            print('Wavelength   :', Wavelength[N, Model])
                            print('Mean Pressure:', PresMean[Iap])
                            print('RayleighOD   :', RayleighOD[N])
        
                    if NumSunTriples[Iap] >= MinSunTriples:
                        SpreadFlag, NumOK[Iap], IndOK = CheckTripletCv(Dbug[2], NumSunTriples[Iap+1], Nstart[Iap],AirMass[:,1,I870],TripletCv[:,I870])
                                  
                        if NumOK[Iap] >= MinSunTriples and SpreadFlag:
                            print(' Passed triplet cv langley filter')
                            N = NumChannels
                            I = 0
                            for K in range(NumOK[Iap]):
                                for J in range(3):
                                    I += 1
                                    X[I, 0:N] = Airmass[IndOK[K, Iap], J, 0:N]
                                    Y[I, 0:N] = LnSignal[IndOK[K, Iap], J, 0:N]
            
                            NumPoints[Iap] = I
            
                            CheckFitQuality(Dbug[2],
                                        NumPoints[Iap],
                                        X, Y,
                                        I870,
                                        MaxSdevFit,
                                        SpreadFlag,
                                        FitFlag)
            
                            if SpreadFlag and FitFlag:
                                LangFlag[Iap] = True
                                print(' Passed fit quality langley filter')
                                NumLangleys += 1
        
                                for N in range(NumLangleyChannels[Model-1]):
                                    print(' Wavelength #', N)
                                    BoxFit(NumPoints[Iap],
                                            X[1, N], Y[1, N],
                                            Intercept[Iap, N],
                                            Slope[Iap, N],
                                            Residual[N],
                                            Erms[N],
                                            DelIntercept[Iap, N],
                                            DelSlope[Iap, N])
        
                                    Intercept[Iap, N] += 2.0 * log(DSunEarth[NStart[0]])
        
                                    AerosolOD[Iap, N] = -Slope[Iap, N] - RayleighOD[N] - OzoneOD[N]
                                    LnV0Glob[NumLangleys, N] = Intercept[Iap, N]
                                    AodGlob[NumLangleys, N] = AerosolOD[Iap, N]
            
                                for N in range(NumLangleyChannels[Model-1]):
                                    for M in range(N, NumLangleyChannels[Model-1]):
                                        if AerosolOD[Iap, N] >= 0 and AerosolOD[Iap, M] >= 0:
                                            Angstrom[Iap, N, M] = -log(AerosolOD[Iap, N] / AerosolOD[Iap, M]) / \
                                                    log(Wavelength[N, Model] / Wavelength[M, Model])
                                        else:
                                            Angstrom[Iap, N, M] = -9.999
        
                                AerosolOD[Iap, IWV] = AerosolOD[Iap, I870] * \
                                        (Wavelength[IWV][Model-1] / Wavelength[I870][Model-1]) ** \
                                        (-Angstrom[Iap, I670, I870])
            
                                for I in range(NumPoints[Iap]):
                                    U[I] = X[I, IWV] ** Bcoef
                                    V[I] = Y[I, IWV] + X[I, IWV] * (RayleighOD[IWV] + AerosolOD[Iap, IWV])
            
            #                        if Dbug[3]:
            #                            print('Signal    :', Y[I, IWV])
            #                            print('Airmass   :', X[I, IWV])
            #                            print('RayleighOD:', RayleighOD[IWV])
            #                            print('AOD 670   :', AerosolOD[Iap, I670])
            #                            print('AOD 870   :', AerosolOD[Iap, I870])
            #                            print('AOD 936   :', AerosolOD[Iap, IWV])
            #                            print('U         :', U[I])
            #                            print('V         :', V[I])
            #                            HALT()
            #
                                Weight[0:NumPoints[Iap]] = 1.0
        
                                ELFIT(NumPoints[Iap],
                                        Weight,
                                        U, V,
                                        Intercept[Iap, IWV],
                                        Slope[Iap, IWV],
                                        Residual[N],
                                        Erms[N],
                                        DelIntercept[Iap, IWV],
                                        DelSlope[Iap, IWV])
        
                                Intercept[Iap, IWV] += 2.0 * log(DSunEarth[NStart[0]])
        
                                if Slope[Iap, IWV] < 0:
                                    NumWVLangleys += 1
                                    WvapMean[Iap] = (-Slope[Iap, IWV] / Acoef) ** (1.0 / Bcoef)
                                    LnV0Glob[NumWVLangleys, IWV] = Intercept[Iap, IWV]
                                    AodGlob[NumWVLangleys, IWV] = WvapMean[Iap]
                                    OD_WvapMean[Iap, I1020] = H2O_k * WvapMean[Iap] ** H2O_e
                                else:
                                    print('Water vapour below detection threshold')
                                    WvapMean[Iap] = 0
                                    OD_WvapMean[Iap, I1020] = 0
        
                                OD_WvapMean[Iap, I1020] = OD_WvapMean[Iap, I1020]
        
                                print('writing record to', FitFile)
                                FitUnit.write(f'# {Wavelength[N, Model]:12.1f} ' +
                                        ' '.join(f'{Wavelength[N, Model]:18.1f}' for N in range(NumLangleyChannels[Model-1])) +
                                        f' {Wavelength[IWV, Model]:18.1f}')
                                FitUnit.write('# ' + '  Airmass   lnV0  ' * 10)
                                for I in range(1, NumPoints[Iap] + 1):
                                    FitUnit.write(' '.join(f'{X[I, N]:9.5f} {Y[I, N]:9.5f}' for N in range(NumLangleyChannels[Model-1])) +
                                            f' {U[I]:9.5f} {V[I]:9.5f}')
        
                                    Langley_Mean_Temperature(Dbug[2],
                                            NumOK[Iap], IndOK[1, Iap],
                                            Temperature,
                                            MeanDetectorTemp[Iap])
                                    DaysSinceEpoch[NumLangleys] = CalDay
                                print(f'{AmPm[Iap]} {PresMean[Iap]:13.1f} {NumOK[Iap]:13d} ' +
                                        ' '.join(f'{Wavelength[N, Model]:13.0f}' for N in range(NumLangleyChannels[Model-1])) +
                                        f' {Wavelength[IWV, Model]:13.0f}')
                                print(' '.join(f'{AerosolOD[Iap, N]:13.4f}' for N in range(NumLangleyChannels[Model-1])) +
                                        f' {WvapMean[Iap]:13.4f}')
                            else:
                                print(' Passed fit quality langley filter')
                        else:
                            print(' Failed triplet cv langley filter')
                    else:
                        print(' Insufficient triplets')
            
            #
            #    if LangFlag[1] or LangFlag[2]:
            #        DayYear = Day_of_Year(Day, Month, Year)
            #        for Iap in range(1, 3):
            #            if LangFlag[Iap]:
            #                LtbUnit.write(FormatString[2].format(
            #                    Year, Month, Day, CalDay, AmPm[Iap],
            #                    PresMean[Iap], NumPoints[Iap],
            #                    *(Intercept[Iap, N], DelIntercept[Iap, N],
            #                        AerosolOD[Iap, N],
            #                        DelSlope[Iap, N],
            #                        for N in range(1, NumLangleyChannels[Model] + 1)),
            #                    Intercept[Iap, IWV], DelIntercept[Iap, IWV],
            #                    WvapMean[Iap],
            #                    DelSlope[Iap, IWV],
            #                    Angstrom[Iap, I440, I870],
            #                    MeanDetectorTemp[Iap]
            #                    ))
            #
            #                AlfUnit.write(f'{Year:4d}{Month:3d}{Day:3d}{DayYear:4d} {AmPm[Iap]} '
            #                        f'{PresMean[Iap]:7.1f} {NumPoints[Iap]:4d} '
            #                        f'{Angstrom[Iap, I440, I670]:9.3f} '
            #                        f'{Angstrom[Iap, I440, I870]:9.3f} '
            #                        f'{Angstrom[Iap, I440, I1020]:9.3f} '
            #                        f'{Angstrom[Iap, I670, I870]:9.3f} '
            #                        f'{Angstrom[Iap, I670, I1020]:9.3f} '
            #                        f'{Angstrom[Iap, I870, I1020]:9.3f}')
            #
            #                LtbUnit.close()
            #    AlfUnit.close()
            #
            #    if NumLangleys >= 2:
            #        for N in range(1, NumLangleyChannels[Model] + 1):
            #            LnV0Mean[N], LnV0Sdev[N] = Stat(LnV0Glob[1, N], NumLangleys)
            #            AodMean[N], AodSdev[N] = Stat(AodGlob[1, N], NumLangleys)
            #    elif NumLangleys == 1:
            #        for N in range(1, NumLangleyChannels[Model] + 1):
            #            LnV0Mean[N] = LnV0Glob[1, N]
            #            LnV0Sdev[N] = 0
            #            AodMean[N] = AodGlob[1, N]
            #            AodSdev[N] = 0
            #
            #    if NumWVLangleys >= 2:
            #        LnV0Mean[IWV], LnV0Sdev[IWV] = Stat(LnV0Glob[1, IWV], NumWVLangleys)
            #        AodMean[IWV], AodSdev[IWV] = Stat(AodGlob[1, IWV], NumWVLangleys)
            #    elif NumWVLangleys == 1:
            #        LnV0Mean[IWV] = LnV0Glob[1, IWV]
            #        LnV0Sdev[IWV] = 0
            #        AodMean[IWV] = AodGlob[1, IWV]
            #        AodSdev[IWV] = 0
            #
            #    if NumLangleys > 0:
            #        LtbUnit.write(FormatString[3].format(
            #            '#Summary', NumLangleys,
            #            *(LnV0Mean[N], LnV0Sdev[N],
            #                AodMean[N], AodSdev[N],
            #                for N in range(1, NumLangleyChannels[Model] + 1)),
            #            LnV0Mean[IWV], LnV, LnV0Sdev[IWV], AodMean[IWV], AodSdev[IWV]))
            #
            #        if FitLineV0 and NumLangleys >= 2:
            #            IOrder = 1
            #        Weight[0:NumLangleys] = 1.0
            #        for N in range(1, NumLangleyChannels[Model] + 1):
            #            Elfit(NumLangleys, Weight,
            #                    DaysSinceEpoch, LnV0Glob[1, N],
            #                    LnV0Coef[0, N], LnV0Coef[1, N],
            #                    Residual[N], Erms[N],
            #                    ELnV0Coef[0, N], ElnV0Coef[1, N])
            #
            #            if Dbug[3]:
            #                print('Channel  number ', N)
            #                print('NumLangleys     ', NumLangleys)
            #                print('Days since epoch', *(DaysSinceEpoch[I] for I in range(1, NumLangleys + 1)))
            #                print('LnV0            ', *(LnV0Glob[I, N] for I in range(1, NumLangleys + 1)))
            #                print('Intercept      ', LnV0Coef[0, N])
            #                print('DelIntercept   ', ELnV0Coef[0, N])
            #                print('Slope          ', LnV0Coef[1, N])
            #                print('DelSlope       ', ELnV0Coef[1, N])
            #                print('Rms error      ', Erms[N])
            #                HALT()
            #
            #        if Dbug[3]:
            #            print('Channel  number ', N)
            #            print('NumLangleys     ', NumLangleys)
            #            print('Days since epoch', *(DaysSinceEpoch[I] for I in range(1, NumLangleys + 1)))
            #            print('LnV0            ', *(LnV0Glob[I, IWV] for I in range(1, NumLangleys + 1)))
            #            print('Intercept      ', LnV0Coef[0, IWV])
            #            print('DelIntercept   ', ELnV0Coef[0, IWV])
            #            print('Slope          ', LnV0Coef[1, IWV])
            #            print('DelSlope       ', ELnV0Coef[1, IWV])
            #            print('Rms error      ', Erms[IWV])
            #            HALT()
            #    elif NumLangleys >= 1:
            #        IOrder = 0
            #        for N in range(1, NumLangleyChannels[Model] + 1):
            #            LnV0Coef[0, N] = LnV0Mean[N]
            #            Erms[N] = LnV0Sdev[N]
            #
            #    if FitLineV0 and NumWVLangleys >= 2:
            #        Elfit(NumWVLangleys, Weight,
            #                DaysSinceEpoch, LnV0Glob[1, IWV],
            #                LnV0Coef[0, IWV], LnV0Coef[1, IWV],
            #                Residual[IWV], Erms[IWV],
            #                ELnV0Coef[0, IWV], ElnV0Coef[1, IWV])
            #    elif NumLangleys >= 1:
            #        LnV0Coef[0, IWV] = LnV0Mean[IWV]
            #        Erms[IWV] = LnV0Sdev[IWV]
            #
            #    if MakeCal and NumLangleys > 0:
            #        if CalFile.exists():
            #            print(f'{CalFile} exists, will be overwritten.')
            #            HALT()
            #
            #        CalUnit = Next_Unit_Number()
            #        open(CalFile, 'w')
            #
            #        NumGeneralCycles = 0
            #        CalUnit.write(f'# Calibration file generated from Langleys between:\n'
            #                f'#{YearMin:5d}{MonthMin:2d}{DayMin:2d} and\n'
            #                f'#{YearMax:5d}{MonthMax:2d}{DayMax:2d}.\n'
            #                f'#{Instrument:5d}        -- Instrument number\n'
            #                f'#{Model:5d}        -- Model number     \n'
            #                f'#{NumLangleyChannels[Model]:5d}        -- Number of Langley wavelengths\n'
            #                f'#{NumLangleys:5d}        -- Number of Langley intervals in period\n'
            #                f'#{NumGeneralCycles:5d}        -- Number of General Method cycles applied\n'
            #                f'#{CalEpoch[0]:5d}{CalEpoch[1]:3d}{CalEpoch[2]:3d}  -- Calibration epoch')
            #
            #        FormatString[5] = f'# Wavel(nm) Order {" ".join(f"LnV0({I})" for I in range(IOrder + 1))} Erms'
            #        CalUnit.write(FormatString[5])
            #
            #        for N in range(1, NumLangleyChannels[Model] + 1):
            #            CalUnit.write(f'{Wavelength[N, Model]:10.1f} {IOrder:6d} ' +
            #                    ' '.join(f'{LnV0Coef[I, N]:13.5e}' for I in range(IOrder + 1)) +
            #                    f' {Erms[N]:13.5e}')
            #
            #            CalUnit.write(f'{Wavelength[IWV, Model]:10.1f} {IOrder:6d} ' +
            #                    ' '.join(f'{LnV0Coef[I, IWV]:13.5e}' for I in range(IOrder + 1)) +
            #                    f' {Erms[IWV]:13.5e}')
            #
            #            CalUnit.close()
            #
            #print('Done.')
