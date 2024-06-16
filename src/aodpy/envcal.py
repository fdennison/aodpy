from cimel_module import *
#from constants_module import *
#from navigation_module import *
from utils_module import *
import time_module as tm
import os as os
import pandas as pd
import math
import datetime as dt

Dbug = [False, False, False]
FileOK = False
Instrument = None
Card = None
Barometer = None
BarometerToUse = None
CardForBarometer = None
Anemometer = None
Year = None
Month = None
Day = None
YearMin = None
YearMax = None
MonthMin = None
MonthMax = None
MMin = None
MMax = None
DayMin = None
DayMax = None
Dmin = None
Dmax = None
InUnit = None
BaroUnit = None
SumUnit = None
HpaUnit = None
EnvUnit = None
ConfigUnit = None
NumDays = None

LoggerVersion = ''
Extension = ''
Station = ''
Year_Code = ''
Month_Code = ''
Day_Code = ''

RootPath = ''
DataPath = ''
AnalogCalPath = ''
ConfigPath = ''
FilePath = ''
SumPath = ''
FileRoot = ''
BaroFile = ''
InFile = ''
SumFile = ''
HpaFile = ''
EnvFile = ''
ConfigFile = ''

Status = None  # Type (Status_)
Env = None  # Type (Raw_Environment_Record)
Pressure = None  # Type (Real_Stat)
Temperature = None  # Type (Real_Stat)
WindV = None  # Type (Real_Stat)
WindD = None  # Type (Real_Stat)

ObsDate = None  # Type (Date)
StartDate = None  # Type (Date)
EndDate = None  # Type (Date)

DateTimeUTC = None  # Type (Date_Time)

I = None
J = None
K = None
N = None
Eof = None
HourIndex = None
NumBaroCals = None
NumWindCals = None
NumTempCals = None
NumVoltCals = None
DTimeSec = None

YMDHMS = [None] * 6
TimeZoneOffset = [None] * 3

MaxCards = 32
MaxBaro = 32
PCalDate = [None] * (MaxCards)
PCalDate0 = None  # Type (Date)
WCalDate = [None] * (MaxCards)  # Type (Date)
TCalDate = [None] * (MaxCards)  # Type (Date)
VCalDate = [None] * (MaxCards)  # Type (Date)

PCoef = [None] * (MaxCards)
PTCoef = [None] * (MaxCards) 
WCoef = [None] * (MaxCards)
VCoef = [None] * (MaxCards) 
TCoef = [None] * (MaxCards) 

Vbatt = None  # Type (LongReal)
DelPressure = None  # Type (LongReal)
Pressure_UnTComp = None  # Type (LongReal)
TZ = None  # Type (Time)
CE = None  # Type (Time)

print("Envcal")
# envcal.men
Dbug[1]=False
Extension='rnv'
RootPath='/home/599/fd0474/AODcode/SampleData'
DataPath='/agsdat'
ConfigPath='/config'
Station='jb1'
AnalogCalPath='/analog.cal'
SumPath='/analysis/tta/sep02'
YearMin=2016
MonthMin=1
DayMin=1
YearMax=2016
MonthMax=1
DayMax=31
CE=dt.timedelta(hours=0,minutes=0,seconds=0)

ConfigPath = RootPath+ConfigPath
print('config path=', ConfigPath)
ConfigFile = ConfigPath+'/'+Station+'.cfn'
print(' Config File=', ConfigFile)

# Get configuration on start date
ObsDate = dt.date(YearMin, MonthMin, DayMin)
Config = Read_Site_Configuration(Dbug[0], ConfigFile, ObsDate)
print(Config)
StartDate = Config.StartDate
EndDate = Config.EndDate
Instrument = Config.CimelNumber
LoggerVersion = 'v' + str(Config.Version).zfill(1)
Card = Config.Card
Barometer = Config.Baro
Anemometer = Config.Wind
TZ = Config.TOffset

    # If there is more than one entry for a given instrument,
    # the last entry is the one used. Ideally there should
    # be discrimination using date of calibration.

if Card >= 0:
    #    Menu.FileName = os.path.join(RootPath, AnalogCalPath, 'tcal.men')
#    # call Menu_Edit(Menu)
#    Menu.Value = read_menu_value()
#    NumTempCals, TempCalData = parse_tcal_menu_value(Menu.Value)
#    for j, tcal_date, tcoef in TempCalData:
#        TCalDate[j] = tcal_date
#      TCoef[0][j], TCoef[1][j] = tcoef
    TCalDate[1]=tm.Date(1998,1,1)    
    TCoef[1]=[-12.9,0.089]
    TCalDate[2]=tm.Date(1998,1,1)    
    TCoef[2]=[-12.9,0.089]
    TCalDate[3]=tm.Date(1998,1,1)    
    TCoef[3]=[-12.9,0.089]
    TCalDate[4]=tm.Date(1998,1,1)    
    TCoef[4]=[-12.9,0.089]
    TCalDate[5]=tm.Date(2009,9,16)    
    TCoef[5]=[-31.35,0.1255]
    TCalDate[6]=tm.Date(1998,1,1)   
    TCoef[6]=[-12.9,0.089]
    TCalDate[7]=tm.Date(1998,1,1)    
    TCoef[7]=[-12.9,0.089]
    TCalDate[8]=tm.Date(1998,1,1)    
    TCoef[8]=[-12.9,0.089]
    TCalDate[9]=tm.Date(1998,1,1)    
    TCoef[9]=[-12.9,0.089]
    TCalDate[10]=tm.Date(1998,1,1)    
    TCoef[10]=[-12.9,0.089]
    TCalDate[11]=tm.Date(1998,1,1)    
    TCoef[11]=[-12.9,0.089]
    TCalDate[12]=tm.Date(1998,1,1)    
    TCoef[12]=[-12.9,0.089]
    TCalDate[21]=tm.Date(1998,1,1)    
    TCoef[21]=[-12.9,0.089]
    TCalDate[22]=tm.Date(1998,1,1)    
    TCoef[22]=[-12.9,0.089]
    TCalDate[23]=tm.Date(1998,1,1)    
    TCoef[23]=[-12.9,0.089]
    TCalDate[24]=tm.Date(1998,1,1)    
    TCoef[24]=[-12.9,0.089]

if Barometer > 0:
    print(' setting baro File')
    BaroFile =RootPath+AnalogCalPath+'/pcal/'+LoggerVersion+'/pcal.dat'
    print(' Baro File=', BaroFile)
    barocolnames=['Year', 'Month', 'Day', 'Baro', 'Card', 'PCoef0', 'PCoef1', 'PTCoef0', 'PTCoef1']
    BaroDF=pd.read_csv(BaroFile, skiprows=4, header=None, delimiter=r'\s+', names=barocolnames) 
    BaroDF['Date'] = [dt.date(y,m,d) for y,m,d in zip(BaroDF.Year, BaroDF.Month, BaroDF.Day)] 
    #print(BaroDF)
    NumBaroCal = len(BaroDF.index)

    # cal with correct card and baro
    baroidx = [i for i,(x,y) in enumerate(zip(BaroDF.Baro, BaroDF.Card)) if (x==Barometer)&(y==Card)]
    if len(baroidx)==1:
        baroidx=baroidx[0]
    elif len(baroidx)>1:
        print('multiple cals with this baro and card?')
    else:
        print('no calabration for Baro = {0} and Card = {1}'.format(Barometer,Card))
        #print(BaroDF)
        #baroidx = int(input(' Choose a cal from the following by inputing the row index: '))
        baroidx=3

    PCoef0 = BaroDF.iloc[baroidx].PCoef0
    PCoef1 = BaroDF.iloc[baroidx].PCoef1
    PTCoef0 = BaroDF.iloc[baroidx].PTCoef0
    PTCoef1 = BaroDF.iloc[baroidx].PTCoef1
    if math.isnan(PTCoef0):
        PTCoef0=0
        PTCoef1=0
    BarometerCal = BaroDF.iloc[baroidx].Baro
    CardforBarometer = BaroDF.iloc[baroidx].Card
    BaroCalDate = BaroDF.iloc[baroidx].Date

#    open_file(BaroUnit, BaroFile, 'read')
#    for i in range(2):
#        for j in range(MaxBaro + 1):
#            for k in range(MaxCards + 1):
#                PCoef[i][j][k] = -1.0  # Missing value flag
#    skip_lines(BaroUnit, 1)
#    NumBaroCals = read_int(BaroUnit)
#    skip_lines(BaroUnit, 1)
#    print(13)
#    for N in range(1, NumBaroCals + 1):
#        PCalDate0, J, K, PCoefData, PTCoefData = read_pcal_data(BaroUnit)
#        PCalDate[J][K] = PCalDate0
#        for i in range(2):
#            PCoef[i][J][K] = PCoefData[i]
#            PTCoef[i][J][K] = PTCoefData[i]
#        print(14, J, K, PCalDate0)
#    BarometerToUse = Barometer
#    CardForBarometer = Card
#    J = Barometer
#    K = Card
#if PCoef[0][J][K] < 0.0:
#    print('No cals for Barometer, card=', J, K)
#    print('Select (Baro,Card) indices from above list')
#    J, K = read_input()
#    BarometerToUse = J
#    CardForBarometer = K
#    print('Cal date=', PCalDate[J][K])
#    print('Pcoef   =', [PCoef[i][J][K] for i in range(2)])
#    print('PTCoef  =', [PTCoef[i][J][K] for i in range(2)])
#    halt()
#else:
#    PCalDate[Barometer][Card] = Date(0, 0, 0)
#
if Anemometer > 0:
    #        Menu.FileName = os.path.join(RootPath, AnalogCalPath, 'wcal.men')
#        # call Menu_Edit(Menu)
#        Menu.Value = read_menu_value()
#        NumWindCals, WindCalData = parse_wcal_menu_value(Menu.Value)
#        for j, wcal_date, wcoef in WindCalData:
#            WCalDate[j] = wcal_date
#            WCoef[0][j], WCoef[1][j] = wcoef
      NumWindCals=5
      WCalDate[1]=tm.Date(1999, 3, 11)
      WCoef[1]=[0.289, 0.2699]
      WCalDate[2]=tm.Date(1999, 3, 11)
      WCoef[2]=[0.337, 0.2673]
      WCalDate[3]=tm.Date(2001, 1, 8)
      WCoef[3]=[0.441, 0.2496]
      WCalDate[4]=tm.Date(2001, 1, 8)
      WCoef[4]=[0.413, 0.2510]
      WCalDate[5]=tm.Date(2001, 1, 8)
      WCoef[5]=[0.413, 0.2510]  
else:
    WCalDate[Anemometer] = tm.Date(0, 0, 0)

if Card >= 0:
    #        Menu.FileName = os.path.join(RootPath, AnalogCalPath, 'vcal.men')
#        # call Menu_Edit(Menu)
#        Menu.Value = read_menu_value()
#        NumVoltCals, VoltCalData = parse_vcal_menu_value(Menu.Value)
#        for j, vcal_date, vcoef in VoltCalData:
#            VCalDate[j] = vcal_date
#            VCoef[0][j], VCoef[1][j] = vcoef
      VCalDate[1]=tm.Date(1999,7,28)   
      VCoef[1]=[0.305,0.015691]
      VCalDate[2]=tm.Date(1999,7,28)   
      VCoef[2]=[0.212,0.01593]
      VCalDate[3]=tm.Date(1999,7,28)   
      VCoef[3]=[0.166,0.015838]
      VCalDate[4]=tm.Date(1999,7,28)   
      VCoef[4]=[0.06,0.015879]
      VCalDate[5]=tm.Date(2000,8,10)   
      VCoef[5]=[-0.01,0.016]
      VCalDate[6]=tm.Date(2001,5,2)   
      VCoef[6]=[0.039,0.015786]
      VCalDate[10]=tm.Date(2018,6,7)   
      VCoef[10]=[0.29,0.017708]
      VCalDate[11]=tm.Date(2005,7,29)   
      VCoef[11]=[0.053,0.01863]
      VCalDate[12]=tm.Date(2005,7,29)   
      VCoef[12]=[0.072,0.01856]
      VCalDate[21]=tm.Date(2009,5,11)   
      VCoef[21]=[0.1466,0.0159]
      VCalDate[22]=tm.Date(2009,5,11)   
      VCoef[22]=[0.0461,0.0159]
      VCalDate[23]=tm.Date(2009,5,11)   
      VCoef[23]=[0.0461,0.0159]
      VCalDate[24]=tm.Date(2009,5,11)   
      VCoef[24]=[0.0461,0.0159]

OutRootPath = RootPath+'/PyOut/'+Station+'/' 
SumPath = OutRootPath + SumPath + '/#' + str(Instrument).zfill(2) +'/'
InRootPath = RootPath + DataPath + '/' + Station + '/#' + str(Instrument).zfill(2) + '/'

def write_env_header(filename, Station, TZ, Instrument, Barometer, Card, BaroCal, CardForBarometer,
                     BaroCalDate, Anemometer, WCalDate, LoggerVersion, VCalDate):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w") as f:
        f.write('# Station              {}\n'.format(Station))
        f.write('# Time offset          {}\n'.format(TZ))
        f.write('# Cimel                 {}\n'.format(Instrument))
        f.write('# Barometer, Card at site      {0}  {1}\n'.format(Barometer,Card))
        f.write('# Barometer, Card used in cal  {0}  {1}   Cal Date {2}\n'.format(BaroCal,CardforBarometer,BaroCalDate))
        f.write('# Anemometer            {0}            Cal Date    {1}\n'.format(Anemometer,WCalDate[Anemometer]))
        f.write('# LoggerVersion         {0}\n'.format(LoggerVersion))
        f.write('# Card                  {0}            Cal Date {1}\n'.format(Card,VCalDate[Card]))
        


for Year in range(YearMin, YearMax + 1):
    Year_Code = str(Year).zfill(4)
    MMin = 1
    MMax = 12
    if Year == YearMin:
        MMin = MonthMin
    if Year == YearMax:
        MMax = MonthMax
    for Month in range(MMin, MMax + 1):
        Month_Code = str(Month).zfill(2)
        Dmin = 1
        Dmax = tm.days_in_month(Month, Year)
        if Month == MonthMin:
            Dmin = DayMin
        if Month == MonthMax:
            Dmax = DayMax
        for Day in range(Dmin, Dmax + 1):
            Day_Code = str(Day).zfill(2)

            print('*************************************')
            InFilePath = os.path.join(InRootPath, '', Year_Code, '', Month_Code, '', Day_Code, '')
            ObsDate = dt.date(Year, Month, Day)

            if tm.datediff(ObsDate,StartDate) < 0:
                print(' Date preceeds first valid date of configuration!')
                print('Start date', StartDate)
                print('End   date', EndDate)
                print('Date      ', ObsDate)
                NumDays = ObsDate - StartDate
                print('Table start to date', NumDays)
                break
            elif tm.datediff(EndDate,ObsDate) < 0:
                print(' Date is after last valid date of configuration!')
                print('End   date', EndDate)
                print('Date      ', ObsDate)
                NumDays = EndDate - ObsDate
                print('Date to Table end', NumDays)
                break

            FileRoot = Year_Code[2:4] + Month_Code + Day_Code
            InFile = InFilePath + FileRoot + '.'+Extension
            print(' Input File=', InFile)
            #FileRoot = FileRoot.strip()
            HpaFile = OutRootPath + FileRoot +  '.hpa'
            print(' HpaFile=', HpaFile)
            write_env_header(HpaFile, Station, TZ, Instrument, Barometer, Card, BarometerCal, CardForBarometer,
                     BaroCalDate, Anemometer, WCalDate, LoggerVersion, VCalDate)
            with open(HpaFile, "a") as f:
                f.write('#YYYY MM dd hh mm ss P[raw]  Tmean   DelP P[comp]\n')

            if Extension == 'sta':
                pass
            elif Extension == 'rnv':
                EnvFile = OutRootPath + FileRoot + '.env'
                print(' EnvFile=', EnvFile)
                write_env_header(EnvFile, Station, TZ, Instrument, Barometer, Card,
                        BarometerCal, CardForBarometer, BaroCalDate, Anemometer,
                        WCalDate, LoggerVersion, VCalDate)
                with open(EnvFile, 'a') as f:
                    f.write('#YYYY MM dd hh mm ss            Pressure               Temperature               Wind speed (m/s)         Wind direction (E of N)   Battery voltage')
                    f.write('#                     Mean    Min    Max   Sdev   Mean    Min    Max   Sdev   Mean    Min    Max   Sdev   Mean    Min    Max   Sdev')
                Env = pd.read_csv(InFile, header=None, delimiter=r'\s+')
                Env.rename(columns={0:'Year', 1:'Month', 2:'Day', 3:'Hour', 4:'Minute', 5:'Second',
                                    6:'Pmean', 7:'Pmin', 8:'Pmax', 9:'Psd',
                                    10:'Tmean', 11:'Tmin', 12:'Tmax', 13:'Tsd',
                                    14:'WVmean', 15:'WVmin', 16:'WVmax', 17:'WVsd',
                                    18:'WDmean', 19:'WDmin', 20:'WDmax', 21:'WDsd',
                                    22:'VBatt'}, inplace=True)
                Env['DateTime'] = [dt.datetime(y,m,d,h,mi,s) for y,m,d,h,mi,s in zip(Env.Year, Env.Month, Env.Day, Env.Hour, Env.Minute, Env.Second)]
                J = Card
                Env['TOutMean'] = TCoef[J][0] + TCoef[J][1] * Env.Tmean
                #Temperature.Min = Temperature.Mean - TCoef[J][1] * Env.Tmin
                #Temperature.Max = Temperature.Mean + TCoef[J][1] * Env.Tmax
                #Temperature.Sdev = TCoef[J][0] * Env.Tsdev
                if Barometer > 0:
                    J = BarometerToUse
                    K = CardForBarometer
                    Env['POutMean'] = PCoef0 + PCoef1 * Env.Pmean
                    #Pressure.Min = Pressure.Mean - PCoef1 * Env.Pmin
                    #Pressure.Max = Pressure.Mean + PCoef1 * Env.Pmax
                    #Pressure.Sdev = PCoef1 * Env.Psdev

                    Env['Pressure_UnTComp'] = Env.POutMean
                    Env['DelPressure'] = PTCoef1 * (Env.TOutMean - PTCoef0)
                    Env['POutMean'] = Env.POutMean - Env.DelPressure
                    #Pressure.Min = Pressure.Min - DelPressure
                    #Pressure.Max = Pressure.Max - DelPressure
                else:
                    Env['POutMean'] = 0
                    #Pressure.Min = 0
                    #Pressure.Max = 0
                    #Pressure.Sdev = 0

                if Anemometer > 0:
                    J = Anemometer
                    Env['WindVOutMean'] = WCoef[J][0] + WCoef[J][1] * Env.WVmean / 100
                    #WindV.Min = WindV.Mean - WCoef[J][1] * Env.WVmin / 100
                    #WindV.Max = WindV.Mean + WCoef[J][1] * Env.WVmax / 100
                    #WindV.Sdev = WCoef[J][1] * Env.WVsdev / 100
                else:
                    Env['WindVOutMean'] = -1.0
                    #WindV.Min = -1.0
                    #WindV.Max = -1.0
                    #WindV.Sdev = -1.0

                if Anemometer > 0:
                    J = Anemometer
                    Env['WindDOutMean'] = Env.WDmean / 100
                    # Boom points west not north, so rotate axes by -90 (+270) degrees.
                    Env['WindDOutMean'] = (Env.WindDOutMean + 270.0) % 360.0
                    #WindD.Min = WindD.Mean - Env.WDmin / 100
                    #WindD.Max = WindD.Mean + Env.WDmax / 100
                    #WindD.Sdev = Env.WDsdev / 100
                else:
                    Env['WindDOutMean'] = 0
                    #WindD.Min = 0
                    #WindD.Max = 0
                    #WindD.Sdev = 0

                J = Card
                Env['VBattOut'] = VCoef[J][0] + VCoef[J][1] * Env.VBatt

                # Add Time zone offset (TZ=-10h for AEST)
                Env['DateTimeUTC'] = Env.DateTime + TZ
                # Add clock error.
                # Negative error: Fast clock
                # Positive error: Slow clock
                # NB This convention differs from that used before October 2000
                Env['DateTimeUTC'] = Env.DateTimeUTC + CE
                
                if Dbug[0]:
                    print('DateTimeLoc:', Env.DateTime)
                    print('DateTimeUTC:', Env.DateTimeUTC)

                #write_hpa_data(HpaUnit, DateTimeUTC, Pressure_UnTComp, Temperature.Mean,
                #        DelPressure, Pressure.Mean)
                with open(HpaFile, 'a') as f:
                    f.write(Env.to_string(header=False, index=False, columns=['DateTimeUTC', 'Pressure_UnTComp', 'TOutMean', 'DelPressure', 'POutMean']))
                    f.write('\n')

                if Dbug[0]:
                    print('DateTimeLoc:', Env.DateTime)
                    print('DateTimeUTC:', Env.DateTimeUTC)
                    print('Pressure:', Env.POutMean)

                    #write_env_data(EnvUnit, DateTimeUTC, Pressure, Temperature, WindV, WindD, VBatt)
                with open(EnvFile, 'a') as f:
                    f.write(Env.to_string(header=False, index=False, columns=['DateTimeUTC', 'POutMean', 'TOutMean', 'WindVOutMean',
                                                                              'WindDOutMean','VBattOut']))
                    f.write('\n')

                    #Env, Eof = read_raw_environment_record(Dbug[0], InUnit)

#            elif Extension == 'onv':
#                EnvFile = make_file_name(FilePath, FileRoot, 'env')
#                print(' EnvFile=', EnvFile)
#                EnvUnit = next_unit_number()
#                open_file(EnvUnit, EnvFile, 'write')
#                write_env_header(EnvUnit, Station, TZ, Instrument, Barometer, Card,
#                        BarometerToUse, CardForBarometer,
#                        PCalDate[BarometerToUse][CardForBarometer], Anemometer,
#                        WCalDate[Anemometer], LoggerVersion, Card, VCalDate[Card])
#                Env, Eof = read_old_environment_record(Dbug[0], InUnit)
#                while Eof == 0:
#                    J = Barometer
#                    K = Card
#                    Pressure.Mean = PCoef0 + PCoef1 * Env.Pmean  # untested
#                    Pressure.Min = Pressure.Mean - PCoef1 * Env.Pmin
#                    Pressure.Max = Pressure.Mean + PCoef1 * Env.Pmax
#                    Pressure.Sdev = PCoef1 * Env.Psdev
#                    J = Card
#                    Temperature.Mean = TCoef[0][J] + TCoef[1][J] * Env.Tmean
#                    Temperature.Min = Temperature.Mean - TCoef[1][J] * Env.Tmin
#                    Temperature.Max = Temperature.Mean + TCoef[1][J] * Env.Tmax
#                    Temperature.Sdev = TCoef[1][J] * Env.Tsdev
#                    J = Card
#                    VBatt = VCoef[J][0] + VCoef[J][1] * Env.VBatt
#
#                    # Add Time zone offset (TZ=-10h for AEST)
#                    DateTimeUTC = Env.DateTime + Date_Time(0, 0, 0, TZ.Hour, TZ.Minute, TZ.Second)
#
#                    # Add clock error.
#                    # Negative error: Fast clock
#                    # Positive error: Slow clock
#                    # NB This convention differs from that used before October 2000
#                    DateTimeUTC = DateTimeUTC + Date_Time(0, 0, 0, CE.Hour, CE.Minute, CE.Second)
#
#                    if Dbug[0]:
#                        print('DateTimeLoc:', Env.DateTime)
#                        print('DateTimeUTC:', DateTimeUTC)
#                        halt()
#
#                    write_hpa_data(HpaUnit, DateTimeUTC, Pressure_UnTComp, Temperature.Mean,
#                            DelPressure, Pressure.Mean)
#
#                    if Dbug[0]:
#                        print('DateTimeLoc:', Env.DateTime)
#                        print('DateTimeUTC:', DateTimeUTC)
#                        print('Pressure:', Env.Pressure)
#                        halt()
#
#                    WindV.Mean = -99.9
#                    WindV.Min = -99.9
#                    WindV.Max = -99.9
#                    WindV.Sdev = -99.9
#
#                    WindD.Mean = -99.9
#                    WindD.Min = -99.9
#                    WindD.Max = -99.9
#                    WindD.Sdev = -99.9
#
#                    write_env_data(EnvUnit, DateTimeUTC, Pressure, Temperature, WindV, WindD, VBatt)
#
#                    Env, Eof = read_old_environment_record(Dbug[0], InUnit)

            else:
                print('Extension ', Extension, ' not supported.')
                halt()


print('Done.')


