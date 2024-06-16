import numpy as np

LnV0Erms = np.zeros(NumChannels)
V0Erms = np.zeros(NumChannels)
WRefPriorGenCals = np.zeros(NumChannels)
USignal = np.zeros(NumChannels)
StdUAod1AM = np.zeros(NumChannels)

LnV0Coef = np.zeros((MaxPolyOrder + 1, NumChannels))

PData = np.zeros(MaxPressurePoints, dtype=object)
DateTime = np.zeros(MaxStack, dtype=object)
CutList = np.zeros(MaxCut, dtype=object)
WVTime = np.zeros(2, dtype=object)
WVDateTime = np.zeros(2, dtype=object)

SunData = np.zeros(MaxStack, dtype=object)
Env = np.zeros(MaxEnvRecs + 1, dtype=object)

IWavelength = np.zeros(NumChannels, dtype=int)
Iorder = np.zeros(NumChannels, dtype=int)
Isort = np.zeros(NumChannels, dtype=int)

JulDay = np.zeros(MaxStack, dtype=int)
YMDHMS = np.zeros((7, MaxStack), dtype=int)

TimeDay = np.zeros(MaxStack)
TimeUT = np.zeros(MaxStack)
TimeAES = np.zeros(MaxStack)
AESTMin = np.zeros(MaxStack)
DSunEarth = np.zeros(MaxStack)
SolarZenithDeg = np.zeros(MaxStack)
SolarZenithApp = np.zeros(MaxStack)
SolarAzimuthDeg = np.zeros(MaxStack)
Pressure = np.zeros(MaxStack)

AirMass = np.zeros((MaxStack, NumChannels))
Volt = np.zeros((MaxStack, NumChannels))
VoltLog = np.zeros((MaxStack, NumChannels))

Aod1 = np.zeros((MaxObsDay, NumChannels))
Aod2 = np.zeros((MaxObsDay, NumChannels))
StdUAod = np.zeros((MaxObsDay, NumChannels))
Tran = np.zeros((MaxObsDay, NumChannels))
WvapOD = np.zeros((MaxObsDay, NumChannels))
SignalCV = np.zeros((MaxObsDay, NumChannels))
AodDay = np.zeros((MaxObsDay, NumChannels))

DAlphaDTau = np.zeros(MaxObsDay)
RCorr = np.zeros(MaxObsDay)

AodMonth = np.zeros((MaxObsMonth, NumChannels))
WvapMonth = np.zeros(MaxObsMonth)
Angstrom48Month = np.zeros(MaxObsMonth)

Angstrom46 = np.zeros(MaxObsDay)
StdUAngstrom46 = np.zeros(MaxObsDay)
Angstrom48 = np.zeros(MaxObsDay)
StdUAngstrom48 = np.zeros(MaxObsDay)
Angstrom68 = np.zeros(MaxObsDay)
StdUAngstrom68 = np.zeros(MaxObsDay)
Angstrom48Day = np.zeros(MaxObsDay)
WvapDay = np.zeros(MaxObsDay)
Wvap = np.zeros(MaxObsDay)
TDet = np.zeros(MaxObsDay)
Sample = np.zeros(MaxObsDay)
CirrusFlag = np.zeros(MaxObsDay + 1, dtype=int)

AodDayMean = np.zeros(NumChannels)
AodDaySdev = np.zeros(NumChannels)
AodDayMin = np.zeros(NumChannels)
AodDayMax = np.zeros(NumChannels)
AodMonthMean = np.zeros(NumChannels)
AodMonthSdev = np.zeros(NumChannels)
AodMonthMin = np.zeros(NumChannels)
AodMonthMax = np.zeros(NumChannels)
AodTemp = np.zeros(NumChannels)

AodGroupMean = np.zeros((MaxGroups, NumChannels))
AodGroupSdev = np.zeros((MaxGroups, NumChannels))

WaveCal = np.zeros((NumChannels, NumModels))
ExtraOD = np.zeros(NumChannels)
VoltMean = np.zeros(NumChannels)
VoltSdev = np.zeros(NumChannels)
VoltCv = np.zeros(NumChannels)
BlackSun = np.zeros(NumChannels)
BlackSkyA = np.zeros(NumSkyA)
BlackSkyK = np.zeros(NumSkyK)
RayleighOD = np.zeros(NumChannels)
OzoneOD = np.zeros(NumChannels)
LnV01AU = np.zeros(NumChannels)
LnV0 = np.zeros(NumChannels)
V01AU = np.zeros(NumChannels)

DtEnv = np.zeros(MaxEnvRecs + 1, dtype=int)

TranBoM = np.zeros((MaxStack, NumChannels))

Ozone = np.zeros(32)
Aureole = np.zeros(NumChannels)




ExtraOD[:] = 0
V01AU[:] = -1
V0Erms[:] = -1

NumArgs, LengthArgs, CommandLineArgs = Get_Arguments()

if NumArgs == 1:
    print('For batch mode: aod2p b xxx m1 m2')
    print('xxx = site code (e.g., la3)')
    print('m1  = number of first menu')
    print('m2  = number of last  menu')
    exit()
elif NumArgs == 0 or NumArgs > 4:
    BatchMode = False
    MenuOpCode = ''
    MenuName[0] = 'aod2p.men'
    FirstInterval = 1
    LastInterval = 1
elif CommandLineArgs[0].strip() == 'b':
    BatchMode = True
    MenuOpCode = '!'
    MenuSubPath = './men/' + CommandLineArgs[1].strip()
    FirstInterval = String_to_Integer(CommandLineArgs[2])
    LastInterval = String_to_Integer(CommandLineArgs[3])
    for Interval in range(FirstInterval, LastInterval + 1):
        MenuName[Interval - 1] = MenuSubPath + '/aod2p_' + str(Interval).zfill(2) + '.men'
        print('First=', FirstInterval)
        print('Last =', LastInterval)
        print('menu name=', MenuName[Interval - 1])
        Halt()

for Interval in range(FirstInterval, LastInterval + 1):
    print(' Processing interval #', Interval)
    print(' Menu file name:      ', MenuName[Interval - 1])
    Halt()
    Menu = {'FileName': MenuName[Interval - 1]}
    Menu = Menu_Edit(Menu, Operation=MenuOpCode)
    Dbug[:4], MakeTran, CirrusOption, Aod_500_Coarse_Threshold, ApplyClockFix, NumStack, Window, CvMax, MinObs, AodDaySdevCrit, Rel_Uaod_440_Threshold, Uangst48_Threshold, Extension, SiteCode, DataRoot, DataSuffix, CalPath, CalName, ResultPathA, ResultPathB, ConfigRoot, YearMin, MonthMin, DayMin, YearMax, MonthMax, DayMax, Epoch_Year, Epoch_Month, Epoch_Day, PressureOption, DefaultSurfacePressure, BoMPresFile, OzoneOption, DefaultOzone, OzoneFile = Menu['Value']

    MenuA = Menu

    Menu = Menu_Edit(Menu, Operation='&')
    WVWindowFlag, WVTime[0], WVTime[1], NWin = Menu['Value']

    MenuB = Menu

    ConfigPath = DataRoot.strip() + ConfigRoot.strip() + '/'

    ConfigFile = Make_File_Name(ConfigPath, SiteCode, 'cfn')
    print(' Config File=', ConfigFile)
    ConfigUnit = Next_Unit_Number()
    with open(ConfigFile, 'r') as file:
        Config = Read_Site_Configuration(Dbug[0], file, Date(YearMin, MonthMin, DayMin))

    Instrument = Config['CimelNumber']
    Model = Config['CimelModel']
    StationLatDeg = Config['Latitude']
    StationLonDeg = Config['Longitude']
    StationLatRad = StationLatDeg * RadiansPerDegree
    StationLonRad = StationLonDeg * RadiansPerDegree
    TZ_hms = Hour2HMS(Config['TimeZone'])
    TimeZone_DT = Date_Time(0, 0, 0, -TZ_hms['Hour'], -TZ_hms['Minute'], -TZ_hms['Second'])
    OzoneFile = DataRoot.strip() + '/ozone/' + OzoneFile.strip()
    BoMPresFile = DataRoot.strip() + '/pressure/' + BoMPresFile.strip()
    ResultPath = DataRoot.strip() + ResultPathA.strip() + ResultPathB.strip()
    CodPath = DataRoot.strip() + ResultPathA.strip() + '/dw2/cod/2005/'
    DataType = Extension
    IWavelength[:] = np.array([int(Wavelength[i, Model]) for i in range(NumChannels)])

    print('OzoneFile=', OzoneFile)
    FileOK = os.path.isfile(OzoneFile)
    if FileOK:
        OzoneUnit = Next_Unit_Number()
        with open(OzoneFile, 'r') as file:
            pass
    else:
        print('Error opening ozone data file', OzoneFile)
        Halt()

    if PressureOption == 3:
        print('PressureFile=', BoMPresFile)
        FileOK = os.path.isfile(BoMPresFile)
        if FileOK:
            BoMPresUnit = Next_Unit_Number()
            with open(BoMPresFile, 'r') as file:
                CheckFileDosUnix(file)
        else:
            print('Error opening Pressure data file', BoMPresFile)
            Halt()

    CalPath = DataRoot.strip() + CalPath.strip() + '/' + '#' + str(Instrument).zfill(2) + '/'
    CalFile = CalPath + CalName.strip()

    if CalName.strip() != Config['CimCalFile'].strip():
        print('Calibration File from menu=', CalName)
        print('Calibration File from config=', Config['CimCalFile'].strip())
        print(' Cal file specified in menu disagrees with')
        print(' Cal file specified in configuration file.')
        print(' Cntl C to abort, enter to continue')
        Halt()

    There = os.path.isfile(CalFile)
    CalUnit = Next_Unit_Number()
    if There:
        with open(CalFile, 'r') as file:
            pass
    else:
        print(' Can''t find calibration file ', CalFile)
        Halt()

    file.readline()
    file.readline()
    file.readline()
    CalInstrument = int(file.readline())
    CalModel = int(file.readline())
    CalLangleyChannels = int(file.readline())
    CalNumLangleys = int(file.readline())
    CalNumGeneralCycles = int(file.readline())
    CalEpoch = [int(x) for x in file.readline().split()]
    file.readline()
    for N in range(CalLangleyChannels):
        line = file.readline().split()
        WaveCal[N, Model] = float(line[0])
        Iorder[N] = int(line[1])
        LnV0Coef[0:Iorder[N]+1, N] = [float(x) for x in line[2:Iorder[N]+3]]
        LnV0Erms[N] = float(line[Iorder[N]+3])
        WRefPriorGenCals[:CalNumGeneralCycles] = [float(x) for x in line[Iorder[N]+4:]]
    line = file.readline().split()
    WaveCal[IWV, Model] = float(line[0])
    Iorder[IWV] = int(line[1])
    LnV0Coef[0:Iorder[IWV]+1, IWV] = [float(x) for x in line[2:Iorder[IWV]+3]]
    LnV0Erms[IWV] = float(line[Iorder[IWV]+3])

    file.close()

    for N in range(CalLangleyChannels):
        USignal[N] = 0.003
        if CalNumLangleys > 0:
            StdUAod1AM[N] = np.sqrt(USignal[N]**2 + LnV0Erms[N]**2 / CalNumLangleys)
        else:
            StdUAod1AM[N] = USignal[N]

    CutFile = Make_File_Name(CalPath, CalName, 'xod')
    print('Cut File Name=', CutFile)
    CutUnit = Next_Unit_Number()
    There = os.path.isfile(CutFile)
    if There:
        with open(CutFile, 'r') as file:
            Skip_Comments(file)
            NumCut = 1
            CutList[NumCut - 1], Eof = file.readline().split()
            while Eof == 0:
                NumCut += 1
                CutList[NumCut - 1], Eof = file.readline().split()
            NumCut -= 1
            print('Number of cut periods=', NumCut)
            Halt()
            print('Num days cut=', NumCut)
    else:
        print('No valid .cut file found')
        NumCut = 0

    if ApplyClockFix:
        ClockFile = Make_File_Name(CalPath, CalName, 'clk')
        print('Clock File Name=', ClockFile)
        ClockUnit = Next_Unit_Number()
        ClockFixFileOK, MaxRecLength = os.path.isfile(ClockFile), 0
        if ClockFixFileOK:
            with open(ClockFile, 'r') as file:
                CheckFileDosUnix(file)
        else:
            print('Can''t apply clockfix because no .clk file found')
            exit()

    FileRoot = '#' + str(Instrument).zfill(1) + str(YearMin % 100).zfill(2) + str(MonthMin).zfill(2) + str(MonthMax).zfill(2)
    Ext4 = 'par'

    ParFile = Make_File_Name(ResultPath, FileRoot, Ext4)
    print('par File Name=', ParFile)
    ParUnit = Next_Unit_Number()
    with open(ParFile, 'w') as file:
        RunDateTime = datetime.datetime.now()
        RunDate = RunDateTime.strftime('%Y-%m-%d')
        RunTime = RunDateTime.strftime('%H:%M:%S')
        file.write('*** Aod2p run on ' + RunDate + ' at ' + RunTime + ' ***\n')
    file.close()

    MenuA['Filename'] = ParFile
    MenuB['Filename'] = ParFile
    Menu_Edit(MenuA, Operation='+')
    Menu_Edit(MenuB, Operation='+')

FileRoot = '#' + str(Instrument).zfill(1) + str(YearMin % 100).zfill(2) + str(MonthMin).zfill(2) + str(MonthMax).zfill(2)
Ext4 = 'day'

if CirrusOption != 0:
    Ext4 += str(CirrusOption)

DayFile = Make_File_Name(ResultPath, FileRoot, Ext4)
print('Daily aod File Name=', DayFile)
There = os.path.isfile(DayFile)
DayUnit = Next_Unit_Number()
with open(DayFile, 'w') as file:
    file.write(FormatString[0].format(NumLangleyChannels(Model)))
    file.write(SiteCode, StationLatDeg, StationLonDeg, Instrument, Model, CalFile, OzoneFile, CirrusOption,
               *[IWavelength[I] for I in range(NumLangleyChannels(Model)) for _ in range(4)],
               IWavelength[IWV], IWavelength[IWV], IWavelength[IWV], IWavelength[IWV],
               4487, 4487, 4487, 4487)

Ext4 = 'mon'
if CirrusOption != 0:
    Ext4 += str(CirrusOption)

MonthFile = Make_File_Name(ResultPath, FileRoot, Ext4)
print('Monthly aod File Name=', MonthFile)
MonthUnit = Next_Unit_Number()
with open(MonthFile, 'w') as file:
    file.write(FormatString[1].format(NumLangleyChannels(Model)))
    file.write(SiteCode, StationLatDeg, StationLonDeg, Instrument, Model, CalFile, OzoneFile, CirrusOption,
               Epoch['Year'], Epoch['Month'], Epoch['Day'],
               *[IWavelength[I] for I in range(NumLangleyChannels(Model)) for _ in range(4)],
               IWavelength[IWV], IWavelength[IWV], IWavelength[IWV], IWavelength[IWV],
               4487, 4487, 4487, 4487)

Ext4 = 'xmy'
XMonthFile = Make_File_Name(ResultPath, FileRoot, Ext4)
print('Exclude monthly File Name=', XMonthFile)
XMonthUnit = Next_Unit_Number()
with open(XMonthFile, 'w') as file:
    file.write(SiteCode + '\n')

Ext4 = 'tex'
if CirrusOption != 0:
    Ext4 += str(CirrusOption)

TeXFile = Make_File_Name(ResultPath, FileRoot, Ext4)
N = NumLangleyChannels(Model)
Mask[:] = False
Mask[:N] = True
Mask[IWV] = False
for I in range(N):
    I1 = np.argmin(IWavelength[Mask])
    Isort[I] = I1
    Mask[Isort[I]] = False

print('Monthly aod File Name=', TeXFile)
TeXUnit = Next_Unit_Number()
with open(TeXFile, 'w') as file:
    file.write(FormatString[2].format(NumLangleyChannels(Model)))
    file.write(*[IWavelength[Isort[I]] for I in range(NumLangleyChannels(Model))])

if CirrusOption == 5:
    Ext4 = 'laod'
    CirrusAodFile = Make_File_Name(ResultPath, FileRoot, Ext4)
    print('CirrusAod File Name=', CirrusAodFile)
    CirrusAodUnit = Next_Unit_Number()
    with open(CirrusAodFile, 'w') as file:
        file.write(FormatString[3])

if CirrusOption >= 4:
    MPLFile = '/data/aerosol/DarwinARCS/lidar_cdf/CloudScreen_2005_whole_year.txt'
    print('MPL file=', MPLFile)
    MPLUnit = Next_Unit_Number()
    with open(MPLFile, 'r') as file:
        Skip(file, 9)

NumDaysOK = 0
NObsMonth = 0
NWvMonth = 0
CurrentMonth = 0

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
        DMin = 1
        DMax = Days_in_Month(Month, Year)

        if Month == MonthMin and Year == YearMin:
            DMin = DayMin
        if Month == MonthMax and Year == YearMax:
            DMax = DayMax

        for Day in range(DMin, DMax + 1):
            Day_Code = str(Day).zfill(2)

            print('*************************************')
            print('Instrument #', Instrument, Year_Code + Month_Code + Day_Code)

            DataPath = DataRoot.strip() + DataSuffix.strip() + '/' + Config['Id'] + '/' + '#' + str(Instrument).zfill(2) + '/' + Year_Code + '/' + Month_Code + '/' + Day_Code + '/'

            DateObs = [Year, Month, Day]
            NDay = DayCount(CalEpoch, DateObs)
            for N in range(CalLangleyChannels):
                LnV01AU[N] = sum(LnV0Coef[I, N] * NDay**I for I in range(Iorder[N] + 1))
                V01AU[N] = np.exp(LnV01AU[N])
                V0Erms[N] = V01AU[N] * LnV0Erms[N]
            N = IWV
            LnV01AU[N] = sum(LnV0Coef[I, N] * NDay**I for I in range(Iorder[N] + 1))
            V01AU[N] = np.exp(LnV01AU[N])
            V0Erms[N] = V01AU[N] * LnV0Erms[N]

            ObsDate = Date(Year, Month, Day)

            if Month != CurrentMonth:
                ThisMonthsDayStamp = Date(Year, Month, 1) - Epoch
                Ozone = Read_BoM_Ozone_2011(Dbug[0], OzoneUnit, ObsDate)
                CurrentMonth = Month

            FileRoot = Year_Code[2:] + Month_Code + Day_Code

            InFile = Make_File_Name(DataPath, FileRoot, Extension)
            print('InFile=', InFile)
            FileOK = os.path.isfile(InFile)

            Include = True
            NObsOK = -1

            ObsDate = Date(Year, Month, Day)
            ObsXDateTime = Date_Time(Year, Month, Day, -1, 0, 0)
            Include = True

            for Ic in range(NumCut):
                if ObsXDateTime == CutList[Ic]:
                    Include = False
                    print("Whole day excluded:", ObsDate)
                    Halt()
                    break

            if FileOK and Include:
                InUnit = Next_Unit_Number()
                with open(InFile, 'r') as file:
                    OutRoot = '#' + str(Instrument).zfill(1) + FileRoot
                    LogFile = Make_File_Name(ResultPath, OutRoot, 'log')
                    print('LogFile=', LogFile)
                    LogUnit = Next_Unit_Number()
                    with open(LogFile, 'w') as log_file:
                        log_file.write(FormatString[4].format(SiteCode, *[Wavelength[I, Model] for I in range(NumLangleyChannels(Model))], Wavelength[IWV, Model]))

                        if PressureOption == 1:
                            PresPath = DataPath
                            if Instrument == 2 and Year == 1998 and 4 <= Month <= 7:
                                I = DataPath.index('tta')
                                PresPath = DataPath[:I] + 'ttb/#03' + DataPath[I+7:]
                            PresFile = Make_File_Name(PresPath, FileRoot, 'hpa')
                            print('Pressure File Name=', PresFile)
                            FileOK = os.path.isfile(PresFile)
                            if FileOK:
                                NumPressurePoints, NumValidPresPoints, PData, DailyMeanPressure = Read_Pressure_File(Dbug[0], PresFile)
                                if NumValidPresPoints <= MinValidPresPoints:
                                    print('Insufficient valid pressure data in file')
                            else:
                                NumPressurePoints = 0
                                print('Pressure file not found (.hpa)-has envcal been run?')
                                Halt()
                            if NumPressurePoints <= 0:
                                print('Assuming Pdefault =', DefaultSurfacePressure)

                        elif PressureOption == 3:
                            PData, NumPressurePoints, NumValidPresPoints, DailyMeanPressure = Read_BoM_Pressure(Dbug[0], BoMPresUnit, ObsDate)

                        if ApplyClockFix:
                            TimeCorr, Eof = Get_Time_Correction(Dbug[0], ClockUnit, ObsDate)

                        EnvFile = Make_File_Name(DataPath, FileRoot, 'env')
                        print('EnvFile=', EnvFile)
                        EnvFound = os.path.isfile(EnvFile)

                        if EnvFound:
                            EnvUnit = Next_Unit_Number()
                            with open(EnvFile, 'r') as env_file:
                                Skip_Comments(env_file)
                                I = 0
                                I += 1
                                if I > MaxEnvRecs:
                                    print('Array overflow impending, env')
                                    print('I=', I, 'Array Size=', MaxEnvRecs)
                                    Halt()
                                Env[I], Eof = Read_Environment_Record(Dbug[0], env_file)
                                while Eof == 0:
                                    I += 1
                                    if I > MaxEnvRecs:
                                        print('Array overflow impending, env')
                                        print('I=', I, 'Array Size=', MaxEnvRecs)
                                        Halt()
                                    Env[I], Eof = Read_Environment_Record(Dbug[0], env_file)
                                EnvRecs = I - 1
                                env_file.close()

                        BlkFile = Make_File_Name(DataPath, FileRoot, 'blk')
                        print('BlackFile=', BlkFile)
                        FileOK = os.path.isfile(BlkFile)

                        if FileOK:
                            BlkUnit = Next_Unit_Number()
                            with open(BlkFile, 'r') as blk_file:
                                I = 0
                                Eof = 0
                                BlackSun[:] = 0
                                BlackSkyA[:] = 0
                                BlackSkyK[:] = 0
                                while Eof == 0:
                                    DateTime[0], Black, Eof = Read_Black_Record(Dbug[0], Model, blk_file)
                                    if Eof == 0:
                                        I += 1
                                        BlackSun[:NumChannels] += Black['Sun']['Chan'][:NumChannels]
                                        BlackSkyA[:NumSkyA] += Black['SkyA'][:NumSkyA]
                                        BlackSkyK[:NumSkyK] += Black['SkyK'][:NumSkyK]
                                BlackSun[:NumChannels] /= I
                                BlackSkyA[:NumSkyA] /= I
                                BlackSkyK[:NumSkyK] /= I
                                if Dbug[1]:
                                    print('Number of Black records:', I)
                                    for J in range(NumChannels):
                                        print('Mean Black Sun     :', BlackSun[J])
                                    Halt()
                                    for J in range(NumSkyA):
                                        print('Mean Black SkyA     :', BlackSkyA[J])
                                    Halt()
                                    for J in range(NumSkyK):
                                        print('Mean Black SkyK     :', BlackSkyK[J])
                                    Halt()
                                blk_file.close()
                        else:
                            BlackSun[:NumChannels] = 0
                            BlackSkyA[:NumSkyA] = 0
                            BlackSkyK[:NumSkyK] = 0

                        if CirrusOption >= 3:
                            CodFile = Make_File_Name(CodPath, OutRoot, 'cod')
                            print('CodFile=', CodFile)
                            CodUnit = Next_Unit_Number()
                            with open(CodFile, 'w') as cod_file:
                                pass

                        if OzoneOption == 0:
                            OzoneColumn = DefaultOzone
                        elif OzoneOption == 1:
                            if Ozone[Day] > 0:
                                OzoneColumn = Ozone[Day] / 1000.0
                            else:
                                print('Ozone file: Missing day, using monthly mean=', Ozone[0])
                                OzoneColumn = Ozone[0] / 1000.0
                        elif OzoneOption == 2:
                            OzoneColumn = Ozone[0] / 1000.0

                        OzoneOD[:NumChannels] = OzoneCoef[:NumChannels, Model] * OzoneColumn

                        FileRoot = str(Config['BoMNumber']).zfill(3) + Year_Code + Month_Code + Day_Code
                        Ext4 = 'bom'
                        BoMFile = Make_File_Name(ResultPath, FileRoot, Ext4)
                        StationNumber = Config['BoMNumber']
                        TypeScatComp = 0
                        TimeZone = Config['TimeZone']
                        YMDHMS[0, 0] = Year
                        YMDHMS[1, 0] = Month
                        YMDHMS[2, 0] = Day
                        YMDHMS[3, 0] = 0
                        YMDHMS[4, 0] = 0
                        YMDHMS[5, 0] = 0
                        YMDHMS[6, 0] = 0

                        JulDay[0], TimeDay[0] = Julian(YMDHMS[:, 0])
                        SunPos = Solar_Position_Almanac(JulDay[0], TimeDay[0])
                        EarthSunRsqFactor = 1 / SunPos['DsunEarth']**2

                        print('BoM transmission File Name=', BoMFile)
                        There = os.path.isfile(BoMFile)
                        BoMUnit = Next_Unit_Number()
                        with open(BoMFile, 'w') as bom_file:
                            bom_file.write(FormatString[5].format(
                                Config['Name'], CR,
                                StationNumber, CR,
                                Year_Code, Month_Code, Day_Code, Day_of_Year(Day, Month, Year), CR,
                                -TimeZone, CR,
                                TypeScatComp, CR,
                                OzoneColumn, CR,
                                EarthSunRsqFactor, CR,
                                StationLatDeg, CR,
                                StationLonDeg, CR,
                                Config['Height'], CR,
                                NumLangleyChannels(Model) + 1, CR
                            ))

                            DateTime[0] = Date_Time(Year, Month, Day, 0, 0, 0)
                            Pressure[0] = GetPressure(DateTime[0], PData,
                                                      PressureOption,
                                                      NumPressurePoints,
                                                      NumValidPresPoints,
                                                      DailyMeanPressure,
                                                      DefaultSurfacePressure)

                            for N in range(NumLangleyChannels(Model)):
                                RayleighOD[N] = Rayleigh(Dbug[2], Wavelength[N, Model], Pressure[0])
                                if Dbug[2]:
                                    print('Wavelength:', Wavelength[N, Model])
                                    print('Pressure  :', Pressure[0])
                                    print('RayleighOD:', RayleighOD[N])
                                bom_file.write(FormatString[6].format(
                                    N, Wavelength[N, Model], RayleighOD[N], OzoneCoef[N, Model],
                                    V01AU[N], 2 * V0Erms[N], 1.0, CR
                                ))

                            N = IWV
                            J = NumLangleyChannels(Model) + 1
                            RayleighOD[N] = Rayleigh(Dbug[2], Wavelength[N, Model], Pressure[0])
                            bom_file.write(FormatString[6].format(
                                J, Wavelength[N, Model], RayleighOD[N], OzoneCoef[N, Model],
                                V01AU[N], 2 * V0Erms[N], 1.0, CR
                            ))

                            bom_file.write(FormatString[7].format(NumLangleyChannels(Model) + 1))
                            bom_file.write(FormatString[8].format(
                                *[IWavelength[I] for I in range(NumLangleyChannels(Model))],
                                IWavelength[IWV], IWavelength[IWV], CR
                            ))

                            bom_file.write(FormatString[9].format(NumLangleyChannels(Model) + 1))

                            Eof = 0
                            print('Processing', DataType, 'records...')
                            print('Pass 1 filtering...')
                            IGroup = 0
                            IGroup1 = 0
                            Is = 0

                            while True:
                                for I in range(NumStack):
                                    DateTime[I], SunData[I], ValidData, Eof = Read_Single_Sun_Record(
                                        Dbug[1], InUnit, Model
                                    )
                                    YMDHMS[:, I] = [
                                        DateTime[I].year,
                                        DateTime[I].month,
                                        DateTime[I].day,
                                        DateTime[I].hour,
                                        DateTime[I].minute,
                                        DateTime[I].second,
                                        0
                                    ]

                                    if ApplyClockFix:
                                        Add_Time(YMDHMS[0, I], -TimeCorr)

                                    JulDay[I], TimeDay[I] = Julian(YMDHMS[:, I])
                                    SunPos = Solar_Position_Almanac(JulDay[I], TimeDay[I])
                                    DSunEarth[I] = SunPos['DsunEarth']
                                    TimeUT[I] = TimeDay[I] - 0.5
                                    TimeAES[I] = TimeUT[I] * 24 + 10.0
                                    AestMin[I] = TimeAES[I] * 60

                                    if Dbug[2]:
                                        print('Year                :', YMDHMS[0, I])
                                        print('Month               :', YMDHMS[1, I])
                                        print('Day                 :', YMDHMS[2, I])
                                        print('Hour                :', YMDHMS[3, I])
                                        print('Minute              :', YMDHMS[4, I])
                                        print('Second              :', YMDHMS[5, I])
                                        print('Julian day number   :', JulDay[I])
                                        print('TimeDay             :', TimeDay[I])
                                        print('Sun Earth distance  :', DSunEarth[I])
                                        print('Time(UT)            :', TimeUT[I])
                                        print('Time(AES)           :', TimeAES[I])
                                        print('AsetMin             :', AestMin[I])
                                        Halt()

                                    LSTCut = Date_Time(Year, Month, Day, 3, 0, 0)
                                    TimeZoneHMS = Hour2HMS(TimeZone)
                                    UTCut = AddDateTime(LSTCut, Date_Time(
                                        0, 0, 0, -TimeZoneHMS['Hour'], -TimeZoneHMS['Minute'], -TimeZoneHMS['Second']
                                    ))
                                    if DateTime[I] < UTCut:
                                        ValidData = 14

                                    for Ic in range(NumCut):
                                        if DateTime[I] == CutList[Ic]:
                                            ValidData = 13

                                    NewStation = True
                                    SolarZenith, SolarAzimuth = Satellite_Position(
                                        JulDay[I], TimeDay[I],
                                        NewStation,
                                        StationLatRad, StationLonRad,
                                        SunElements
                                    )
                                    SolarZenithDeg[I] = SolarZenith * DegreesPerRadian
                                    SolarZenithApp[I] = Apparent_Zenith(SolarZenithDeg[I])
                                    SolarAzimuthDeg[I] = SolarAzimuth * DegreesPerRadian

                                    if Dbug[2]:
                                        print('Julian day    :', JulDay[I])
                                        print('TimeDay       :', TimeDay[I])
                                        print('Latitude      :', StationLatDeg)
                                        print('Longitude     :', StationLonDeg)
                                        print('Sol Zen (true):', SolarZenithDeg[I])
                                        print('Sol Zen (app) :', SolarZenithApp[I])
                                        print('Solar Azimuth :', SolarAzimuthDeg[I])
                                        Halt()

                                    RejectCode[I] = ValidData
                                    while ValidData != 0 and Eof == 0:
                                        log_file.write(FormatString[10].format(
                                            *YMDHMS[:6, I],
                                            SolarZenithApp[I],
                                            IGroup,
                                            RejectCode[I],
                                            *[SunData[I]['Signal']['Chan'][J] for J in range(NumLangleyChannels(Model))],
                                            SunData[I]['Signal']['Chan'][IWV]
                                        ))
                                        DateTime[I], SunData[I], ValidData, Eof = Read_Single_Sun_Record(
                                            Dbug[1], InUnit, Model
                                        )
                                        YMDHMS[:, I] = [
                                            DateTime[I].year,
                                            DateTime[I].month,
                                            DateTime[I].day,
                                            DateTime[I].hour,
                                            DateTime[I].minute,
                                            DateTime[I].second,
                                            0
                                        ]

                                        if ApplyClockFix:
                                            Add_Time(YMDHMS[0, I], -TimeCorr)

                                        JulDay[I], TimeDay[I] = Julian(YMDHMS[:, I])
                                        SunPos = Solar_Position_Almanac(JulDay[I], TimeDay[I])
                                        DSunEarth[I] = SunPos['DsunEarth']
                                        TimeUT[I] = TimeDay[I] - 0.5
                                        TimeAES[I] = TimeUT[I] * 24 + 10.0
                                        AestMin[I] = TimeAES[I] * 60

                                        if Dbug[2]:
                                            print('Year                :', YMDHMS[0, I])
                                            print('Month               :', YMDHMS[1, I])
                                            print('Day                 :', YMDHMS[2, I])
                                            print('Hour                :', YMDHMS[3, I])
                                            print('Minute              :', YMDHMS[4, I])
                                            print('Second              :', YMDHMS[5, I])
                                            print('Julian day number   :', JulDay[I])
                                            print('TimeDay             :', TimeDay[I])
                                            print('Sun Earth distance  :', DSunEarth[I])
                                            print('Time(UT)            :', TimeUT[I])
                                            print('Time(AES)           :', TimeAES[I])
                                            print('AsetMin             :', AestMin[I])
                                            Halt()

                                        LSTCut = Date_Time(Year, Month, Day, 3, 0, 0)
                                        TimeZoneHMS = Hour2HMS(TimeZone)
                                        UTCut = AddDateTime(LSTCut, Date_Time(
                                            0, 0, 0, -TimeZoneHMS['Hour'], -TimeZoneHMS['Minute'], -TimeZoneHMS['Second']
                                        ))
                                        if DateTime[I] < UTCut:
                                            ValidData = 14

                                        for Ic in range(NumCut):
                                            if DateTime[I] == CutList[Ic]:
                                                ValidData = 13

                                        NewStation = True
                                        SolarZenith, SolarAzimuth = Satellite_Position(
                                            JulDay[I], TimeDay[I],
                                            NewStation,
                                            StationLatRad, StationLonRad,
                                            SunElements
                                        )
                                        SolarZenithDeg[I] = SolarZenith * DegreesPerRadian
                                        SolarZenithApp[I] = Apparent_Zenith(SolarZenithDeg[I])
                                        SolarAzimuthDeg[I] = SolarAzimuth * DegreesPerRadian

                                        if Dbug[2]:
                                            print('Julian day    :', JulDay[I])
                                            print('TimeDay       :', TimeDay[I])
                                            print('Latitude      :', StationLatDeg)
                                            print('Longitude     :', StationLonDeg)
                                            print('Sol Zen (true):', SolarZenithDeg[I])
                                            print('Sol Zen (app) :', SolarZenithApp[I])
                                            print('Solar Azimuth :', SolarAzimuthDeg[I])
                                            Halt()

                                        RejectCode[I] = ValidData

                                    if Eof != 0:
                                        break

                                    for I in range(NumStack):
                                        YMDHMS[:, I] = [
                                            DateTime[I].year,
                                            DateTime[I].month,
                                            DateTime[I].day,
                                            DateTime[I].hour,
                                            DateTime[I].minute,
                                            DateTime[I].second,
                                            0
                                        ]

                                        if ApplyClockFix:
                                            Add_Time(YMDHMS[0, I], -TimeCorr)

                                        JulDay[I], TimeDay[I] = Julian(YMDHMS[:, I])
                                        SunPos = Solar_Position_Almanac(JulDay[I], TimeDay[I])
                                        DSunEarth[I] = SunPos['DsunEarth']
                                        TimeUT[I] = TimeDay[I] - 0.5
                                        TimeAES[I] = TimeUT[I] * 24 + 10.0
                                        AestMin[I] = TimeAES[I] * 60

                                        if Dbug[2]:
                                            print('Year                :', YMDHMS[0, I])
                                            print('Month               :', YMDHMS[1, I])
                                            print('Day                 :', YMDHMS[2, I])
                                            print('Hour                :', YMDHMS[3, I])
                                            print('Minute              :', YMDHMS[4, I])
                                            print('Second              :', YMDHMS[5, I])
                                            print('Julian day number   :', JulDay[I])
                                            print('TimeDay             :', TimeDay[I])
                                            print('Sun Earth distance  :', DSunEarth[I])
                                            print('Time(UT)            :', TimeUT[I])
                                            print('Time(AES)           :', TimeAES[I])
                                            print('AsetMin             :', AestMin[I])
                                            Halt()

                                    RejectCode[0] = 0

                                while AestMin[NumStack - 1] - AestMin[0] > Window:
                                    I = 0
                                    SolarZenithApp[I] = -1.0
                                    RejectCode[I] = 2
                                    log_file.write(FormatString[10].format(
                                        *YMDHMS[:6, I],
                                        SolarZenithApp[I],
                                        IGroup,
                                        RejectCode[I],
                                        *[SunData[I]['Signal']['Chan'][J] for J in range(NumLangleyChannels(Model))],
                                        SunData[I]['Signal']['Chan'][IWV]
                                    ))

                                    DateTime[:NumStack - 1] = DateTime[1:NumStack]
                                    SunData[:NumStack - 1] = SunData[1:NumStack]

                                    I = NumStack - 1
                                    DateTime[I], SunData[I], ValidData, Eof = Read_Single_Sun_Record(
                                        Dbug[1], InUnit, Model
                                    )
                                    YMDHMS[:, I] = [
                                        DateTime[I].year,
                                        DateTime[I].month,
                                        DateTime[I].day,
                                        DateTime[I].hour,
                                        DateTime[I].minute,
                                        DateTime[I].second,
                                        0
                                    ]

                                    if ApplyClockFix:
                                        Add_Time(YMDHMS[0, I], -TimeCorr)

                                    JulDay[I], TimeDay[I], DSunEarth[I], TimeUT[I], TimeAES[I], AestMin[I] = Julian(YMDHMS[:,I])
                                    SunPos = Solar_Position_Almanac(JulDay[I], TimeDay[I])
                                    DSunEarth[I] = SunPos['DsunEarth']
                                    TimeUT[I] = TimeDay[I] - 0.5
                                    TimeAES[I] = TimeUT[I] * 24 + 10.0
                                    AestMin[I] = TimeAES[I] * 60
                                    if Dbug[2]:
                                        print('Year                :', YMDHMS[0, I])
                                        print('Month               :', YMDHMS[1, I])
                                        print('Day                 :', YMDHMS[2, I])
                                        print('Hour                :', YMDHMS[3, I])
                                        print('Minute              :', YMDHMS[4, I])
                                        print('Second              :', YMDHMS[5, I])
                                        print('Julian day number   :', JulDay[I])
                                        print('TimeDay             :', TimeDay[I])
                                        print('Sun Earth distance  :', DSunEarth[I])
                                        print('Time(UT)            :', TimeUT[I])
                                        print('Time(AES)           :', TimeAES[I])
                                        print('AsetMin             :', AestMin[I])
                                        Halt()

                                    LSTCut = Date_Time(Year, Month, Day, 3, 0, 0)
                                    TimeZoneHMS = Hour2HMS(TimeZone)
                                    UTCut = AddDateTime(LSTCut, Date_Time(
                                        0, 0, 0, -TimeZoneHMS['Hour'], -TimeZoneHMS['Minute'], -TimeZoneHMS['Second']
                                    ))
                                    if DateTime[I] < UTCut:
                                        ValidData = 14

                                    for Ic in range(NumCut):
                                        if DateTime[I] == CutList[Ic]:
                                            ValidData = 13

                                    NewStation = True
                                    SolarZenith, SolarAzimuth = Satellite_Position(
                                        JulDay[I], TimeDay[I],
                                        NewStation,
                                        StationLatRad, StationLonRad,
                                        SunElements
                                    )
                                    SolarZenithDeg[I] = SolarZenith * DegreesPerRadian
                                    SolarZenithApp[I] = Apparent_Zenith(SolarZenithDeg[I])
                                    SolarAzimuthDeg[I] = SolarAzimuth * DegreesPerRadian

                                    if Dbug[2]:
                                        print('Julian day    :', JulDay[I])
                                        print('TimeDay       :', TimeDay[I])
                                        print('Latitude      :', StationLatDeg)
                                        print('Longitude     :', StationLonDeg)
                                        print('Sol Zen (true):', SolarZenithDeg[I])
                                        print('Sol Zen (app) :', SolarZenithApp[I])
                                        print('Solar Azimuth :', SolarAzimuthDeg[I])
                                        Halt()
        
                                    RejectCode[I] = ValidData
                                    while ValidData != 0 and Eof == 0:
                                        log_file.write(FormatString[10].format(
                                            *YMDHMS[:6, I],
                                            SolarZenithApp[I],
                                            IGroup,
                                            RejectCode[I],
                                            *[SunData[I]['Signal']['Chan'][J] for J in range(NumLangleyChannels(Model))],
                                            SunData[I]['Signal']['Chan'][IWV]
                                        ))
                                        DateTime[I], SunData[I], ValidData, Eof = Read_Single_Sun_Record(
                                            Dbug[1], InUnit, Model
                                        )
                                        YMDHMS[:, I] = [
                                            DateTime[I].year,
                                            DateTime[I].month,
                                            DateTime[I].day,
                                            DateTime[I].hour,
                                            DateTime[I].minute,
                                            DateTime[I].second,
                                            0
                                        ]

                                        if ApplyClockFix:
                                            Add_Time(YMDHMS[0, I], -TimeCorr)

                                        JulDay[I], TimeDay[I] = Julian(YMDHMS[:, I])
                                        SunPos = Solar_Position_Almanac(JulDay[I], TimeDay[I])
                                        DSunEarth[I] = SunPos['DsunEarth']
                                        TimeUT[I] = TimeDay[I] - 0.5
                                        TimeAES[I] = TimeUT[I] * 24 + 10.0
                                        AestMin[I] = TimeAES[I] * 60

                                        if Dbug[2]:
                                            print('Year                :', YMDHMS[0, I])
                                            print('Month               :', YMDHMS[1, I])
                                            print('Day                 :', YMDHMS[2, I])
                                            print('Hour                :', YMDHMS[3, I])
                                            print('Minute              :', YMDHMS[4, I])
                                            print('Second              :', YMDHMS[5, I])
                                            print('Julian day number   :', JulDay[I])
                                            print('TimeDay             :', TimeDay[I])
                                            print('Sun Earth distance  :', DSunEarth[I])
                                            print('Time(UT)            :', TimeUT[I])
                                            print('Time(AES)           :', TimeAES[I])
                                            print('AsetMin             :', AestMin[I])
                                            Halt()

                                        LSTCut = Date_Time(Year, Month, Day, 3, 0, 0)
                                        TimeZoneHMS = Hour2HMS(TimeZone)
                                        UTCut = AddDateTime(LSTCut, Date_Time(
                                            0, 0, 0, -TimeZoneHMS['Hour'], -TimeZoneHMS['Minute'], -TimeZoneHMS['Second']
                                        ))
                                        if DateTime[I] < UTCut:
                                            ValidData = 14

                                        for Ic in range(NumCut):
                                            if DateTime[I] == CutList[Ic]:
                                                ValidData = 13

                                        NewStation = True
                                        SolarZenith, SolarAzimuth = Satellite_Position(
                                            JulDay[I], TimeDay[I],
                                            NewStation,
                                            StationLatRad, StationLonRad,
                                            SunElements
                                        )
                                        SolarZenithDeg[I] = SolarZenith * DegreesPerRadian
                                        SolarZenithApp[I] = Apparent_Zenith(SolarZenithDeg[I])
                                        SolarAzimuthDeg[I] = SolarAzimuth * DegreesPerRadian

                                        if Dbug[2]:
                                            print('Julian day    :', JulDay[I])
                                            print('TimeDay       :', TimeDay[I])
                                            print('Latitude      :', StationLatDeg)
                                            print('Longitude     :', StationLonDeg)
                                            print('Sol Zen (true):', SolarZenithDeg[I])
                                            print('Sol Zen (app) :', SolarZenithApp[I])
                                            print('Solar Azimuth :', SolarAzimuthDeg[I])
                                            Halt()

                                        RejectCode[I] = ValidData

                                    if Eof != 0:
                                        break

                                    for I in range(NumStack):
                                        YMDHMS[:, I] = [
                                            DateTime[I].year,
                                            DateTime[I].month,
                                            DateTime[I].day,
                                            DateTime[I].hour,
                                            DateTime[I].minute,
                                            DateTime[I].second,
                                            0
                                        ]

                                        if ApplyClockFix:
                                            Add_Time(YMDHMS[0, I], -TimeCorr)

                                        JulDay[I], TimeDay[I] = Julian(YMDHMS[:, I])
                                        SunPos = Solar_Position_Almanac(JulDay[I], TimeDay[I])
                                        DSunEarth[I] = SunPos['DsunEarth']
                                        TimeUT[I] = TimeDay[I] - 0.5
                                        TimeAES[I] = TimeUT[I] * 24 + 10.0
                                        AestMin[I] = TimeAES[I] * 60

                                        if Dbug[2]:
                                            print('Year                :', YMDHMS[0, I])
                                            print('Month               :', YMDHMS[1, I])
                                            print('Day                 :', YMDHMS[2, I])
                                            print('Hour                :', YMDHMS[3, I])
                                            print('Minute              :', YMDHMS[4, I])
                                            print('Second              :', YMDHMS[5, I])
                                            print('Julian day number   :', JulDay[I])
                                            print('TimeDay             :', TimeDay[I])
                                            print('Sun Earth distance  :', DSunEarth[I])
                                            print('Time(UT)            :', TimeUT[I])
                                            print('Time(AES)           :', TimeAES[I])
                                            print('AsetMin             :', AestMin[I])
                                            Halt()

                                    RejectCode[0] = 0

                                IGroup += 1
                                for N in range(NumChannels):
                                    Volt[:NumStack, N] = np.real(SunData[:NumStack]['Signal']['Chan'][N] - BlackSun[N])
                                    VoltMean[N], VoltSdev[N] = Stat(Volt[:, N], NumStack)
                                    if VoltMean[N] > 0.0:
                                        VoltCv[N] = VoltSdev[N] / VoltMean[N]
                                    else:
                                        VoltCv[N] = 0.9999
                                    if VoltCv[N] > 0.9999:
                                        VoltCv[N] = 0.9999
                                    for I in range(NumStack):
                                        if Volt[I, N] > 0.0:
                                            VoltLog[I, N] = np.log(Volt[I, N])
                                        else:
                                            VoltLog[I, N] = -9.99

                                if Dbug[1]:
                                    print('IGroup  =', IGroup)
                                    print('NumStack=', NumStack)
                                    print('Raw Volts(870)=', [SunData[I]['Signal']['Chan'][I870] for I in range(3)])
                                    print('1AU Volts(870)=', Volt[:NumStack, I870])
                                    print('1AU Mean (870)=', VoltMean[I870])
                                    print('1AU sdev (870)=', VoltSdev[I870])
                                    print('Volt CV  (870)=', VoltCv[I870])
                                    print('Volt Log (1,870)=', VoltLog[0, I870])
                                    Halt()

                                CvOK = VoltCv[I870] <= CvMax
                                SunZenOK = True
                                for I in range(NumStack):
                                    NewStation = True
                                    SolarZenith, SolarAzimuth = Satellite_Position(
                                        JulDay[I], TimeDay[I],
                                        NewStation,
                                        StationLatRad, StationLonRad,
                                        SunElements
                                    )
                                    SolarZenithDeg[I] = SolarZenith * DegreesPerRadian
                                    SolarZenithApp[I] = Apparent_Zenith(SolarZenithDeg[I])
                                    SolarAzimuthDeg[I] = SolarAzimuth * DegreesPerRadian

                                    if Dbug[2]:
                                        print('Julian day    :', JulDay[I])
                                        print('TimeDay       :', TimeDay[I])
                                        print('Latitude      :', StationLatDeg)
                                        print('Longitude     :', StationLonDeg)
                                        print('Sol Zen (true):', SolarZenithDeg[I])
                                        print('Sol Zen (app) :', SolarZenithApp[I])
                                        print('Solar Azimuth :', SolarAzimuthDeg[I])
                                        Halt()

                                    if SolarZenithDeg[I] > SolarZenithMax:
                                        SunZenOK = False

                                    Pressure[I] = GetPressure(
                                        DateTime[I], PData,
                                        PressureOption,
                                        NumPressurePoints,
                                        NumValidPresPoints,
                                        DailyMeanPressure,
                                        DefaultSurfacePressure
                                    )

                                if SolarAzimuthDeg[1] < 180.0:
                                    Iampm = 0
                                else:
                                    Iampm = 1

                                if Dbug[1]:
                                    print('Am/Pm (1/2) :', Iampm)
                                    Halt()

                                for N in range(NumLangleyChannels(Model)):
                                    LnV0[N] = LnV01AU[N] - 2.0 * np.log(DSunEarth[0])
                                    for I in range(NumStack):
                                        TranBoM[I, N] = np.exp(VoltLog[I, N] - LnV0[N])
                                        AirMass[I, N] = (
                                            GetAirMass(1, SolarZenithApp[I]) * RayleighOD[N] +
                                            GetAirMass(2, SolarZenithApp[I]) * 0.03 +
                                            GetAirMass(3, SolarZenithApp[I]) * OzoneOD[N]
                                        ) / (RayleighOD[N] + 0.03 + OzoneOD[N])

                                        if Dbug[3]:
                                            print('Stack index:', I)
                                            print('Channel    :', N)
                                            print('Sun zen app:', SolarZenithApp[I])
                                            print('Rayleigh am:', GetAirMass(1, SolarZenithApp[I]))
                                            print('Aerosol  am:', GetAirMass(2, SolarZenithApp[I]))
                                            print('Ozone    am:', GetAirMass(3, SolarZenithApp[I]))
                                            print('Weighted am:', AirMass[I, N])
                                
                                N = IWV
                                LnV0[N] = LnV01AU[N] - 2.0 * np.log(DSunEarth[0])
                                for I in range(NumStack):
                                    TranBoM[I, N] = np.exp(VoltLog[I, N] - LnV0[N])
                                    AirMass[I, N] = (
                                        GetAirMass(1, SolarZenithApp[I]) * RayleighOD[N] +
                                        GetAirMass(2, SolarZenithApp[I]) * 0.03 +
                                        GetAirMass(3, SolarZenithApp[I]) * OzoneOD[N]
                                    ) / (RayleighOD[N] + 0.03 + OzoneOD[N])

                                    if Dbug[3]:
                                        print('Stack index:', I)
                                        print('Channel    :', N)
                                        print('Sun zen app:', SolarZenithApp[I])
                                        print('Rayleigh am:', GetAirMass(1, SolarZenithApp[I]))
                                        print('Aerosol  am:', GetAirMass(2, SolarZenithApp[I]))
                                        print('Ozone    am:', GetAirMass(3, SolarZenithApp[I]))
                                        print('Weighted am:', AirMass[I, N])
                                        Halt()

                                if SunZenOK:
                                    if CvOK:
                                        IGroup1 += 1
                                        for I in range(NumStack):
                                            Is += 1
                                            for N in range(NumLangleyChannels(Model)):
                                                RayleighOD[N] = Rayleigh(Dbug[2], Wavelength[N, Model], Pressure[I])
                                                if Dbug[2]:
                                                    print('Wavelength:', Wavelength[N, Model])
                                                    print('Pressure  :', Pressure[I])
                                                    print('RayleighOD:', RayleighOD[N])

                                                AirMass[I, N] = (
                                                    GetAirMass(1, SolarZenithApp[I]) * RayleighOD[N] +
                                                    GetAirMass(2, SolarZenithApp[I]) * 0.03 +
                                                    GetAirMass(3, SolarZenithApp[I]) * OzoneOD[N]
                                                ) / (RayleighOD[N] + 0.03 + OzoneOD[N])

                                                if Dbug[3]:
                                                    print('Stack index:', I)
                                                    print('Channel    :', N)
                                                    print('Sun zen app:', SolarZenithApp[I])
                                                    print('Rayleigh am:', GetAirMass(1, SolarZenithApp[I]))
                                                    print('Aerosol  am:', GetAirMass(2, SolarZenithApp[I]))
                                                    print('Ozone    am:', GetAirMass(3, SolarZenithApp[I]))
                                                    print('Weighted am:', AirMass[I, N])
                                                    Halt()

                                                LnV0[N] = LnV01AU[N] - 2.0 * np.log(DSunEarth[0])
                                                Tran[Is, N] = np.exp(VoltLog[I, N] - LnV0[N])
                                                Aod1[Is, N] = (
                                                    (LnV0[N] - VoltLog[I, N]) / AirMass[I, N] -
                                                    RayleighOD[N] - OzoneOD[N] - ExtraOD[N]
                                                )

                                                if Dbug[3]:
                                                    print('Stack index ', I)
                                                    print('Channel     ', N)
                                                    print('LnV0        ', LnV0[N])
                                                    print('LnV         ', VoltLog[I, N])
                                                    print('Airmass     ', AirMass[I, N])
                                                    print('RayleighOD  ', RayleighOD[N])
                                                    print('OzoneOD     ', OzoneOD[N])
                                                    print('ExtraOD     ', ExtraOD[N])
                                                    print('AerosolOD   ', Aod1[Is, N])
                                                    Halt()

                                                StdUAod[Is, N] = StdUAod1AM[N] / AirMass[I, N]
                                                SignalCV[Is, N] = VoltCv[N]

                                            Ratio = Aod1[Is, I670] / Aod1[Is, I440]
                                            if Aod1[Is, I670] > 0.0 and Aod1[Is, I440] > 0.0:
                                                WlRatioLog = np.log(Wavelength[I670, Model] / Wavelength[I440, Model])
                                                Angstrom46[Is] = -np.log(Ratio) / WlRatioLog
                                                StdUAngstrom46[Is] = (
                                                    StdUAod[Is, I440] / Aod1[Is, I440] +
                                                    StdUAod[Is, I670] / Aod1[Is, I670]
                                                ) / WlRatioLog
                                            else:
                                                Angstrom46[Is] = -9.999
                                                StdUAngstrom46[Is] = -9.999

                                            Angstrom46[Is] = max(Angstrom46[Is], -9.99)
                                            Angstrom46[Is] = min(Angstrom46[Is], +9.99)

                                            Ratio = Aod1[Is, I870] / Aod1[Is, I440]
                                            if Aod1[Is, I870] > 0.0 and Aod1[Is, I440] > 0.0:
                                                WlRatioLog = np.log(Wavelength[I870, Model] / Wavelength[I440, Model])
                                                Angstrom48[Is] = -np.log(Ratio) / WlRatioLog
                                                StdUAngstrom48[Is] = (
                                                    StdUAod[Is, I440] / Aod1[Is, I440] +
                                                    StdUAod[Is, I870] / Aod1[Is, I870]
                                                ) / WlRatioLog
                                            else:
                                                Angstrom48[Is] = -9.999
                                                StdUAngstrom48[Is] = -9.999

                                            Angstrom48[Is] = max(Angstrom48[Is], -9.99)
                                            Angstrom48[Is] = min(Angstrom48[Is], +9.99)

                                            Ratio = Aod1[Is, I870] / Aod1[Is, I670]
                                            if Aod1[Is, I670] > 0.0 and Aod1[Is, I870] > 0.0:
                                                WlRatioLog = np.log(Wavelength[I870, Model] / Wavelength[I670, Model])
                                                Angstrom68[Is] = -np.log(Ratio) / WlRatioLog
                                                StdUAngstrom68[Is] = (
                                                    StdUAod[Is, I670] / Aod1[Is, I670] +
                                                    StdUAod[Is, I870] / Aod1[Is, I870]
                                                ) / WlRatioLog
                                            else:
                                                Angstrom68[Is] = -9.999
                                                StdUAngstrom68[Is] = -9.999

                                            Angstrom68[Is] = max(Angstrom68[Is], -9.99)
                                            Angstrom68[Is] = min(Angstrom68[Is], +9.99)

                                            if CirrusOption == 0:
                                                CirrusFlag[Is] = -1
                                            elif CirrusOption == 1:
                                                if Aod1[Is, I440] > 0:
                                                    CirrusAngstromThreshold = 0.862 + 0.556 * np.log10(Aod1[Is, I440])
                                                else:
                                                    CirrusAngstromThreshold = 0.5
                                                if Angstrom48[Is] < CirrusAngstromThreshold:
                                                    CirrusFlag[Is] = 1
                                                else:
                                                    CirrusFlag[Is] = 0
                                            elif CirrusOption == 2:
                                                if Aod1[Is, I1020] > Aod1[Is, I870]:
                                                    CirrusFlag[Is] = 1
                                                else:
                                                    CirrusFlag[Is] = 0
                                            elif CirrusOption >= 3 and CirrusOption <= 5:
                                                CirrusFlag[Is] = -1

                                            LnV0[IWV] = LnV01AU[IWV] - 2.0 * np.log(DSunEarth[0])
                                            RayleighOD[IWV] = Rayleigh(Dbug[2], Wavelength[IWV, Model], Pressure[I])

                                            if Dbug[2]:
                                                print('Wavelength:', Wavelength[IWV, Model])
                                                print('Pressure  :', Pressure[I])
                                                print('RayleighOD:', RayleighOD[IWV])

                                            AirMass[I, IWV] = (
                                                GetAirMass(1, SolarZenithApp[I]) * RayleighOD[IWV] +
                                                GetAirMass(2, SolarZenithApp[I]) * 0.03 +
                                                GetAirMass(3, SolarZenithApp[I]) * OzoneOD[IWV]
                                            ) / (RayleighOD[IWV] + 0.03 + OzoneOD[IWV])

                                            if Dbug[2] and I == 0:
                                                print('Rayleigh am:', GetAirMass(1, SolarZenithApp[I]))
                                                print('Aerosol  am:', GetAirMass(2, SolarZenithApp[I]))
                                                print('Ozone    am:', GetAirMass(3, SolarZenithApp[I]))
                                                print('Weighted am:', AirMass[I, IWV])
                                                Halt()

                                            Aod1[Is, IWV] = (
                                                Aod1[Is, I870] *
                                                (Wavelength[IWV, Model] / Wavelength[I870, Model]) **
                                                (-Angstrom68[Is])
                                            )
                                            Tran[Is, IWV] = np.exp(VoltLog[I, IWV] - LnV0[IWV])
                                            TRalAsl = np.exp(-AirMass[I, IWV] * (RayleighOD[IWV] + Aod1[Is, IWV]))
                                            Tran[Is, IWV] = Tran[Is, IWV] * TRalAsl

                                            Y = (
                                                VoltLog[I, IWV] +
                                                AirMass[I, IWV] * (RayleighOD[IWV] + Aod1[Is, IWV])
                                            )
                                            WVarg = (LnV0[IWV] - Y) / Acoef
                                            if WVarg > 0:
                                                Wvap[Is] = WVarg ** (1.0 / Bcoef) / AirMass[I, IWV]
                                            else:
                                                Wvap[Is] = -9.99

                                            if Dbug[3]:
                                                print('Stack index ', I)
                                                print('WV Channel     ')
                                                print('LnV0        ', LnV0[IWV])
                                                print('LnV         ', VoltLog[I, IWV])
                                                print('Airmass     ', AirMass[I, IWV])
                                                print('RayleighOD  ', RayleighOD[IWV])
                                                print('AOD670      ', Aod1[Is, I670])
                                                print('AOD870      ', Aod1[Is, I870])
                                                print('AOD936      ', Aod1[Is, IWV])
                                                print('LnV+m(tr+ta)', Y)
                                                print('Acoef       ', Acoef)
                                                print('Bcoef       ', Bcoef)
                                                print('WVarg       ', WVarg)
                                                print('Wvap        ', Wvap[Is])
                                                Halt()

                                            WvapOD[Is, I1020] = H2O_k * Wvap[Is] ** H2O_e
                                            N = NumLangleyChannels(Model)
                                            OutRecord[Is]['DateTime']['Year'] = YMDHMS[0, I]
                                            OutRecord[Is]['DateTime']['Month'] = YMDHMS[1, I]
                                            OutRecord[Is]['DateTime']['Day'] = YMDHMS[2, I]
                                            OutRecord[Is]['DateTime']['Hour'] = YMDHMS[3, I]
                                            OutRecord[Is]['DateTime']['Minute'] = YMDHMS[4, I]
                                            OutRecord[Is]['DateTime']['Second'] = YMDHMS[5, I]
                                            OutRecord[Is]['SunZen'] = SolarZenithDeg[I]
                                            OutRecord[Is]['AirMass'] = AirMass[I, I870]
                                            OutRecord[Is]['Pressure'] = Pressure[I]
                                            OutRecord[Is]['Temperature'] = SunData[I]['Temperature']
                                            OutRecord[Is]['Aod'][:NumChannels] = -1.0
                                            OutRecord[Is]['Aod'][:N] = Aod1[Is, :N]
                                            OutRecord[Is]['StdUAod'][:NumChannels] = -1.0
                                            OutRecord[Is]['StdUAod'][:N] = StdUAod[Is, :N]
                                            OutRecord[Is]['Wvap'] = Wvap[Is]
                                            OutRecord[Is]['StdUWvap'] = -1.0
                                            OutRecord[Is]['Aod1020c'] = (
                                                Aod1[Is, I1020] - WvapOD[Is, I1020]
                                            )
                                            OutRecord[Is]['Angstrom46'] = Angstrom46[Is]
                                            OutRecord[Is]['StdUAngstrom46'] = StdUAngstrom46[Is]
                                            OutRecord[Is]['Angstrom48'] = Angstrom48[Is]
                                            OutRecord[Is]['StdUAngstrom48'] = StdUAngstrom48[Is]
                                            OutRecord[Is]['Angstrom68'] = Angstrom68[Is]
                                            OutRecord[Is]['StdUAngstrom68'] = StdUAngstrom68[Is]
                                            OutRecord[Is]['Tran'][:NumChannels] = -1.0
                                            OutRecord[Is]['Tran'][:N] = Tran[Is, :N]
                                            OutRecord[Is]['Tran'][IWV] = Tran[Is, IWV]
                                            OutRecord[Is]['CirrusFlag'] = CirrusFlag[Is]
                                    for N in range(NumLangleyChannels(Model)):
                                        AodGroupMean[IGroup1, N], AodGroupSdev[IGroup1, N] = Stat(
                                            Aod1[Is - NumStack + 1:Is + 1, N], NumStack
                                        )
                                else:
                                    RejectCode[:NumStack] = 3
                            else:
                                RejectCode[:NumStack] = 4

                            for I in range(NumStack):
                                ObsTime_0UT = (DateTime[I] - Date_Time(Year, Month, Day, 0, 0, 0)).total_seconds()
                                ObsTime_0LST = ObsTime_0UT / 60.0 + TimeZone * 60
                                ObsTime_UTDay = ObsTime_0UT / 3600 / 24

                                FLT = 0
                                StdDev = VoltSdev[I870] / np.exp(LnV0[I870]) * 1000
                                Aureole[:NumChannels] = 0

                                bom_file.write(FormatString[9].format(
                                    ObsDate,
                                    Day_of_Year(Day, Month, Year),
                                    ObsTime_0LST,
                                    ObsTime_UTDay,
                                    Pressure[I],
                                    SolarZenithDeg[I],
                                    AirMass[I, I870],
                                    FLT,
                                    StdDev,
                                    *[TranBoM[I, J] for J in range(NumLangleyChannels(Model))],
                                    *Aureole[:NumLangleyChannels(Model)],
                                    TranBoM[I, IWV],
                                    Aureole[IWV],
                                    CR
                                ))

                            for I in range(NumStack):
                                log_file.write(FormatString[10].format(
                                    *YMDHMS[:6, I],
                                    SolarZenithApp[I],
                                    IGroup,
                                    RejectCode[I],
                                    *[SunData[I]['Signal']['Chan'][J] for J in range(NumLangleyChannels(Model))],
                                    SunData[I]['Signal']['Chan'][IWV]
                                ))

                        NObsOK = Is
                        file.close()
                        log_file.close()
                        bom_file.close()

                        if NObsOK >= MinObs:
                            print('Pass 2 filtering...')
                            for N in range(NumLangleyChannels(Model)):
                                AodDayMean[N], AodDaySdev[N] = Stat(Aod1[:, N], Is)

                            Aod2[:NObsOK, :NumLangleyChannels(Model)] = Aod1[:NObsOK, :NumLangleyChannels(Model)]

                            if AodDaySdev[I870] > AodDaySdevCrit:
                                NGroupOK = IGroup1
                                Reject = False
                                Ns = 0
                                Ng = 0
                                for Ig in range(1, NGroupOK + 1):
                                    Is = (Ig - 1) * NumStack
                                    if np.abs(AodGroupMean[Ig, I870] - AodDayMean[I870]) < 3 * AodDaySdev[I870]:
                                        Ng += 1
                                        AodGroupMean[Ng, I870] = AodGroupMean[Ig, I870]
                                        for I in range(1, NumStack + 1):
                                            Is += 1
                                            Ns += 1
                                            Aod2[Ns, :NumLangleyChannels(Model)] = Aod2[Is, :NumLangleyChannels(Model)]
                                            OutRecord[Ns] = OutRecord[Is]
                                    else:
                                        Reject = True

                                print('After filt:Num groups ok=', Ng)
                                NGroupOK = Ng
                                NObsOK = Ns

                                while Reject:
                                    for N in range(NumLangleyChannels(Model)):
                                        AodDayMean[N], AodDaySdev[N] = Stat(Aod2[:, N], NObsOK)

                                    print('Revised aod870mean=', AodDayMean[I870])
                                    print('Revised aod870sdev=', AodDaySdev[I870])
                                    print('Revised 3sdev     =', 3 * AodDaySdev[I870])

                                    Reject = False
                                    Ns = 0
                                    Ng = 0
                                    for Ig in range(1, NGroupOK + 1):
                                        Is = (Ig - 1) * NumStack
                                        if np.abs(AodGroupMean[Ig, I870] - AodDayMean[I870]) < 3 * AodDaySdev[I870]:
                                            Ng += 1
                                            AodGroupMean[Ng, I870] = AodGroupMean[Ig, I870]
                                            for I in range(1, NumStack + 1):
                                                Is += 1
                                                Ns += 1
                                                Aod2[Ns, :NumLangleyChannels(Model)] = Aod2[Is, :NumLangleyChannels(Model)]
                                                OutRecord[Ns] = OutRecord[Is]
                                        else:
                                            Reject = True

                                    print('After filt:Num groups ok=', Ng)
                                    NGroupOK = Ng
                                    NObsOK = Ns

                            J = 0
                            K = 0
                            for I in range(1, NObsOK + 1):
                                N = NumLangleyChannels(Model)
                                AodTemp[:N] = Aod2[I, :N]

                                if (
                                    np.abs(OutRecord[I]['StdUAod'][I440] / OutRecord[I]['Aod'][I440]) > Rel_Uaod_440_Threshold or
                                    np.abs(OutRecord[I]['StdUAngstrom48']) > Uangst48_Threshold
                                ):
                                    continue

                                J += 1
                                OutRecord[J] = OutRecord[I]
                                AodDay[J, :N] = OutRecord[J]['Aod'][:N]
                                Angstrom48Day[J] = OutRecord[J]['Angstrom48']

                                if OutRecord[J]['Wvap'] < 10.0:
                                    K += 1
                                    WvapDay[K] = OutRecord[J]['Wvap']
                                    NWvMonth += 1
                                    WvapMonth[NWvMonth] = WvapDay[K]
                                else:
                                    OutRecord[J]['Wvap'] = 9.999

                                NObsMonth += 1
                                AodMonth[NObsMonth, :N] = AodDay[J, :N]
                                Angstrom48Month[NObsMonth] = Angstrom48Day[J]

                            NObsOK = J
                            NWvOK = K

                            if CirrusOption == 3:
                                NumCodPoints, Cod = Read_Cod_File(Dbug[0], CodFile, NumLangleyChannels(Model))
                                for Is in range(1, NObsOK + 1):
                                    CirrusFlag[Is] = -1
                                    for I in range(1, NumCodPoints + 1):
                                        if Cod[I]['DateTime'] == OutRecord[Is]['DateTime']:
                                            if Cod[I]['Aod_500_Coarse'] > Aod_500_Coarse_Threshold:
                                                CirrusFlag[Is] = 1
                                            else:
                                                CirrusFlag[Is] = 0
                                            OutRecord[Is]['CirrusFlag'] = CirrusFlag[Is]
                                    if CirrusFlag[Is] == -1:
                                        print('Trouble, CirrusFlag not set')
                                        Halt()

                            if CirrusOption >= 4:
                                for Is in range(1, NObsOK + 1):
                                    CirrusFlag[Is] = -1
                                    MPL_Current = False
                                    MPL, Eof = Read_MPL_Record(Dbug[0], MPLUnit)
                                    MPL_Current = MPL['CloudBackgroundRatio'] > 0.0
                                    while Eof == 0 and (MPL['DateTime'] - OutRecord[Is]['DateTime']).total_seconds() < 0:
                                        MPL, Eof = Read_MPL_Record(Dbug[0], MPLUnit)
                                        MPL_Previous = MPL_Current
                                        MPL_Current = MPL['CloudBackgroundRatio'] > 0.0
                                    if Eof != 0:
                                        print('Trouble, No time match found in MPL file')
                                        Halt()
                                    if MPL_Previous or MPL_Current:
                                        CirrusFlag[Is] = 1
                                    else:
                                        CirrusFlag[Is] = 0
                                    OutRecord[Is]['CirrusFlag'] = CirrusFlag[Is]

                            if CirrusOption == 5:
                                NumCodPoints, Cod = Read_Cod_File(Dbug[0], CodFile, NumLangleyChannels(Model))
                                for Is in range(1, NObsOK + 1):
                                    if CirrusFlag[Is] == 1:
                                        print('Numcodpoints=', NumCodPoints)
                                        print('cirrus flag set, now looking for cod match')
                                        Nmatch = 0
                                        for I in range(1, NumCodPoints + 1):
                                            if Cod[I]['DateTime'] == OutRecord[Is]['DateTime']:
                                                print('cod match found')
                                                Nmatch += 1
                                                CirrusAodUnit.write(FormatString[11].format(
                                                    OutRecord[Is]['DateTime'],
                                                    OutRecord[Is]['SunZen'],
                                                    OutRecord[Is]['AirMass'],
                                                    OutRecord[Is]['Pressure'],
                                                    OutRecord[Is]['Temperature'],
                                                    OutRecord[Is]['Aod'][I500[Model]],
                                                    Cod[I]['Aod_500_Fine'],
                                                    Cod[I]['Aod_500_Coarse']
                                                ))
                                        if Nmatch == 0:
                                            print('NO cod match found')
                                            Halt()

                            OutFile = Make_File_Name(ResultPath, OutRoot, 'aod')
                            print('OutFile=', OutFile)
                            OutUnit = Next_Unit_Number()
                            with open(OutFile, 'w') as out_file:
                                out_file.write(FormatString[12].format(
                                    'Aerosol optical depth',
                                    SiteCode,
                                    StationLatDeg,
                                    StationLonDeg,
                                    Instrument,
                                    Model,
                                    CalFile,
                                    OzoneFile,
                                    NumLangleyChannels(Model) + 1,
                                    OzoneColumn,
                                    NumLangleyChannels(Model) + 1
                                ))
                                out_file.write(FormatString[13].format(NumLangleyChannels(Model)))
                                out_file.write(FormatString[14].format(NumLangleyChannels(Model)))
                                out_file.write(FormatString[15].format(
                                    *[IWavelength[I] for I in range(NumLangleyChannels(Model))],
                                    IWavelength[IWV],
                                    *[OzoneOD[I] for I in range(NumLangleyChannels(Model))],
                                    OzoneOD[IWV]
                                ))
                                out_file.write(FormatString[16].format(NumLangleyChannels(Model)))
                                out_file.write(FormatString[17].format(
                                    *[IWavelength[I] for I in range(NumLangleyChannels(Model))],
                                    IWavelength[IWV],
                                    IWavelength[I1020]
                                ))

                                N = NumLangleyChannels(Model)
                                for Is in range(1, NObsOK + 1):
                                    out_file.write(FormatString[14].format(
                                        OutRecord[Is]['DateTime'],
                                        OutRecord[Is]['SunZen'],
                                        OutRecord[Is]['AirMass'],
                                        OutRecord[Is]['Pressure'],
                                        OutRecord[Is]['Temperature'],
                                        *[OutRecord[Is]['Aod'][J] for J in range(N)],
                                        *[OutRecord[Is]['StdUAod'][J] for J in range(N)],
                                        OutRecord[Is]['Wvap'],
                                        OutRecord[Is]['Aod1020c'],
                                        OutRecord[Is]['Angstrom46'],
                                        OutRecord[Is]['StdUAngstrom46'],
                                        OutRecord[Is]['Angstrom48'],
                                        OutRecord[Is]['StdUAngstrom48'],
                                        OutRecord[Is]['Angstrom68'],
                                        OutRecord[Is]['StdUAngstrom68'],
                                        OutRecord[Is]['CirrusFlag']
                                    ))
                                out_file.close()

                            if MakeTran:
                                if EnvFound:
                                    for Is in range(1, NObsOK + 1):
                                        ObsTime = OutRecord[Is]['DateTime']
                                        EnvOK = False
                                        DtEnv[0] = ObsTime - Env[0]['DateTime']
                                        if DtEnv[0] >= 0:
                                            for J in range(1, EnvRecs):
                                                DtEnv[J] = ObsTime - Env[J]['DateTime']
                                                DtEnv[J + 1] = ObsTime - Env[J + 1]['DateTime']
                                                if np.abs(DtEnv[J]) > DtMax:
                                                    DtEnv[J] = np.sign(DtEnv[J]) * DtMax
                                                if np.abs(DtEnv[J + 1]) > DtMax:
                                                    DtEnv[J + 1] = np.sign(DtEnv[J + 1]) * DtMax
                                                if DtEnv[J] * DtEnv[J + 1] <= 0:
                                                    if np.abs(DtEnv[J]) < np.abs(DtEnv[J + 1]):
                                                        IEnv = J
                                                    else:
                                                        IEnv = J + 1
                                                    EnvOK = True
                                                    break
                                            if not EnvOK:
                                                print('No match reading Env data, time after end')
                                        else:
                                            print('No match reading Env data, time prior to start')
                                            IEnv = 0
                                            Env[IEnv]['Temperature']['Mean'] = 0
                                        OutRecord[Is]['Temperature'] = Env[IEnv]['Temperature']['Mean']
                                else:
                                    EnvOK = False

                                if not EnvOK or not EnvFound:
                                    IEnv = 0
                                    DtEnv[IEnv] = 0
                                    Env[IEnv]['Pressure']['Mean'] = 0
                                    Env[IEnv]['Temperature']['Mean'] = 0
                                    Env[IEnv]['WindV']['Mean'] = 0
                                    Env[IEnv]['WindV']['Sdev'] = 0
                                    Env[IEnv]['WindD']['Mean'] = 0
                                    Env[IEnv]['WindD']['Sdev'] = 0

                                TrnFile = Make_File_Name(ResultPath, OutRoot, 'trn')
                                print('TrnFile=', TrnFile)
                                TrnUnit = Next_Unit_Number()
                                with open(TrnFile, 'w') as trn_file:
                                    trn_file.write(FormatString[18].format(
                                        'Transmittance',
                                        SiteCode,
                                        StationLatDeg,
                                        StationLonDeg,
                                        Instrument,
                                        Model,
                                        CalFile,
                                        NumLangleyChannels(Model) + 1,
                                        NumLangleyChannels(Model) + 1,
                                        NumLangleyChannels(Model) + 1
                                    ))
                                    trn_file.write(FormatString[19].format(
                                        *[IWavelength[I] for I in range(NumLangleyChannels(Model))],
                                        IWavelength[IWV],
                                        *[V01AU[I] for I in range(NumLangleyChannels(Model))],
                                        V01AU[IWV],
                                        *[V0Erms[I] for I in range(NumLangleyChannels(Model))],
                                        V0Erms[IWV]
                                    ))
                                    trn_file.write(FormatString[20].format(NumLangleyChannels(Model) + 1))
                                    trn_file.write(FormatString[21].format(
                                        *[IWavelength[I] for I in range(NumLangleyChannels(Model))],
                                        IWavelength[IWV]
                                    ))
                                    N = NumLangleyChannels(Model)
                                    for Is in range(1, NObsOK + 1):
                                        trn_file.write(FormatString[14].format(
                                            OutRecord[Is]['DateTime'],
                                            OutRecord[Is]['SunZen'],
                                            OutRecord[Is]['AirMass'],
                                            OutRecord[Is]['Pressure'],
                                            OutRecord[Is]['Temperature'],
                                            *[OutRecord[Is]['Tran'][I] for I in range(N)],
                                            OutRecord[Is]['Tran'][IWV]
                                        ))
                                    trn_file.close()

                            if NObsOK >= MinObs:
                                for N in range(NumLangleyChannels(Model)):
                                    AodDayMean[N], AodDaySdev[N], AodDayMin[N], AodDayMax[N] = Stat(
                                        AodDay[:, N], NObsOK
                                    )
                                if Dbug[3]:
                                    print('NObsOK         =', NObsOK)
                                    print('AodDayMean(870)=', AodDayMean[I870])
                                    print('AodDaySdev(870)=', AodDaySdev[I870])
                                    print('AodDayMin (870)=', AodDayMin[I870])
                                    print('AodDayMax (870)=', AodDayMax[I870])
                                    Halt()

                                Angstrom48DayMean, Angstrom48DaySdev, Angstrom48DayMin, Angstrom48DayMax = Stat(
                                    Angstrom48Day[:NObsOK], NObsOK
                                )
                                WvapDayMean, WvapDaySdev, WvapDayMin, WvapDayMax = Stat(
                                    WvapDay[:NWvOK], NWvOK
                                )
                                if Dbug[3]:
                                    print('WvapDayMean=', WvapDayMean)
                                    print('WvapDaySdev=', WvapDaySdev)
                                    print('WvapDayMin =', WvapDayMin)
                                    print('WvapDayMax =', WvapDayMax)
                                    Halt()

                                TDet[:NObsOK] = [OutRecord[I]['Temperature'] for I in range(1, NObsOK + 1)]
                                TDetDayMean, TDetDaySdev, TDetDayMin, TDetDayMax = Stat(
                                    TDet[:NObsOK], NObsOK
                                )
                                if Dbug[3]:
                                    print('Is  =', IGroup)
                                    print('TDetDayMean=', TDetDayMean)
                                    print('TDetDaySdev=', TDetDaySdev)
                                    print('TDetDayMin =', TDetDayMin)
                                    print('TDetDayMax =', TDetDayMax)
                                    Halt()

                                DayUnit.write(FormatString[22].format(
                                    Year, Month, Day,
                                    Day_of_Year(Day, Month, Year),
                                    NObsOK,
                                    *[AodDayMean[N] for N in range(NumLangleyChannels(Model))],
                                    *[AodDayMax[N] for N in range(NumLangleyChannels(Model))],
                                    *[AodDayMin[N] for N in range(NumLangleyChannels(Model))],
                                    *[AodDaySdev[N] for N in range(NumLangleyChannels(Model))],
                                    WvapDayMean, WvapDayMax,
                                    WvapDayMin, WvapDaySdev,
                                    Angstrom48DayMean,
                                    Angstrom48DayMax,
                                    Angstrom48DayMin,
                                    Angstrom48DaySdev
                                ))

                        if CirrusOption != 0:
                            SsxFile = Make_File_Name(ResultPath, OutRoot, 'ssx')
                            print('SsxFile=', SsxFile)
                            SsxUnit = Next_Unit_Number()
                            with open(SsxFile, 'w') as ssx_file:
                                ssx_file.write(FormatString[23])
                                for Is in range(2, NObsOK + 1):
                                    if np.abs(OutRecord[Is]['Aod'][I440] - OutRecord[Is - 1]['Aod'][I440]) >= 0.001:
                                        DAlphaDTau[Is] = (
                                            (OutRecord[Is]['Angstrom48'] - OutRecord[Is - 1]['Angstrom48']) /
                                            (OutRecord[Is]['Aod'][I440] - OutRecord[Is - 1]['Aod'][I440])
                                        )
                                    else:
                                        DAlphaDTau[Is] = -9.999
                                    RCorr[Is] = -9.999
                                RCorrHist[:] = 0
                                for Is in range(NWin // 2 + 1, NObsOK - NWin // 2 + 1):
                                    RCorr[Is] = Correlation_Coefficient(
                                        OutRecord[Is - NWin // 2:Is + NWin // 2]['Aod'][I440],
                                        OutRecord[Is - NWin // 2:Is + NWin // 2]['Angstrom48']
                                    )
                                    RCorrHist[int(10 * RCorr[Is] + 11)] += 1
                                for Is in range(2, NObsOK + 1):
                                    ssx_file.write(FormatString[24].format(
                                        OutRecord[Is]['DateTime'],
                                        DAlphaDTau[Is], RCorr[Is]
                                    ))
                                ssx_file.write('\n')
                                ssx_file.write('\n')
                                ssx_file.write(FormatString[25])
                                for I in range(1, RCorrHistLength + 1):
                                    ssx_file.write(FormatString[26].format(
                                        I, (I - 11) / 10, (I - 10) / 10, RCorrHist[I - 1]
                                    ))
                                ssx_file.close()

                    if NObsOK >= MinObs:
                        NumDaysOK += 1

        PossibleDayFrac = (DMax - DMin + 1) / Days_in_Month(Month, Year)

        NumDaysOKMin[0] = 8 * PossibleDayFrac
        NumDaysOKMin[1] = 6 * PossibleDayFrac
        NObsMonthMin[0] = 300 * PossibleDayFrac
        NObsMonthMin[1] = 600 * PossibleDayFrac

        MonthStatsFlag = (
            (NumDaysOK >= NumDaysOKMin[0] and NObsMonth >= NObsMonthMin[0]) or
            (NumDaysOK >= NumDaysOKMin[1] and NObsMonth >= NObsMonthMin[1])
        )

        print('Month          =', Month)
        print('PossibleDayFrac=', PossibleDayFrac)
        print('NumDaysOK, NObsMonth:', NumDaysOK, NObsMonth)
        print('Must equal or exceed:', NumDaysOKMin[0], NObsMonthMin[0])
        print('OR                  :', NumDaysOKMin[1], NObsMonthMin[1])
        print('MonthStatsFlag=', MonthStatsFlag)

        if MonthStatsFlag:
            for N in range(NumLangleyChannels(Model)):
                AodMonthMean[N], AodMonthSdev[N], AodMonthMin[N], AodMonthMax[N] = Stat(
                    AodMonth[:NObsMonth, N], NObsMonth
                )
            if Dbug[3]:
                print('Month            =', Month)
                print('NObsMonth        =', NObsMonth)
                print('AodMonthMean(870)=', AodMonthMean[I870])
                print('AodMonthSdev(870)=', AodMonthSdev[I870])
                print('AodMonthMin (870)=', AodMonthMin[I870])
                print('AodMonthMax (870)=', AodMonthMax[I870])
                Halt()

            Angstrom48MonthMean, Angstrom48MonthSdev, Angstrom48MonthMin, Angstrom48MonthMax = Stat(
                Angstrom48Month[:NObsMonth], NObsMonth
            )
            if Dbug[3]:
                print('Month        =', Month)
                print('NObsMonth    =', NObsMonth)
                print('Angstrom48MonthMean=', Angstrom48MonthMean)
                print('Angstrom48MonthSdev=', Angstrom48MonthSdev)
                print('Angstrom48MonthMin =', Angstrom48MonthMin)
                print('Angstrom48MonthMax =', Angstrom48MonthMax)
                Halt()

            WvapMonthMean, WvapMonthSdev, WvapMonthMin, WvapMonthMax = Stat(
                WvapMonth[:NWvMonth], NWvMonth
            )
            if Dbug[3]:
                print('Month        =', Month)
                print('NWvMonth     =', NWvMonth)
                print('WvapMonthMean=', WvapMonthMean)
                print('WvapMonthSdev=', WvapMonthSdev)
                print('WvapMonthMin =', WvapMonthMin)
                print('WvapMonthMax =', WvapMonthMax)
                Halt()

            MonthUnit.write(FormatString[27].format(
                Year, Month,
                ThisMonthsDayStamp,
                NumDaysOK,
                NObsMonth,
                *[AodMonthMean[N] for N in range(NumLangleyChannels(Model))],
                *[AodMonthMax[N] for N in range(NumLangleyChannels(Model))],
                *[AodMonthMin[N] for N in range(NumLangleyChannels(Model))],
                *[AodMonthSdev[N] for N in range(NumLangleyChannels(Model))],
                WvapMonthMean, WvapMonthMax,
                WvapMonthMin, WvapMonthSdev,
                Angstrom48MonthMean, Angstrom48MonthMax,
                Angstrom48MonthMin, Angstrom48MonthSdev
            ))
            TeXUnit.write(FormatString[28].format(
                Year, Month,
                NObsMonth,
                *[AodMonthMean[Isort[N]] for N in range(NumLangleyChannels(Model))],
                *[AodMonthSdev[Isort[N]] for N in range(NumLangleyChannels(Model))],
                Angstrom48MonthMean, Angstrom48MonthSdev,
                WvapMonthMean, WvapMonthSdev
             ))
        else:
            XMonthUnit.write('{:4d}{:02d}\n'.format(Year, Month))
            NObsMonth = 0
            NWvMonth = 0
            NumDaysOK = 0

DayUnit.close()
MonthUnit.close()
XMonthUnit.close()
TeXUnit.close()
OzoneUnit.close()
print('Done.')
