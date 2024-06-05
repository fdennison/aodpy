import math
import datetime

# Constants
MaxLangleys = 2048
MaxPolyOrder = 8
NumChannels = 8
NumModels = 3
MaxPressurePoints = 100
MaxTriplets = 1000
MaxSinglets = 3000
MaxCutPeriods = 100
MinValidPresPoints = 4
RadiansPerDegree = math.pi / 180.0
DegreesPerRadian = 180.0 / math.pi
I670 = 3
I870 = 4
IWV = 8
NumSkyA = 4
NumSkyK = 4

# Function to calculate apparent zenith angle
def Apparent_Zenith(true_zenith):
    return true_zenith

# Function to calculate air mass
def GetAirMass(type, zenith_angle):
    if type == 1:
        return 1.0 / math.cos(zenith_angle * RadiansPerDegree)
    elif type == 2:
        return 1.0 / math.cos(zenith_angle * RadiansPerDegree)
    elif type == 3:
        return 1.0 / math.cos(zenith_angle * RadiansPerDegree)

# Function to calculate pressure
def GetPressure(datetime, pdata, pressure_option, num_pressure_points, num_valid_pres_points, daily_mean_pressure, default_surface_pressure):
    if pressure_option == 1:
        # Logger pressure
        if num_valid_pres_points > 0:
            return daily_mean_pressure
        else:
            return default_surface_pressure
    elif pressure_option == 3:
        # BoM pressure
        if num_valid_pres_points > 0:
            return daily_mean_pressure
        else:
            return default_surface_pressure

# Function to calculate Rayleigh optical depth
def Rayleigh(dbug, wavelength, pressure, rayleigh_od):
    pass

# Function to calculate statistics
def Stat(data, n, mean, sdev):
    pass

# Function to perform linear regression
def Elfit(n, weight, x, y, intercept, slope, residual, erms, del_intercept, del_slope):
    pass

# Main program
def main():
    # Initialize variables
    dbug = [False] * 5
    lang_flag = [False] * 2
    new_station = False
    make_cal = False
    apply_clock_fix = False
    clock_fix_file_ok = False
    include = False
    use_ref_channel4qa = False
#    instrument = 0
#    model = 0
#    iref = 0
#    year = 0
#    month = 0
#    day = 0
#    year_min = 0
#    year_max = 0
#    month_min = 0
#    month_max = 0
#    mmin = 0
#    mmax = 0
#    day_min = 0
#    day_max = 0
#    dmin = 0
#    dmax = 0
#    day_year = 0
#    in_unit = 0
#    blk_unit = 0
#    ozone_unit = 0
#    oz_tab_unit = 0
#    bom_pres_unit = 0
#    gen_unit = 0
#    cal_unit = 0
#    cut_unit = 0
#    config_unit = 0
#    par_unit = 0
#    clock_unit = 0
#    pressure_option = 0
#    ozone_option = 0
#    num_pressure_points = 0
#    num_valid_pres_points = 0
#    cal_day = 0
#    valid_data = 0
#    time_corr = 0
#    qa_ref = 0
#    date_obs = [0] * 3
#    cal_epoch = [0] * 3
#    cal_instrument = 0
#    cal_model = 0
#    cal_langley_channels = 0
#    cal_num_langleys = 0
#    cal_num_general_cycles = 0
#    iwavelength = [0] * NumChannels
#    iorder = [0] * NumChannels
#    erms_in = [0.0] * NumChannels
#    wref_prior_gen_cals = [0.0] * NumChannels
#    lnv0coef_in = [[0.0] * (MaxPolyOrder + 1) for _ in range(NumChannels)]
#    lnv0coef = [[0.0] * (MaxPolyOrder + 1) for _ in range(NumChannels)]
#    elnv0coef = [[0.0] * (MaxPolyOrder + 1) for _ in range(NumChannels)]
#    wave_cal = [[0.0] * NumModels for _ in range(NumChannels)]
#    extension = ""
#    data_type = ""
#    site_code = ""
#    cal_in_ext = ""
#    cal_out_ext = ""
#    year_code = ""
#    month_code = ""
#    day_code = ""
#    data_root = ""
#    data_path = ""
#    data_suffix = ""
#    work_root = ""
#    work_path = ""
#    result_path = ""
#    config_root = ""
#    config_path = ""
#    pres_path = ""
#    bom_pres_path = ""
#    file_root = ""
#    config_file = ""
#    pres_file = ""
#    bom_pres_file = ""
#    in_file = ""
#    ozone_file = ""
#    blk_file = ""
#    gen_file = ""
#    cut_file = ""
#    cal_file = ""
#    par_file = ""
#    clock_file = ""
#    format_string = [""] * 10
#    config = {}
#    tz_hms = {}
#    obs_date = datetime.date(1, 1, 1)
#    date_time = {}
#    start_local_day = {}
#    end_local_day = {}
#    timezone_dt = {}
#    pdata = [{}] * MaxPressurePoints
#    sun3data = {}
#    sundata = {}
#    black = {}
#    sunpos = {}
#    cut_list = [{}] * MaxCutPeriods
#    num_cut = 0
#    jul_day = [0] * 3
#    ymdhms = [[0] * 7 for _ in range(3)]
#    station_lat_deg = 0.0
#    station_lon_deg = 0.0
#    station_lat_rad = 0.0
#    station_lon_rad = 0.0
#    solar_zenith = 0.0
#    solar_azimuth = 0.0
#    default_surface_pressure = 0.0
#    daily_mean_pressure = 0.0
#    default_ozone = 0.0
#    ozone_column = 0.0
#    lnv0ref = 0.0
#    lnv0ref1au = 0.0
#    sum_m = 0.0
#    sum_maod = 0.0
#    time_day = [[0.0] * 3 for _ in range(MaxTriplets)]
#    solar_zenith_deg = [[0.0] * 3 for _ in range(MaxTriplets)]
#    solar_zenith_app = [[0.0] * 3 for _ in range(MaxTriplets)]
#    solar_azimuth_deg = [[0.0] * 3 for _ in range(MaxTriplets)]
#    pressure = [0.0] * MaxTriplets
#    temperature = [0.0] * MaxTriplets
#    dsun_earth = [0.0] * MaxTriplets
#    air_mass_od_rayleigh = [[[0.0] * NumChannels for _ in range(3)] for _ in range(MaxTriplets)]
#    air_mass_od_aerosol = [[[0.0] * NumChannels for _ in range(3)] for _ in range(MaxTriplets)]
#    air_mass_od_ozone = [[[0.0] * NumChannels for _ in range(3)] for _ in range(MaxTriplets)]
#    air_mass = [[[0.0] * NumChannels for _ in range(3)] for _ in range(MaxTriplets)]
#    signal = [[[0.0] * NumChannels for _ in range(3)] for _ in range(MaxTriplets)]
#    ln_signal = [[[0.0] * NumChannels for _ in range(3)] for _ in range(MaxTriplets)]
#    triplet = [0.0] * 3
#    triplet_mean = [[0.0] * NumChannels for _ in range(MaxTriplets)]
#    triplet_sdev = [[0.0] * NumChannels for _ in range(MaxTriplets)]
#    triplet_cv = [[0.0] * NumChannels for _ in range(MaxTriplets)]
#    triplet_log = [[0.0] * NumChannels for _ in range(MaxTriplets)]
#    ratio = [[0.0] * NumChannels for _ in range(2)]
#    black_sun = [0.0] * NumChannels
#    black_skya = [0.0] * NumSkyA
#    black_skyk = [0.0] * NumSkyK
#    rayleigh_od = [0.0] * NumChannels
#    ozone_od = [0.0] * NumChannels
#    days_since_epoch = [0.0] * MaxLangleys
#    weight = [0.0] * MaxSinglets
#    lnv0glob = [[0.0] * NumChannels for _ in range(MaxLangleys)]
#    ratioglob = [[0.0] * NumChannels for _ in range(MaxLangleys)]
#    lnv0mean = [0.0] * NumChannels
#    lnv0sdev = [0.0] * NumChannels
#    ratiomean = [0.0] * NumChannels
#    ratiosdev = [0.0] * NumChannels
    file_ok = False
    there = False
    spread_flag = False
    fit_flag = False
#    max_sdev_fit = 0.0
#    max_rec_length = 0
#    iap = 0
#    num_langleys = 0
#    min_sun_triples = 10
#    nstart = [0] * 2
#    nend = [0] * 2
#    num_ok = [0] * 2
#    num_points = [0] * 2
#    ind_ok = [[0] * MaxSinglets for _ in range(2)]
#    num_sun_triples = [0] * 3
#    num_sun_singles = [0] * 3
#    x = [[0.0] * NumChannels for _ in range(MaxSinglets)]
#    y = [[0.0] * NumChannels for _ in range(MaxSinglets)]
#    z = [[0.0] * NumChannels for _ in range(MaxSinglets)]
#    v = [0.0] * MaxSinglets
#    w = [0.0] * MaxSinglets
#    intercept = [[0.0] * NumChannels for _ in range(2)]
#    del_intercept = [[0.0] * NumChannels for _ in range(2)]
#    slope = [[0.0] * NumChannels for _ in range(2)]
#    del_slope = [[0.0] * NumChannels for _ in range(2)]
    am_pm = ["Am", "Pm"]
#    pres_mean = [0.0] * 2
#    od_aerosol_ref = [0.0] * 2
#    current_month = 0
#    ozone = [0.0] * 32
#    residual = [0.0] * NumChannels
#    erms = [0.0] * NumChannels

    # Read input parameters and initialize variables
    # ...

    # Main processing loop
CurrentMonth = 0
for Year in range(YearMin, YearMax+1):
    Year_Code = Integer_to_String(Year, 4)
    MMin = 1
    MMax = 12
    if Year == YearMin:
        MMin = MonthMin
    if Year == YearMax:
        MMax = MonthMax
    for Month in range(MMin, MMax+1):
        Month_Code = Integer_to_String(Month, 2)
        Dmin = 1
        Dmax = Days_in_Month(Month, Year)
        if Month == MonthMin and Year == YearMin:
            Dmin = DayMin
        if Month == MonthMax and Year == YearMax:
            Dmax = DayMax
        for Day in range(Dmin, Dmax+1):
            Day_Code = Integer_to_String(Day, 2)

            print('*************************************')
            print(f' Instrument #{Instrument}, {Year_Code}{Month_Code}{Day_Code}')
            DataPath = f"{DataRoot}{DataSuffix}/{Config.Id}/#" \
                       f"{Integer_to_String(Instrument, 2)}/" \
                       f"{Year_Code}/{Month_Code}/{Day_Code}/"

            DateObs[0] = Year
            DateObs[1] = Month
            DateObs[2] = Day
            CalDay = DayCount(CalEpoch, DateObs)
            LnV0Ref1AU = 0
            for I in range(Iorder[Iref]+1):
                LnV0Ref1AU += LnV0Coef_In[Iref][I] * CalDay**I

            ObsDate = Date(Year, Month, Day)
            if Month != CurrentMonth:
                Ozone = Read_BoM_Ozone_2011(Dbug[0], OzoneUnit, ObsDate)
                CurrentMonth = Month

            FileRoot = f"{Year_Code[2:]}_{Month_Code}_{Day_Code}"

            InFile = Make_File_Name(DataPath, FileRoot, Extension)
            print(f"InFile={InFile}")

            FileOK = os.path.isfile(InFile)
            if FileOK:
                InUnit = Next_Unit_Number()
                with open(InFile, 'r') as file:
                    if PressureOption == 1:
                        PresPath = DataPath
                        if Instrument == 2 and Year == 1998 and 4 <= Month <= 7:
                            i = PresPath.index('tta')
                            PresPath = f"{PresPath[:i]}ttb/#03{PresPath[i+7:]}"
                        PresFile = Make_File_Name(PresPath, FileRoot, 'hpa')
                        print(f"Pressure File Name={PresFile}")
                        FileOK = os.path.isfile(PresFile)
                        if FileOK:
                            NumPressurePoints, NumValidPresPoints, PData, DailyMeanPressure = \
                                Read_Pressure_File(Dbug[0], PresFile)
                            if NumValidPresPoints <= MinValidPresPoints:
                                print('Insufficient valid pressure data in file')
                        else:
                            NumPressurePoints = 0
                            print('Pressure file not found (.hpa)-has envcal been run?')
                            HALT()
                        if NumPressurePoints <= 0:
                            print(f"Assuming Pdefault ={DefaultSurfacePressure}")
                    elif PressureOption == 3:
                        NumPressurePoints, NumValidPresPoints, PData, DailyMeanPressure = \
                            Read_BoM_Pressure(Dbug[0], BoMPresUnit, ObsDate)

                if OzoneOption == 0:
                    OzoneColumn = DefaultOzone
                else:
                    if OzoneOption == 1:
                        if Ozone[Day] > 0:
                            OzoneColumn = Ozone[Day] / 1000.0
                        else:
                            print(f"Ozone file: Missing day, using monthly mean={Ozone[0]}")
                            OzoneColumn = Ozone[0] / 1000.0
                    elif OzoneOption == 2:
                        OzoneColumn = Ozone[0] / 1000.0

                OzoneOD[:NumChannels] = [OzoneCoef[i][Model] * OzoneColumn for i in range(NumChannels)]
                print(f"Ozone(atm-cm), od(670): {OzoneColumn}, {OzoneOD[I670]}")

                if ApplyClockFix:
                    TimeCorr = Get_Time_Correction(Dbug[0], ClockUnit, ObsDate)

                BlkFile = Make_File_Name(DataPath, FileRoot, 'blk')
                print(f"BlackFile={BlkFile}")
                FileOK = os.path.isfile(BlkFile)

                if FileOK:
                    BlkUnit = Next_Unit_Number()
                    with open(BlkFile, 'r') as file:
                        I = 0
                        Eof = 0
                        BlackSun[:NumChannels] = [0] * NumChannels
                        BlackSkyA[:NumSkyA] = [0] * NumSkyA
                        BlackSkyK[:NumSkyK] = [0] * NumSkyK
                        while Eof == 0:
                            DateTime, Black, Eof = Read_Black_Record(Dbug[0], Model, file)
                            if Eof == 0:
                                I += 1
                                BlackSun[:NumChannels] = [BlackSun[j] + Black.Sun.Chan[j] for j in range(NumChannels)]
                                BlackSkyA[:NumSkyA] = [BlackSkyA[j] + Black.SkyA[j] for j in range(NumSkyA)]
                                BlackSkyK[:NumSkyK] = [BlackSkyK[j] + Black.SkyK[j] for j in range(NumSkyK)]
                        BlackSun[:NumChannels] = [BlackSun[j] / I for j in range(NumChannels)]
                        BlackSkyA[:NumSkyA] = [BlackSkyA[j] / I for j in range(NumSkyA)]
                        BlackSkyK[:NumSkyK] = [BlackSkyK[j] / I for j in range(NumSkyK)]
                        if Dbug[0]:
                            print(f"Number of Black records: {I}")
                            for i in range(NumChannels):
                                print(f"Mean Black Sun     : {BlackSun[i]}")
                            HALT()
                            for i in range(NumSkyA):
                                print(f"Mean Black SkyA     : {BlackSkyA[i]}")
                            HALT()
                            for i in range(NumSkyK):
                                print(f"Mean Black SkyK     : {BlackSkyK[i]}")
                            HALT()
                else:
                    BlackSun[:NumChannels] = [0] * NumChannels
                    BlackSkyA[:NumSkyA] = [0] * NumSkyA
                    BlackSkyK[:NumSkyK] = [0] * NumSkyK

                Eof = 0
                NumSunTriples[:3] = [0] * 3
                K = 0
                NumSunSingles[:3] = [0] * 3
                L = 0
                print(f"Processing {DataType} records...")

                if DataType == 'sun':
                    DateTime, Sun3Data, ValidData, Eof = Read_Triple_Sun_Record(Dbug[1], InUnit, Model)
                elif DataType == 'lsu':
                    DateTime, SunData, ValidData, Eof = Read_Single_Sun_Record(Dbug[1], InUnit, Model)
                else:
                    print(f" Data type {DataType} not supported yet")

                while Eof == 0:
                    StartLocalDay = Date_Time(Year, Month, Day, 4, 0, 0) + TimeZone_DT
                    EndLocalDay = Date_Time(Year, Month, Day, 21, 0, 0) + TimeZone_DT
                    while (DateTime < StartLocalDay or DateTime > EndLocalDay) and Eof == 0:
                        print('Time outside local day, skipping')
                        print(f"Startday: {StartLocalDay}")
                        print(f"Endday  : {EndLocalDay}")
                        print(f"ObsDT   : {DateTime}")
                        if DataType == 'sun':
                            DateTime, Sun3Data, ValidData, Eof = Read_Triple_Sun_Record(Dbug[1], InUnit, Model)
                        elif DataType == 'lsu':
                            DateTime, SunData, ValidData, Eof = Read_Single_Sun_Record(Dbug[1], InUnit, Model)
                        else:
                            print(f" Data type {DataType} not supported yet")
                    if Eof != 0:
                        break
                    if DataType == 'sun':
                        YMDHMS[0][0] = DateTime.Year
                        YMDHMS[1][0] = DateTime.Month
                        YMDHMS[2][0] = DateTime.Day
                        YMDHMS[3][0] = DateTime.Hour
                        YMDHMS[4][0] = DateTime.Minute
                        YMDHMS[5][0] = DateTime.Second
                        YMDHMS[6][0] = 0
                        if ApplyClockFix:
                            Add_Time(YMDHMS[0], -TimeCorr)
                        K += 1
                        NumSunTriples[0] = K
                        JulDay[0], TimeDay[K-1][0] = Julian(YMDHMS[:, 0])
                        SunPos = Solar_Position_Almanac(JulDay[0], TimeDay[K-1][0])
                        DSunEarth[K-1] = SunPos.DsunEarth
                        if Dbug[1]:
                            print(f"Year             : {YMDHMS[0][0]}")
                            print(f"Month            : {YMDHMS[1][0]}")
                            print(f"Day              : {YMDHMS[2][0]}")
                            print(f"Hour             : {YMDHMS[3][0]}")
                            print(f"Minute           : {YMDHMS[4][0]}")
                            print(f"Second           : {YMDHMS[5][0]}")
                            print(f"Julian Day       : {JulDay[0]}")
                            print(f"TimeDay          : {TimeDay[K-1][0]}")
                            print(f"DSunEarth        : {DSunEarth[K-1]}")
                            HALT()
                        I = 0
                        NewStation = True
                        SolarZenith, SolarAzimuth = Satellite_Position(JulDay[I], TimeDay[K-1][I], NewStation,
                                                                       StationLatRad, StationLonRad)
                        SolarZenithDeg[K-1][I] = SolarZenith * DegreesPerRadian
                        SolarZenithApp[K-1][I] = Apparent_Zenith(SolarZenithDeg[K-1][I])
                        SolarAzimuthDeg[K-1][I] = SolarAzimuth * DegreesPerRadian
                        if Dbug[1]:
                            print(f"Latitude      : {StationLatDeg}")
                            print(f"Longitude     : {StationLonDeg}")
                            print(f"Sol Zen (true): {SolarZenithDeg[K-1][I]}")
                            print(f"Sol Zen (app) : {SolarZenithApp[K-1][I]}")
                            print(f"Solar Azimuth : {SolarAzimuthDeg[K-1][I]}")
                            HALT()

                        for I in range(1, 3):
                            YMDHMS[:, I] = YMDHMS[:, I-1]
                            Add_Time(YMDHMS[0, I], 30)
                            JulDay[I], TimeDay[K-1][I] = Julian(YMDHMS[:, I])
                            if Dbug[1]:
                                print(f"Day                 : {DateTime.Day}")
                                print(f"Month               : {DateTime.Month}")
                                print(f"Year                : {DateTime.Year}")
                                print(f"Julian Day number   : {JulDay[I]}")
                                print(f"Fractional Day      : {TimeDay[K-1][I]}")
                                HALT()
                            NewStation = True
                            SolarZenith, SolarAzimuth = Satellite_Position(JulDay[I], TimeDay[K-1][I], NewStation,
                                                                           StationLatRad, StationLonRad)
                            SolarZenithDeg[K-1][I] = SolarZenith * DegreesPerRadian
                            SolarZenithApp[K-1][I] = Apparent_Zenith(SolarZenithDeg[K-1][I])
                            SolarAzimuthDeg[K-1][I] = SolarAzimuth * DegreesPerRadian
                            if Dbug[1]:
                                print(f"Latitude      : {StationLatDeg}")
                                print(f"Longitude     : {StationLonDeg}")
                                print(f"Sol Zen (true): {SolarZenithDeg[K-1][I]}")
                                print(f"Sol Zen (app) : {SolarZenithApp[K-1][I]}")
                                print(f"Solar Azimuth : {SolarAzimuthDeg[K-1][I]}")
                                HALT()

                        if SolarAzimuthDeg[K-1][1] < 180.0:
                            NumSunTriples[1] += 1
                            NumSunSingles[1] += 3
                        else:
                            NumSunTriples[2] += 1
                            NumSunSingles[2] += 3

                        Pressure[K-1] = GetPressure(DateTime, PData, PressureOption, NumPressurePoints,
                                                    NumValidPresPoints, DailyMeanPressure, DefaultSurfacePressure)

                        for N in range(NumChannels):
                            Rayleigh(Dbug[2], Wavelength[N][Model], Pressure[K-1], RayleighOD[N])
                            if Dbug[2]:
                                print(f"Wavelength: {Wavelength[N][Model]}")
                                print(f"Pressure  : {Pressure[K-1]}")
                                print(f"RayleighOD: {RayleighOD[N]}")

                        for N in range(NumChannels):
                            for I in range(3):
                                AirMassODRayleigh[K-1][I][N] = GetAirMass(1, SolarZenithApp[K-1][I]) * RayleighOD[N]
                                AirMassODAerosol[K-1][I][N] = GetAirMass(2, SolarZenithApp[K-1][I]) * 0.03
                                AirMassODOzone[K-1][I][N] = GetAirMass(3, SolarZenithApp[K-1][I]) * OzoneOD[N]
                                AirMass[K-1][I][N] = (AirMassODRayleigh[K-1][I][N] +
                                                      AirMassODAerosol[K-1][I][N] +
                                                      AirMassODOzone[K-1][I][N]) / (RayleighOD[N] + 0.03 + OzoneOD[N])
                                if Dbug[2]:
                                    print(f"Time          : {YMDHMS[:6, I]}")
                                    print(f"Sol Zen (true): {SolarZenithDeg[K-1][I]}")
                                    print(f"Sol Zen (app) : {SolarZenithApp[K-1][I]}")
                                    print(f"Rayleigh am   : {GetAirMass(1, SolarZenithApp[K-1][I])}")
                                    print(f"Aerosol  am   : {GetAirMass(2, SolarZenithApp[K-1][I])}")
                                    print(f"Ozone    am   : {GetAirMass(3, SolarZenithApp[K-1][I])}")
                                    print(f"Weighted am   : {AirMass[K-1][I][N]}")
                                    HALT()

                        for N in range(NumChannels):
                            Triplet[:] = [Sun3Data.Signal[i].Chan[N] - BlackSun[N] for i in range(3)]
                            TripletMean[K-1][N], TripletSdev[K-1][N] = Stat(Triplet, 3)
                            if TripletSdev[K-1][N] >= 1000.0:
                                TripletSdev[K-1][N] = 999.99
                            if TripletMean[K-1][N] > 0.0:
                                TripletCv[K-1][N] = 100 * TripletSdev[K-1][N] / TripletMean[K-1][N]
                            else:
                                TripletCv[K-1][N] = 99.99
                            if TripletCv[K-1][N] >= 99.99:
                                TripletCv[K-1][N] = 99.99
                            if TripletMean[K-1][N] > 0.0:
                                TripletLog[K-1][N] = math.log(TripletMean[K-1][N])
                            else:
                                TripletLog[K-1][N] = -9.99
                            if Dbug[1]:
                                print(f"In Channel loop, Cha #= {N}")
                                print(f"Triplet(1)= {Triplet[0]}")
                                print(f"TripletMean= {TripletMean[K-1][N]}")
                                print(f"TripletSdev= {TripletSdev[K-1][N]}")
                                HALT()

                        for I in range(3):
                            L += 1
                            if L > MaxTriplets:
                                print("Array bounds exceeded")
                                print("Increase MaxTriplets!")
                                HALT()
                            NumSunSingles[0] = L
                            for N in range(NumChannels):
                                Signal[K-1][I][N] = Sun3Data.Signal[I].Chan[N] - BlackSun[N]
                                if Signal[K-1][I][N] > 0.0:
                                    LnSignal[K-1][I][N] = math.log(Signal[K-1][I][N])
                                else:
                                    LnSignal[K-1][I][N] = -9.99

                        Temperature[K-1] = Sun3Data.Temperature
                    else:
                        print(f" Data type {DataType} not supported yet")

                    if DataType == 'sun':
                        DateTime, Sun3Data, ValidData, Eof = Read_Triple_Sun_Record(Dbug[1], InUnit, Model)
                    elif DataType == 'lsu':
                        DateTime, SunData, ValidData, Eof = Read_Single_Sun_Record(Dbug[1], InUnit, Model)
                    else:
                        print(f" Data type {DataType} not supported yet")

                LangFlag[:] = [False] * 2
                if DataType == 'sun' or DataType == 'lsu':
                    print('General method, processing...')
                    for Iap in range(2):
                        Nstart[Iap] = 1 + NumSunTriples[1] * Iap
                        Nend[Iap] = Nstart[Iap] + NumSunTriples[Iap] - 1

                    LnV0Ref = LnV0Ref1AU - 2.0 * math.log(DSunEarth[Nstart[0]-1])
                    if Dbug[2]:
                        print(f"Number of sun triples(am)  : {NumSunTriples[1]}")
                        print(f"Number of sun triples(pm)  : {NumSunTriples[2]}")
                        print(f"Number of sun triples(tot) : {NumSunTriples[0]}")
                        print(f"Number of sun singles(am)  : {NumSunSingles[1]}")
                        print(f"Number of sun singles(pm)  : {NumSunSingles[2]}")
                        print(f"Number of sun singles(tot) : {NumSunSingles[0]}")
                        print(f"Nstart (am)                : {Nstart[0]}")
                        print(f"Nend   (am)                : {Nend[0]}")
                        print(f"Nstart (pm)                : {Nstart[1]}")
                        print(f"Nend   (pm)                : {Nend[1]}")
                        HALT()

                    for Iap in range(2):
                        ObsDate = Date(Year, Month, Day)
                        Include = True
                        for Ic in range(NumCut):
                            if ObsDate == CutList[Ic].Date and AmPm[Iap] == CutList[Ic].AP:
                                Include = False
                                break
                        if Include:
                            if NumSunTriples[Iap] >= MinSunTriples:
                                PresMean[Iap] = 0
                                N = 0
                                for K in range(Nstart[Iap]-1, Nend[Iap]):
                                    if Pressure[K] > 0.0:
                                        N += 1
                                        PresMean[Iap] += Pressure[K]
                                PresMean[Iap] /= N
                                SpreadFlag, NumOK[Iap], IndOK[:, Iap] = CheckTripletCv(Dbug[2], NumSunTriples[Iap],
                                                                                        Nstart[Iap]-1,
                                                                                        AirMass[:, 1, I870],
                                                                                        TripletCv[:, I870])
                                if NumOK[Iap] >= MinSunTriples and SpreadFlag:
                                    print(" Passed triplet cv langley filter")
                                    N = NumChannels
                                    I = 0
                                    for K in range(NumOK[Iap]):
                                        for J in range(3):
                                            I += 1
                                            X[I-1, :N] = AirMass[IndOK[K, Iap]-1, J, :N]
                                            Y[I-1, :N] = LnSignal[IndOK[K, Iap]-1, J, :N]
                                            V[I-1] = GetAirMass(2, SolarZenithApp[IndOK[K, Iap]-1, J])
                                            W[I-1] = (LnV0Ref -
                                                      LnSignal[IndOK[K, Iap]-1, J, Iref] -
                                                      AirMassODRayleigh[IndOK[K, Iap]-1, J, Iref] -
                                                      AirMassODOzone[IndOK[K, Iap]-1, J, Iref])
                                            Z[I-1, :N] = (LnSignal[IndOK[K, Iap]-1, J, :N] +
                                                          AirMassODRayleigh[IndOK[K, Iap]-1, J, :N] +
                                                          AirMassODOzone[IndOK[K, Iap]-1, J, :N])
                                    NumPoints[Iap] = I

                                    if UseRefChannel4QA:
                                        QARef = IRef
                                    else:
                                        QARef = I870
                                    SpreadFlag, FitFlag, IndOK[:, Iap] = CheckFitQuality(Dbug[2], NumPoints[Iap],
                                                                                         X, Y, QARef, MaxSdevFit)
                                    N = NumChannels
                                    for I in range(NumPoints[Iap]):
                                        V[I] = V[IndOK[I, Iap]-1]
                                        W[I] = W[IndOK[I, Iap]-1]
                                        Z[I, :N] = Z[IndOK[I, Iap]-1, :N]

                                    if SpreadFlag and FitFlag:
                                        LangFlag[Iap] = True
                                        print(" Passed fit quality langley filter")
                                        NumLangleys += 1
                                        SumM = 0
                                        SumMAod = 0
                                        for I in range(NumPoints[Iap]):
                                            SumMAod += W[I]
                                            SumM += V[I]
                                        ODAerosolRef[Iap] = SumMAod / SumM
                                        Weight[:NumPoints[Iap]] = [1.0] * NumPoints[Iap]
                                        for N in range(NumLangleyChannels[Model]):
                                            Intercept[Iap, N], Slope[Iap, N], Residual[N], Erms[N], \
                                            DelIntercept[Iap, N], DelSlope[Iap, N] = ELFIT(NumPoints[Iap],
                                                                                           Weight[:NumPoints[Iap]],
                                                                                           W[:NumPoints[Iap]],
                                                                                           Z[:NumPoints[Iap], N])
                                            Intercept[Iap, N] += 2.0 * math.log(DSunEarth[Nstart[0]-1])
                                            Ratio[Iap, N] = -Slope[Iap, N]
                                            LnV0Glob[NumLangleys-1, N] = Intercept[Iap, N]
                                            RatioGlob[NumLangleys-1, N] = Ratio[Iap, N]
                                        DaysSinceEpoch[NumLangleys-1] = CalDay
                                    else:
                                        print(" Failed fit quality filter")
                                else:
                                    print(" Failed triplet cv filter")
                            else:
                                print(" Insufficient triplets")

                if LangFlag[0] or LangFlag[1]:
                    DayYear = Day_of_Year(Day, Month, Year)
                    for Iap in range(2):
                        if LangFlag[Iap]:
                            print(GenUnit, FormatString[1],
                                  Year, Month, Day, CalDay, AmPm[Iap],
                                  PresMean[Iap], NumOK[Iap],
                                  Iref, LnV0Ref1AU,
                                  ODAerosolRef[Iap],
                                  [(Intercept[Iap, N], DelIntercept[Iap, N],
                                    -Slope[Iap, N]) for N in range(NumLangleyChannels[Model])])

            else:
                print(f"InFile {InFile} not found")

    if NumLangleys >= 2:
        for N in range(NumLangleyChannels[Model]):
            LnV0Mean[N], LnV0Sdev[N] = Stat(LnV0Glob[:NumLangleys, N], NumLangleys)
            RatioMean[N], RatioSdev[N] = Stat(RatioGlob[:NumLangleys, N], NumLangleys)
    elif NumLangleys == 1:
        for N in range(NumLangleyChannels[Model]):
            LnV0Mean[N] = LnV0Glob[0, N]
            RatioMean[N] = RatioGlob[0, N]

    print(GenUnit, FormatString[2],
          "#Summary", NumLangleys,
          [(LnV0Mean[N], LnV0Sdev[N], RatioMean[N]) for N in range(NumLangleyChannels[Model])])

    if Iorder[I870] == 1:
        Weight[:NumLangleys] = [1.0] * NumLangleys
        for N in range(NumLangleyChannels[Model]):
            LnV0Coef[0, N], LnV0Coef[1, N], Residual[N], Erms[N], \
            ELnV0Coef[0, N], ELnV0Coef[1, N] = ELFIT(NumLangleys, Weight[:NumLangleys],
                                                     DaysSinceEpoch[:NumLangleys], LnV0Glob[:NumLangleys, N])
            if Dbug[3]:
                print(f"Channel  number  {N}")
                print(f"NumLangleys      {NumLangleys}")
                print(f"Days since epoch {DaysSinceEpoch[:NumLangleys]}")
                print(f"LnV0             {LnV0Glob[:NumLangleys, N]}")
                print(f"Intercept        {LnV0Coef[0, N]}")
                print(f"DelIntercept     {ELnV0Coef[0, N]}")
                print(f"Slope            {LnV0Coef[1, N]}")
                print(f"DelSlope         {ELnV0Coef[1, N]}")
                print(f"Rms error        {Erms[N]}")
                HALT()

        LnV0Coef[0, IWV] = LnV0Coef_In[0, IWV]
        LnV0Coef[1, IWV] = LnV0Coef_In[1, IWV]
        Erms[IWV] = Erms_In[IWV]
    elif NumLangleys >= 1:
        for N in range(NumLangleyChannels[Model]):
            LnV0Coef[0, N] = LnV0Mean[N]
            Erms[N] = LnV0Sdev[N]
        LnV0Coef[0, IWV] = LnV0Coef_In[0, IWV]
        Erms[IWV] = Erms_In[IWV]

    Erms[IRef] = Erms_In[IRef]

    if MakeCal:
        There = os.path.isfile(CalFile)
        if There:
            print(f"{CalFile} exists, will be overwritten.")
            HALT()
        CalUnit = Next_Unit_Number()
        with open(CalFile, 'w') as file:
            file.write(f"# Calibration file generated from General method between:\n")
            file.write(f"# {YearMin:04d}{MonthMin:02d}{DayMin:02d} and\n")
            file.write(f"# {YearMax:04d}{MonthMax:02d}{DayMax:02d}.\n")
            file.write(f"# {Instrument:05d}        -- Instrument number\n")
            file.write(f"# {Model:05d}        -- Model number\n")
            file.write(f"# {NumLangleyChannels[Model]:05d}        -- Number of Langley wavelengths\n")
            file.write(f"# {NumLangleys:05d}        -- Number of Langley intervals in period\n")
            file.write(f"# {CalNumGeneralCycles+1:05d}        -- Number of General Method cycles applied\n")
            file.write(f"# {CalEpoch[0]:05d}{CalEpoch[1]:02d}{CalEpoch[2]:02d}  -- Calibration epoch\n")
            FormatString[4] = f"# Wavel(nm) Order {' '.join([f'LnV0({i})' for i in range(Iorder[IRef]+1)])} Erms WlRef\n"
            file.write(FormatString[4])
            for N in range(NumLangleyChannels[Model]):
                file.write(f"{Wavelength[N][Model]:10.1f} {Iorder[N]:6d} {' '.join([f'{LnV0Coef[I, N]:13.5e}' for I in range(Iorder[N]+1)])} {Erms[N]:13.5e} {WRefPriorGenCals[:CalNumGeneralCycles]} {Wavelength[IRef][Model]:13.5e}\n")
            file.write(f"{Wavelength[IWV][Model]:10.1f} {Iorder[IWV]:6d} {' '.join([f'{LnV0Coef[I, IWV]:13.5e}' for I in range(Iorder[IWV]+1)])} {Erms[IWV]:13.5e}\n")

    print("Done.")

if __name__ == "__main__":
    General_Calibration()
