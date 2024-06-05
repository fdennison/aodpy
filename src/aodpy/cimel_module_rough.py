import math

class Env_Stat:
    def __init__(self):
        self.Mean = 0
        self.Min = 0
        self.Max = 0
        self.Sdev = 0

class Real_Stat:
    def __init__(self):
        self.Mean = 0.0
        self.Min = 0.0
        self.Max = 0.0
        self.Sdev = 0.0

class Site_Configuration:
    def __init__(self):
        self.Name = ""
        self.BoMNumber = 0
        self.Id = ""
        self.Latitude = 0.0
        self.Longitude = 0.0
        self.Height = 0.0
        self.TimeZone = 0.0
        self.Date = [Date(), Date()]
        self.TOffset = Time()
        self.CimelNumber = 0
        self.CimelModel = 0
        self.Version = 0
        self.Card = 0
        self.Baro = 0
        self.Wind = 0
        self.Neph = ""
        self.CimCalFile = ""
        self.SkyCalDate = Date()

class CutCalPeriod:
    def __init__(self):
        self.Date = Date()
        self.AP = ""

class Status_Record:
    def __init__(self):
        self.Hour = 0
        self.Vext = 0.0
        self.Vint = 0.0
        self.Text = 0.0
        self.Tint = 0.0
        self.Vbaro = 0.0
        self.NumLogSamples = 0

class Status_:
    def __init__(self):
        self.ID = 0
        self.DateTime = Date_Time()
        self.Data = [Status_Record() for _ in range(NumLogTimes)]

class Raw_Environment_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.Pressure = Env_Stat()
        self.Temperature = Env_Stat()
        self.WindV = Env_Stat()
        self.WindD = Env_Stat()
        self.VBatt = 0

class Environment_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.Pressure = Real_Stat()
        self.Temperature = Real_Stat()
        self.WindV = Real_Stat()
        self.WindD = Real_Stat()
        self.VBatt = 0.0

class Pressure_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.Pressure = 0.0

class BoM_Pressure_Record:
    def __init__(self):
        self.Id = ""
        self.Station = 0
        self.DateTime = Date_Time()
        self.MSLP = 0.0
        self.Qual_MSLP = ""
        self.Station_Level_Pressure = 0.0
        self.Qual_SLP = ""
        self.QNHP = 0.0
        self.Qual_QNHP = ""
        self.AWS_Flag = ""
        self.EOR_Hash = ""

class BoM_Ozone_Record:
    def __init__(self):
        self.Station = 0
        self.LineNumber = 0
        self.Month = 0
        self.Year = 0
        self.StartHourUT = [0] * 32
        self.EndHourUT = [0] * 32
        self.WavelengthCode = [0] * 32
        self.ObsCode = [0] * 32
        self.DU = [0] * 32

class BoM_Ozone_Record_2011:
    def __init__(self):
        self.Station = 0
        self.Year = 0
        self.Month = 0
        self.Day = 0
        self.StartHourUT = 0
        self.EndHourUT = 0
        self.WavelengthCode = 0
        self.ObsCode = 0
        self.DU = 0
        self.StdErrDU = 0
        self.Instrument = 0
        self.SerialNumber = 0

class BoM_Transmission_Record:
    def __init__(self):
        self.Date = Date()
        self.DayOfYear = 0
        self.Minute = 0.0
        self.Day = 0.0
        self.Pressure = 0.0
        self.SolarZD = 0.0
        self.Airmass = 0.0
        self.FLT = 0
        self.StDev = 0.0
        self.Tran = [0.0] * MaxBomChannels
        self.Diffuse = [0.0] * MaxBomChannels

class Neph_Record:
    def __init__(self):
        self.RecordNumber = 0
        self.DateTime = Date_Time()
        self.ScatCoef = 0.0
        self.CalCoef = 0.0
        self.Pressure = 0.0
        self.Temperature = 0.0
        self.RH = 0.0

class Neph_Cal_Record:
    def __init__(self):
        self.Date = Date()
        self.Gain = 0.0
        self.Offset = 0.0

class Aod_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.SunZen = 0.0
        self.AirMass = 0.0  # 870nm
        self.Pressure = 0.0
        self.Temperature = 0.0
        self.Aod = [0.0] * NumChannels
        self.StdUAod = [0.0] * NumChannels
        self.Wvap = 0.0
        self.StdUWvap = 0.0
        self.Aod1020c = 0.0
        self.Angstrom46 = 0.0
        self.StdUAngstrom46 = 0.0
        self.Angstrom48 = 0.0
        self.StdUAngstrom48 = 0.0
        self.Angstrom68 = 0.0
        self.StdUAngstrom68 = 0.0
        self.Tran = [0.0] * NumChannels
        self.CirrusFlag = 0

class Cod_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.SunZen = 0.0
        self.AirMass = 0.0  # 870nm
        self.Pressure = 0.0
        self.Temperature = 0.0
        self.Aod = [0.0] * NumChannels
        self.Wvap = 0.0
        self.Aod1020c = 0.0
        self.Angstrom46 = 0.0
        self.Angstrom48 = 0.0
        self.Angstrom68 = 0.0
        self.Aod_500_Fine = 0.0
        self.Aod_500_Coarse = 0.0
        self.SeparationError = 0.0
        self.CirrusFlag = 0

class MPL_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.LayerDepthkm = 0.0
        self.CloudBackgroundRatio = 0.0

class Sun_Record:
    def __init__(self):
        self.Chan = [0] * NumChannels

class Single_Sun_Record:
    def __init__(self):
        self.Signal = Sun_Record()
        self.Temperature = 0.0

class Triple_Sun_Record:
    def __init__(self):
        self.Signal = [Sun_Record() for _ in range(3)]
        self.Temperature = 0.0

class Almucantar_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.AerosolChannel = 0
        self.Signal = [0] * NumAzimuths
        self.Temperature = 0.0

class Almucantar_Split_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.AerosolChannel = 0
        self.SignalRight = [0] * NumNAzimuths
        self.SignalLeft = [0] * NumNAzimuths
        self.Temperature = 0.0

class PP_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.AerosolChannel = 0
        self.Signal = [0] * NumPPZeniths
        self.Temperature = 0.0

class PPP_Record:
    def __init__(self):
        self.DateTime = Date_Time()
        self.Signal = [[0] * NumPPPZeniths for _ in range(3)]
        self.Temperature = 0.0

class Black_Record:
    def __init__(self):
        self.Sun = Sun_Record()
        self.SkyA = [0] * NumSkyA
        self.SkyK = [0] * NumSkyK

class Sky_Record:
    def __init__(self):
        self.SkyA = [0] * NumSkyA
        self.SkyK = [0] * NumSkyK
        self.Temperature = 0.0

class Sphere_Cal:
    def __init__(self):
        self.Instrument = 0
        self.DateTime = Date_Time()
        self.Waveff = [0] * NumSkyA
        self.GainA = [0.0] * NumSkyA
        self.SdevA = [0.0] * NumSkyA
        self.GainK = [0.0] * NumSkyK
        self.SdevK = [0.0] * NumSkyK

# Constants
MaxTriplets = 1024
MaxSinglets = 3 * MaxTriplets
MaxCutPeriods = 512
NumLogTimes = 24
MaxPressurePoints = 512
MinValidPresPoints = 3
MinSignal = 50
MaxSignal = 65000
NumModels = 8
NumChannels = 10
NumAerosolChannels = 4
NumLangleyChannels = [4, 4, 7, 7, 9, 9, 7, 8]
I500 = [7, 6, 7, 7, 6, 6]
Wavelength_Aerosol = [440.0, 670.0, 870.0, 1020.0]
Wavelength = [
    [440.0, 670.0, 870.0, 1020.0, 0.0, 0.0, 0.0, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 870.1, 870.2, 870.3, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 340.0, 380.0, 500.0, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 380.0, 500.0, 778.0, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 340.0, 380.0, 500.0, 1021.0, 1246.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 340.0, 380.0, 500.0, 1021.0, 1640.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 412.0, 500.0, 550.0, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 412.0, 500.0, 531.0, 550.0, 0.0, 936.0]
]
SkyChanIndex = [
    [4, 3, 2, 1, 0, 0, 0, 0, 0, 0],
    [4, 0, 2, 1, 0, 3, 0, 0, 0, 0],
    [4, 3, 2, 1, 0, 0, 0, 0, 0, 0],
    [4, 3, 2, 1, 0, 0, 0, 0, 0, 0],
    [4, 3, 2, 1, 0, 0, 0, 0, 0, 0],
    [4, 3, 2, 1, 0, 0, 0, 0, 0, 0],
    [4, 3, 2, 1, 0, 0, 0, 0, 0, 0],
    [4, 3, 2, 1, 0, 0, 0, 0, 0, 0]
]
OzoneCoef = [
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0447, 0.0002, 0.0327, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0004, 0.0322, 0.0087, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0447, 0.0002, 0.0327, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0447, 0.0002, 0.0327, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0002, 0.0327, 0.0327, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0002, 0.0327, 0.0327, 0.0327, 0.0000, 0.0000]
]
NumAzimuths = 76
AlmAzimuth = [0.0, -6.0, -5.0, -4.0, -3.5, -3.0, -2.5, -2.0, 2.0, 2.5,
              3.0, 3.5, 4.0, 5.0, 6.0, 6.0, 7.0, 8.0, 10.0, 12.0,
              14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
              60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0,
              220.0, 240.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 315.0, 320.0,
              325.0, 330.0, 335.0, 340.0, 342.0, 344.0, 346.0, 348.0, 350.0, 352.0,
              353.0, 354.0, 354.0, 355.0, 356.0, 356.5, 357.0, 357.5, 358.0, 362.0,
              362.5, 363.0, 363.5, 364.0, 365.0, 366.0]

CollGainAlm = ["SunL", "SunH", "SunH", "SunH", "SunH", "SunH", "SunH", "SunH", "SunH", "SunH",
               "SunH", "SunH", "SunH", "SunH", "SunH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
               "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
               "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
               "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
               "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
               "SkyH", "SkyH", "SunH", "SunH", "SunH", "SunH", "SunH", "SunH", "SunH", "SunH",
               "SunH", "SunH", "SunH", "SunH", "SunH", "SunH"]

NumNAzimuths = 30

AlmRightAzimuth = [0.0,
                   3.0, 3.5, 4.0, 5.0, 6.0, 6.0, 7.0, 8.0, 10.0, 12.0,
                   14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
                   60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0, 160.0, 180.0]

AlmLeftAzimuth = [0.0,
                  357.0, 356.5, 356.0, 355.0, 354.0, 354.0, 353.0, 352.0, 350.0, 348.0,
                  346.0, 344.0, 342.0, 340.0, 335.0, 330.0, 325.0, 320.0, 315.0, 310.0,
                  300.0, 290.0, 280.0, 270.0, 260.0, 240.0, 220.0, 200.0, 180.0]

CollGainAlmN = ["SunL", "SunH", "SunH", "SunH", "SunH", "SunH",
                "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
                "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
                "SkyH", "SkyH", "SkyH", "SkyH"]

NumPPZeniths = 40

PPZenith = [-6.0, -5.0, -4.0, -3.5, -3.0, -2.5, -2.0, 0.0, 2.0, 2.5,
            3.0, 3.5, 4.0, 5.0, 6.0, 6.0, 8.0, 10.0, 12.0, 14.0,
            16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0,
            65.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0]

CollGainPP = ["SunH", "SunH", "SunH", "SunH", "SunH", "SunH", "SunH", "SunL", "SunH", "SunH",
              "SunH", "SunH", "SunH", "SunH", "SunH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
              "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
              "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH"]

NumPPPZeniths = 35

PPPZenith = [95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0,
             145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0,
             195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0,
             245.0, 250.0, 255.0, 260.0, 265.0]

CollGainPPP = ["SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
               "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
               "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH", "SkyH",
               "SkyH", "SkyH", "SkyH", "SkyH", "SkyH"]

MaxBomChannels = 8
MaxBomObs = 1024

Acoef = 0.6548
Bcoef = 0.574

TempMin = -10.0
TempMax = 60.0
TempLastValid = 25.0

I440 = 1
I670 = 2
I870 = 3
I1020 = 4
IWV = 10

MaxNephCalRecords = 32


def Read_Single_Sun_Record(Dbug, filename, Model, DateTime, Data, ValidData):
    with open(filename) as f:
        if Model in [1, 3]:
            Date, Time, Data.Signal.Chan[4], Data.Signal.Chan[3], Data.Signal.Chan[2], Data.Signal.Chan[1], Data.Signal.Chan[7], Data.Signal.Chan[10], Data.Signal.Chan[6], Data.Signal.Chan[5], Data.Temperature = read_data(filename)
            Data.Signal.Chan[8:10] = [0, 0]
        elif Model == 2:
            Date, Time, Data.Signal.Chan[4], Data.Signal.Chan[5], Data.Signal.Chan[2], Data.Signal.Chan[1], Data.Signal.Chan[6], Data.Signal.Chan[3], Data.Signal.Chan[10], Data.Signal.Chan[7], Data.Temperature = read_data(filename)
            Data.Signal.Chan[8:10] = [0, 0]
        elif Model == 4:
            Date, Time, Data.Signal.Chan[4], Data.Signal.Chan[3], Data.Signal.Chan[2], Data.Signal.Chan[1], Data.Signal.Chan[6], Data.Signal.Chan[10], Data.Signal.Chan[7], Data.Signal.Chan[5], Data.Temperature = read_data(filename)
            Data.Signal.Chan[8:10] = [0, 0]
        elif Model in [5, 6]:
            Date, Time, Data.Signal.Chan[4], Data.Signal.Chan[9], Data.Signal.Chan[3], Data.Signal.Chan[2], Data.Signal.Chan[1], Data.Signal.Chan[7], Data.Signal.Chan[8], Data.Signal.Chan[10], Data.Signal.Chan[6], Data.Signal.Chan[5], Data.Temperature = read_data(filename)
        elif Model == 7:
            Date, Time, Data.Signal.Chan[4], Data.Signal.Chan[3], Data.Signal.Chan[2], Data.Signal.Chan[1], Data.Signal.Chan[5], Data.Signal.Chan[10], Data.Signal.Chan[6], Data.Signal.Chan[7], Data.Temperature = read_data(filename, Eof)
            Data.Signal.Chan[8:10] = [0, 0]
        elif Model == 8:
            Date, Time, Data.Signal.Chan[4], Data.Signal.Chan[3], Data.Signal.Chan[2], Data.Signal.Chan[1], Data.Signal.Chan[7], Data.Signal.Chan[5], Data.Signal.Chan[10], Data.Signal.Chan[6], Data.Signal.Chan[8], Data.Temperature = read_data(filename)
            Data.Signal.Chan[9] = 0

    ValidData = 0  # 0 = Valid
    for N in range(1, NumLangleyChannels[Model] + 1):
        if Data.Signal.Chan[N] < MinSignal:
            ValidData = 10
        elif Data.Signal.Chan[N] > MaxSignal:
            ValidData = 11

    if Data.Signal.Chan[IWV] < MinSignal:
        ValidData = 101
    elif Data.Signal.Chan[IWV] > MaxSignal:
        ValidData = 111

    if Eof == 0:
        if len(Date) != 8 or len(Time) != 8:
            print("Read_Sky_Record: Date or time string ne 8")
        DateTime.Day = int(Date[0:2])
        DateTime.Month = int(Date[3:5])
        DateTime.Year = int(Date[6:8])
        DateTime.Hour = int(Time[0:2])
        DateTime.Minute = int(Time[3:5])
        DateTime.Second = int(Time[6:8])
        if DateTime.Year >= 90:
            DateTime.Year += 1900
        else:
            DateTime.Year += 2000
        RawSignal1020 = Data.Signal.Chan[I1020]
        Data.Signal.Chan[I1020] = Data.Signal.Chan[I1020] / TempDep1020(Data.Temperature)
        if Data.Temperature < TempMin or Data.Temperature > TempMax:
            Data.Temperature = -99.99
            ValidData = 12
        if Dbug:
            print("Day        ", DateTime.Day)
            print("Month      ", DateTime.Month)
            print("Year       ", DateTime.Year)
            print("Hour       ", DateTime.Hour)
            print("Minute     ", DateTime.Minute)
            print("Second     ", DateTime.Second)
            print("Signal     ", [Data.Signal.Chan[I] for I in range(1, NumChannels + 1)])
            print("Raw 1020   ", RawSignal1020)
            print("Temperature", Data.Temperature)
            halt()

def Read_Triple_Sun_Record(Dbug, filename, Model, DateTime, Data, ValidData, Eof):
    if Model in [1, 3, 4]:
        Date, Time, *Data.Signal[0].Chan[1:], Data.Temperature = read_data(filename, Eof)
        for I in range(3):
            Data.Signal[I].Chan[8:10] = [0, 0]
    elif Model == 2:
        Date, Time, *Data.Signal[0].Chan[1:], Data.Temperature = read_data(filename, Eof)
        for I in range(3):
            Data.Signal[I].Chan[8:10] = [0, 0]
    elif Model in [5, 6]:
        Date, Time, *Data.Signal[0].Chan[1:], Data.Temperature = read_data(filename, Eof)
    elif Model == 7:
        Date, Time, *Data.Signal[0].Chan[1:], Data.Temperature = read_data(filename, Eof)
        for I in range(3):
            Data.Signal[I].Chan[8:10] = [0, 0]
    elif Model == 8:
        Date, Time, *Data.Signal[0].Chan[1:], Data.Temperature = read_data(filename, Eof)
        for I in range(3):
            Data.Signal[I].Chan[9] = 0

    ValidData = 0  # 0 = Valid
    for N in range(1, NumLangleyChannels[Model] + 1):
        if any(Data.Signal[I].Chan[N] < MinSignal for I in range(3)):
            ValidData = 10
        elif any(Data.Signal[I].Chan[N] > MaxSignal for I in range(3)):
            ValidData = 11

    if any(Data.Signal[I].Chan[IWV] < MinSignal for I in range(3)):
        ValidData = 101
    elif any(Data.Signal[I].Chan[IWV] > MaxSignal for I in range(3)):
        ValidData = 111

    if Eof == 0:
        if len(Date) != 8 or len(Time) != 8:
            print("Read_Sky_Record: Date or time string ne 8")
        DateTime.Day = int(Date[0:2])
        DateTime.Month = int(Date[3:5])
        DateTime.Year = int(Date[6:8])
        DateTime.Hour = int(Time[0:2])
        DateTime.Minute = int(Time[3:5])
        DateTime.Second = int(Time[6:8])
        if DateTime.Year >= 90:
            DateTime.Year += 1900
        else:
            DateTime.Year += 2000
        Data.Signal[0].Chan[I1020], Data.Signal[1].Chan[I1020], Data.Signal[2].Chan[I1020] = [
            Signal.Chan[I1020] / TempDep1020(Data.Temperature) for Signal in Data.Signal[:3]
        ]
        if Data.Temperature < TempMin or Data.Temperature > TempMax:
            Data.Temperature = -99.99
            ValidData = 12
        if Dbug:
            print("Day        ", DateTime.Day)
            print("Month      ", DateTime.Month)
            print("Year       ", DateTime.Year)
            print("Hour       ", DateTime.Hour)
            print("Minute     ", DateTime.Minute)
            print("Second     ", DateTime.Second)
            print("Signal(1)  ", [Data.Signal[0].Chan[I] for I in range(1, NumChannels + 1)])
            print("Signal(2)  ", [Data.Signal[1].Chan[I] for I in range(1, NumChannels + 1)])
            print("Signal(3)  ", [Data.Signal[2].Chan[I] for I in range(1, NumChannels + 1)])
            print("Temperature", Data.Temperature)
            halt()

def Read_Almucantar_Record(Dbug, UnitNumber, Model, Data, Eof):
    if Model in [1, 3, 4]:
        Date, Time, Data.AerosolChannel, Data.Signal, Data.Temperature = read_data(UnitNumber, Eof, NumAerosolChannels, NumAzimuths)
    elif Model == 2:
        Date, Time, Data.AerosolChannel, Data.Signal, Data.Temperature = read_data(UnitNumber, Eof, NumAerosolChannels, NumAzimuths)

    if Eof == 0:
        for I in range(NumAerosolChannels):
            Data[I].DateTime.Day = int(Date[I][0:2])
            Data[I].DateTime.Month = int(Date[I][3:5])
            Data[I].DateTime.Year = int(Date[I][6:8])
            Data[I].DateTime.Hour = int(Time[I][0:2])
            Data[I].DateTime.Minute = int(Time[I][3:5])
            Data[I].DateTime.Second = int(Time[I][6:8])
            if Data[I].DateTime.Year >= 90:
                Data[I].DateTime.Year += 1900
            else:
                Data[I].DateTime.Year += 2000
            if Dbug:
                print("Day        ", Data[I].DateTime.Day)
                print("Month      ", Data[I].DateTime.Month)
                print("Year       ", Data[I].DateTime.Year)
                print("Hour       ", Data[I].DateTime.Hour)
                print("Minute     ", Data[I].DateTime.Minute)
                print("Second     ", Data[I].DateTime.Second)
                print("Temperature", Data[I].Temperature)
                print("Channel    ", Data[I].AerosolChannel)
                halt()
                for J in range(0, NumAzimuths - 3, 6):
                    print([AlmAzimuth[K], Data[I].Signal[K] for K in range(J, min(NumAzimuths, J + 6))])
                halt()

def Read_Almucantar_Split_Record(Dbug, UnitNumber, Model, Data, Eof):
    if Model in [1, 3, 4]:
        Date, Time, Data.AerosolChannel, Data.SignalRight, Data.SignalLeft, Data.Temperature = read_data(UnitNumber, Eof, NumAerosolChannels, NumNAzimuths)
    else:
        print("This model did not do split almucantar scans!")
        halt()

    if Eof == 0:
        for I in range(NumAerosolChannels):
            Data[I].DateTime.Day = int(Date[I][0:2])
            Data[I].DateTime.Month = int(Date[I][3:5])
            Data[I].DateTime.Year = int(Date[I][6:10])
            Data[I].DateTime.Hour = int(Time[I][0:2])
            Data[I].DateTime.Minute = int(Time[I][3:5])
            Data[I].DateTime.Second = int(Time[I][6:8])
            if Data[I].DateTime.Year <= 1910:
                Data[I].DateTime.Year += 100
            if Dbug:
                print("Day        ", Data[I].DateTime.Day)
                print("Month      ", Data[I].DateTime.Month)
                print("Year       ", Data[I].DateTime.Year)
                print("Hour       ", Data[I].DateTime.Hour)
                print("Minute     ", Data[I].DateTime.Minute)
                print("Second     ", Data[I].DateTime.Second)
                print("Temperature", Data[I].Temperature)
                print("Channel    ", Data[I].AerosolChannel)
                halt()
                print([AlmRightAzimuth[J], Data[I].SignalRight[J] for J in range(NumNAzimuths)])
                print([AlmLeftAzimuth[J], Data[I].SignalLeft[J] for J in range(NumNAzimuths)])
                halt()
            
 def Read_PP_Record(Dbug, UnitNumber, Model, Data, Eof):
    if Model in [1, 3, 4]:
        Date, Time, Data.AerosolChannel, Data.Signal, Data.Temperature = read_data(UnitNumber, Eof, NumAerosolChannels, NumPPZeniths)
    elif Model == 2:
        Date, Time, Data.AerosolChannel, Data.Signal, Data.Temperature = read_data(UnitNumber, Eof, NumAerosolChannels, NumPPZeniths)

    if Eof == 0:
        for I in range(NumAerosolChannels):
            Data[I].DateTime.Day = int(Date[I][0:2])
            Data[I].DateTime.Month = int(Date[I][3:5])
            Data[I].DateTime.Year = int(Date[I][6:8])
            Data[I].DateTime.Hour = int(Time[I][0:2])
            Data[I].DateTime.Minute = int(Time[I][3:5])
            Data[I].DateTime.Second = int(Time[I][6:8])
            if Data[I].DateTime.Year >= 90:
                Data[I].DateTime.Year += 1900
            else:
                Data[I].DateTime.Year += 2000
            if Dbug:
                print("Day        ", Data[I].DateTime.Day)
                print("Month      ", Data[I].DateTime.Month)
                print("Year       ", Data[I].DateTime.Year)
                print("Hour       ", Data[I].DateTime.Hour)
                print("Minute     ", Data[I].DateTime.Minute)
                print("Second     ", Data[I].DateTime.Second)
                print("Temperature", Data[I].Temperature)
                print("Channel    ", Data[I].AerosolChannel)
                halt()
                for J in range(0, NumPPZeniths - 3, 6):
                    print([PPZenith[K], Data[I].Signal[K] for K in range(J, min(NumPPZeniths, J + 6))])
                halt()

def Read_PPP_Record(Dbug, UnitNumber, Data, Eof):
    Date, Time, Data.Signal, Data.Temperature = read_data(UnitNumber, Eof)

    if Eof == 0:
        Data.DateTime.Day = int(Date[0:2])
        Data.DateTime.Month = int(Date[3:5])
        Data.DateTime.Year = int(Date[6:8])
        Data.DateTime.Hour = int(Time[0:2])
        Data.DateTime.Minute = int(Time[3:5])
        Data.DateTime.Second = int(Time[6:8])
        if Data.DateTime.Year >= 90:
            Data.DateTime.Year += 1900
        else:
            Data.DateTime.Year += 2000

        if Dbug:
            print("Day        ", Data.DateTime.Day)
            print("Month      ", Data.DateTime.Month)
            print("Year       ", Data.DateTime.Year)
            print("Hour       ", Data.DateTime.Hour)
            print("Minute     ", Data.DateTime.Minute)
            print("Second     ", Data.DateTime.Second)
            print("Temperature", Data.Temperature)
            halt()

def Read_Black_Record(Dbug, Model, UnitNumber, DateTime, Black, Eof):
    Eof = 0
    if Model in [1, 3]:
        Date, Time, Black.Sun.Chan[4], Black.Sun.Chan[3], Black.Sun.Chan[2], Black.Sun.Chan[1], Black.Sun.Chan[7], Black.Sun.Chan[8], Black.Sun.Chan[6], Black.Sun.Chan[5], Black.SkyA[4], Black.SkyA[3], Black.SkyA[2], Black.SkyA[1], Black.SkyK[4], Black.SkyK[3], Black.SkyK[2], Black.SkyK[1] = read_data(UnitNumber, Eof)
        Black.Sun.Chan[9:11] = [0, 0]  # Needed for initialization
    elif Model == 2:
        Date, Time, Black.Sun.Chan[4], Black.Sun.Chan[3], Black.Sun.Chan[2], Black.Sun.Chan[1], Black.Sun.Chan[5], Black.SkyA[4], Black.SkyA[3], Black.SkyA[2], Black.SkyA[1], Black.SkyK[4], Black.SkyK[3], Black.SkyK[2], Black.SkyK[1] = read_data(UnitNumber, Eof)
        Black.Sun.Chan[6:NumChannels+1] = [-1] * (NumChannels - 5)  # Otherwise uninitialized
    elif Model == 4:
        Date, Time, Black.Sun.Chan[4], Black.Sun.Chan[3], Black.Sun.Chan[2], Black.Sun.Chan[1], Black.Sun.Chan[6], Black.Sun.Chan[8], Black.Sun.Chan[7], Black.Sun.Chan[5], Black.SkyA[4], Black.SkyA[3], Black.SkyA[2], Black.SkyA[1], Black.SkyK[4], Black.SkyK[3], Black.SkyK[2], Black.SkyK[1] = read_data(UnitNumber, Eof)
        Black.Sun.Chan[9:11] = [0, 0]  # Needed for initialization

    if Eof == 0:
        if len(Date) != 8 or len(Time) != 8:
            print("Read_Sky_Record: Date or time string ne 8")
        DateTime.Day = int(Date[0:2])
        DateTime.Month = int(Date[3:5])
        DateTime.Year = int(Date[6:8])
        DateTime.Hour = int(Time[0:2])
        DateTime.Minute = int(Time[3:5])
        DateTime.Second = int(Time[6:8])
        if DateTime.Year >= 90:
            DateTime.Year += 1900
        else:
            DateTime.Year += 2000

        if Dbug:
            print("Day        ", DateTime.Day)
            print("Month      ", DateTime.Month)
            print("Year       ", DateTime.Year)
            print("Hour       ", DateTime.Hour)
            print("Minute     ", DateTime.Minute)
            print("Second     ", DateTime.Second)
            print("BlackSun   ", [Black.Sun.Chan[I] for I in range(1, NumChannels + 1)])
            print("BlackSkyA  ", [Black.SkyA[I] for I in range(1, NumSkyA + 1)])
            print("BlackSkyK  ", [Black.SkyK[I] for I in range(1, NumSkyK + 1)])
            halt()

def Read_Sky_Record(Dbug, UnitNumber, DateTime, Data, Eof):
    Date, Time, Data.SkyA[4], Data.SkyA[3], Data.SkyA[2], Data.SkyA[1], Data.SkyK[4], Data.SkyK[3], Data.SkyK[2], Data.SkyK[1], Data.Temperature = read_data(UnitNumber, Eof)

    if Eof == 0:
        if len(Date) != 8 or len(Time) != 8:
            print("Read_Sky_Record: Date or time string ne 8")
        DateTime.Day = int(Date[0:2])
        DateTime.Month = int(Date[3:5])
        DateTime.Year = int(Date[6:8])
        DateTime.Hour = int(Time[0:2])
        DateTime.Minute = int(Time[3:5])
        DateTime.Second = int(Time[6:8])
        if DateTime.Year >= 90:
            DateTime.Year += 1900
        else:
            DateTime.Year += 2000

        if Dbug:
            print("Day        ", DateTime.Day)
            print("Month      ", DateTime.Month)
            print("Year       ", DateTime.Year)
            print("Hour       ", DateTime.Hour)
            print("Minute     ", DateTime.Minute)
            print("Second     ", DateTime.Second)
            print("SkyA (LG)  ", [Data.SkyA[I] for I in range(1, NumSkyA + 1)])
            print("SkyK (HG)  ", [Data.SkyK[I] for I in range(1, NumSkyK + 1)])
            halt()

 def Read_Sphere_Cal(Dbug, UnitNumber, CalDateTime, Data, Eof):
    Data.Instrument, Data.DateTime = read_data(UnitNumber, Eof)
    if Dbug:
        print("Inst #        ", Data.Instrument)
        print("Specified date", CalDateTime)
        print("YMDHMS        ", Data.DateTime)

    while Data.DateTime != CalDateTime:
        skip(UnitNumber, 9)
        Data.Instrument, Data.DateTime = read_data(UnitNumber, Eof)

    skip(UnitNumber, 1)
    for I in range(NumSkyA):
        GainCode, Data.Waveff[I], Data.GainA[I], Data.SdevA[I] = read_data(UnitNumber, Eof)
        if GainCode != 'A':
            print("Error reading sphere cal A data!")
    for I in range(NumSkyK):
        GainCode, Data.Waveff[I], Data.GainK[I], Data.SdevK[I] = read_data(UnitNumber, Eof)
        if GainCode != 'K':
            print("Error reading sphere cal K data!")

    if Eof == 0:
        if Dbug:
            print("Inst #     ", Data.Instrument)
            print("YMDHMS     ", Data.DateTime)
            for I in range(NumSkyA):
                print(Data.Waveff[I], Data.GainA[I], Data.SdevA[I], Data.GainK[I], Data.SdevK[I])
            halt()

def Read_Site_Configuration(Dbug, filename, ObsDate, Config):
    with open(filename,'r') as f:
        next(f)
        next(f)
        Config.Name = f.line
        Config.BoMNumber = int(f.line)
        Config.Id, lat, lon= f.line.split()
        Config.Latitude = float(lat)
        Config.Longitude = float(lon)
        Config.Height = float(f.line.split()[0])
        Config.TimeZone = float(f.line.split()[0])
        next(f)
        if Dbug:
            print("Config.Name          =", Config.Name)
            print("Config.Id            =", Config.Id)
            print("Config.Latitude      =", Config.Latitude)
            print("Config.Longitude     =", Config.Longitude)
            print("Config.Height        =", Config.Height)
            print("Config.TimeZone      =", Config.TimeZone)
            halt()
        Config.Date[0], Config.Date[1], Config.TOffset, Config.CimelNumber, Config.CimelModel, Config.Version, Config.Card, Config.Baro, Config.Wind, Config.Neph, Config.CimCalFile, Config.SkyCalDate = f.line.split()
        if Dbug:
            print("Config.Date(0).Year  =", Config.Date[0].Year)
            print("Config.Date(0).Month =", Config.Date[0].Month)
            print("Config.Date(0).Day   =", Config.Date[0].Day)
            print("Config.Date(1).Year  =", Config.Date[1].Year)
            print("Config.Date(1).Month =", Config.Date[1].Month)
            print("Config.Date(1).Day   =", Config.Date[1].Day)
            print("Config.TOffset.Hour  =", Config.TOffset.Hour)
            print("Config.TOffset.Minute=", Config.TOffset.Minute)
            print("Config.TOffset.Second=", Config.TOffset.Second)
            print("Config.CimelNumber   =", Config.CimelNumber)
            print("Config.CimelModel    =", Config.CimelModel)
            print("Config.Version       =", Config.Version)
            print("Config.Card          =", Config.Card)
            print("Config.Baro          =", Config.Baro)
            print("Config.Wind          =", Config.Wind)
            print("Config.Neph          =", Config.Neph)
            print("Config.CimCalFile    =", Config.CimCalFile)
            print("Config.SkyCalDate    =", Config.SkyCalDate)
            halt()

        if ObsDate - Config.Date[0] < 0:
            if Dbug:
                print("Date(tab)  =", Config.Date[0])
                print("Date(entry)=", ObsDate)
                print("entry-tab=", ObsDate - Config.Date[0])
            print("Error: specified time is earlier than first entry on configuration file")
            halt()
            return

        while ObsDate - Config.Date[0] > 0 and ObsDate - Config.Date[1] > 0:
            if Dbug:
                print("Read_Config, in loop")
                print("Start Date(tab)  =", Config.Date[0])
                print("End Date  (tab)  =", Config.Date[1])
                print("Test Date        =", ObsDate)
                print("test-start=", ObsDate - Config.Date[0])
                print("test-end  =", ObsDate - Config.Date[1])
                halt()

            Config.Date[0], Config.Date[1], Config.TOffset, Config.CimelNumber, Config.CimelModel, Config.Version, Config.Card, Config.Baro, Config.Wind, Config.Neph, Config.CimCalFile, Config.SkyCalDate = read_data(filename)

            if Eof != 0:
                print("Error: specified time is later than last entry on configuration file")
                halt()
                return

    if Dbug:
        print("Config.Id            =", Config.Id)
        print("Config.Latitude      =", Config.Latitude)
        print("Config.Longitude     =", Config.Longitude)
        print("Config.Date(0).Year  =", Config.Date[0].Year)
        print("Config.Date(0).Month =", Config.Date[0].Month)
        print("Config.Date(0).Day   =", Config.Date[0].Day)
        print("Config.Date(1).Year  =", Config.Date[1].Year)
        print("Config.Date(1).Month =", Config.Date[1].Month)
        print("Config.Date(1).Day   =", Config.Date[1].Day)
        print("Config.TOffset.Hour  =", Config.TOffset.Hour)
        print("Config.TOffset.Minute=", Config.TOffset.Minute)
        print("Config.TOffset.Second=", Config.TOffset.Second)
        print("Config.CimelNumber   =", Config.CimelNumber)
        print("Config.CimelModel    =", Config.CimelModel)
        print("Config.Version       =", Config.Version)
        print("Config.Card          =", Config.Card)
        print("Config.Baro          =", Config.Baro)
        print("Config.Wind          =", Config.Wind)
        print("Config.Neph          =", Config.Neph)
        print("Config.CimCalFile    =", Config.CimCalFile)
        print("Config.SkyCalDate    =", Config.SkyCalDate)
        halt()

def Read_Cut_Periods(Dbug, UnitNumber, NumCut, CutList):
    rewind(UnitNumber)
    skip_comments(UnitNumber)
    N = 0
    Line = read_line(UnitNumber, Eof)
    while Eof == 0:
        N += 1
        CutList[N] = read_data(Line)
        if Dbug:
            print("CutList.Date  =", CutList[N].Date)
            print("CutList.AP    =", CutList[N].AP)
            halt()
        Line = read_line(UnitNumber, Eof)
    NumCut = N

def Read_Status(Dbug, UnitNumber, Status, Eof):
    String = read_line(UnitNumber, Eof)
    if Eof == 0:
        if String[0:8] != "# STATUS":
            print("Old or invalid status file")
            print("First line of file=", String[0:8])
            Eof = -1
            return
        Status.ID = read_integer(UnitNumber)
        Date, Time = read_data(UnitNumber)
        Status.DateTime.Year = int(Date[0:2])
        Status.DateTime.Month = int(Date[3:5])
        Status.DateTime.Day = int(Date[6:8])
        Status.DateTime.Hour = int(Time[0:2])
        Status.DateTime.Minute = int(Time[3:5])
        Status.DateTime.Second = int(Time[6:8])
        if Dbug:
            print("Status ID  ", Status.ID)
            print("Day        ", Status.DateTime.Day)
            print("Month      ", Status.DateTime.Month)
            print("Year       ", Status.DateTime.Year)
            print("Hour       ", Status.DateTime.Hour)
            print("Minute     ", Status.DateTime.Minute)
            print("Second     ", Status.DateTime.Second)
            halt()
        String = read_line(UnitNumber)
        while String[0] == '#':
            String = read_line(UnitNumber)
        for I in range(NumLogTimes):
            Status.Data[I] = read_data(String)
            if Dbug:
                print(Status.Data[I])
                halt()
            String = read_line(UnitNumber)

def Read_Cod_File(Dbug, Name, NumAerosolChannels, NumDataPoints, Data):
    SkipLines, NumDataPoints = scan_file(Name, 1)
    print("NumDataPoints=", NumDataPoints)
    print("NumAerosolChannels=", NumAerosolChannels)
    if NumDataPoints > MaxTriplets:
        print("Pending array overflow in Read_Cod_File")
        halt()
    UnitNumber = next_unit_number()
    open_file(UnitNumber, Name, 'read')
    skip(UnitNumber, SkipLines)
    for I in range(NumDataPoints):
        Data[I].DateTime, Data[I].SunZen, Data[I].AirMass, Data[I].Pressure, Data[I].Temperature, Data[I].Aod[0:NumAerosolChannels], Data[I].Wvap, Data[I].Aod1020c, Data[I].Angstrom46, Data[I].Angstrom48, Data[I].Angstrom68, Data[I].Aod_500_Fine, Data[I].Aod_500_Coarse, Data[I].SeparationError, Data[I].CirrusFlag = read_data(UnitNumber)
    close_file(UnitNumber)
    if Dbug:
        print("YMDHMS(1)         ", Data[0].DateTime)
        print("Aod_500_Coarse    ", Data[0].Aod_500_Coarse)
        halt()

def Read_MPL_Record(Dbug, UnitNumber, Data, Eof):
    Data = read_data(UnitNumber, Eof)
    if Dbug:
        print("YMDHMS            ", Data.DateTime)
        print("Depth             ", Data.LayerDepthkm)
        print("Cld/Bkg           ", Data.CloudBackgroundRatio)
        halt()

def Read_Pressure_File(Dbug, Name, NumPressurePoints, NumValidPresPoints, Data, MeanPressure):
    SkipLines, NumPressurePoints = scan_file(Name, 1)
    UnitNumber = next_unit_number()
    open_file(UnitNumber, Name, 'read')
    skip(UnitNumber, SkipLines)
    MeanPressure = 0
    NumValidPresPoints = 0
    for I in range(NumPressurePoints):
        Data[I] = read_data(UnitNumber)
        if Data[I].Pressure > 0:
            MeanPressure += Data[I].Pressure
            NumValidPresPoints += 1
    close_file(UnitNumber)
    if NumValidPresPoints >= MinValidPresPoints:
        MeanPressure /= NumValidPresPoints
    else:
        MeanPressure = -999.9
    if Dbug:
        print("YMDHMS(1)         ", Data[0].DateTime)
        print("Pressure(1)       ", Data[0].Pressure)
        print("No. pts           ", NumPressurePoints)
        print("YMDHMS(N)         ", Data[NumPressurePoints-1].DateTime)
        print("Pressure(N)       ", Data[NumPressurePoints-1].Pressure)
        print("Mean Pressure     ", MeanPressure)
        halt()

def Read_Raw_Environment_Record(Dbug, UnitNumber, Env, Eof):
    Env = read_data(UnitNumber, Eof)
    if Env.DateTime.Year == 99:
        Env.DateTime.Year = 1999
    if Dbug:
        print("YMDHMS           ", Env.DateTime)
        print("Pmean,min,max,sd ", Env.Pressure)
        print("Tmean,min,max,sd ", Env.Temperature)
        print("WVmean,min,max,sd", Env.WindV)
        print("WDmean,min,max,sd", Env.WindD)
        print("VBatt            ", Env.VBatt)
        halt()

def Read_Old_Environment_Record(Dbug, UnitNumber, Env, Eof):
    Env.DateTime.Year, Env.DateTime.Month, Env.DateTime.Day, Env.DateTime.Hour, Env.DateTime.Minute, Env.DateTime.Second, Env.WindV, Env.WindD, Env.Pressure, Env.Temperature, BrickFlow, Env.Vbatt = read_data(UnitNumber, Eof)
    Env.DateTime.Year = Env.DateTime.Year
    if Dbug:
        print("YMDHMS           ", Env.DateTime)
        print("Pmean,min,max,sd ", Env.Pressure)
        print("Tmean,min,max,sd ", Env.Temperature)
        print("WVmean,min,max,sd", Env.WindV)
        print("WDmean,min,max,sd", Env.WindD)
        print("BrickFlow(1)     ", BrickFlow[0])
        print("VBatt            ", Env.VBatt)
        halt()

def Read_Environment_Record(Dbug, UnitNumber, Env, Eof):
    Env = read_data(UnitNumber, Eof)
    if Dbug and Eof == 0:
        print("YMDHMS           ", Env.DateTime)
        print("Pmean,min,max,sd ", Env.Pressure)
        print("Tmean,min,max,sd ", Env.Temperature)
        print("WVmean,min,max,sd", Env.WindV)
        print("WDmean,min,max,sd", Env.WindD)
        print("VBatt            ", Env.VBatt)
        halt()

def Read_Neph_Record(Dbug, UnitNumber, Type, Data, Eof):
    Line = read_line(UnitNumber, Eof)
    if Dbug:
        print("neph read, line=", Line)
        print("neph read, LEN_TRIM(Line)=", len(Line.strip()))
    for I in range(len(Line)):
        if ord(Line[I]) > ord('e'):
            print("Record =", Line[0:4])
            print("Dud char at pos ", I, "=", ord(Line[I]))
            halt()
            Line = Line[:I] + ' ' + Line[I+1:]
    if Type.upper() == 'NEF':
        RecordLength = 57
    elif Type.upper() == 'NNF':
        RecordLength = 51
    else:
        print("Unknown file type in Read_Neph_Record", Type)
        halt()
    while len(Line.strip()) != RecordLength and Eof == 0:
        Line = read_line(UnitNumber, Eof)
    if Eof == 0:
        if Type.upper() == 'NEF':
            Data = read_data(Line)
        elif Type.upper() == 'NNF':
            Data.DateTime.Year, Data.DateTime.Month, Data.DateTime.Day, Data.DateTime.Hour, Data.DateTime.Minute, Data.ScatCoef, Data.CalCoef, Data.Pressure, Data.Temperature, Data.RH = read_data(Line)
            Data.DateTime.Second = 0
            Data.RecordNumber = 0
        else:
            print("Unknown file type in Read_Neph_Record", Type)
            halt()
        while (Data.CalCoef == 0.0 or Data.Pressure <= 400) and Eof == 0:
            print("Zero calibrator or p<400, skipping")
            Line = read_line(UnitNumber, Eof)
            if Type.upper() == 'NEF':
                Data = read_data(Line)
            elif Type.upper() == 'NNF':
                Data.DateTime.Year, Data.DateTime.Month, Data.DateTime.Day, Data.DateTime.Hour, Data.DateTime.Minute, Data.ScatCoef, Data.CalCoef, Data.Pressure, Data.Temperature, Data.RH = read_data(Line)
                Data.DateTime.Second = 0
                Data.RecordNumber = 0
            else:
                print("Unknown file type in Read_Neph_Record", Type)
                halt()
        if Dbug:
            print("RecordNumber ", Data.RecordNumber)
            print("Year         ", Data.DateTime.Year)
            print("Month        ", Data.DateTime.Month)
            print("Day          ", Data.DateTime.Day)
            print("Hour         ", Data.DateTime.Hour)
            print("Minute       ", Data.DateTime.Minute)
            print("Second       ", Data.DateTime.Second)
            print("ScatCoef     ", Data.ScatCoef)
            print("CalCoef      ", Data.CalCoef)
            print("Pressure     ", Data.Pressure)
            print("Temperature  ", Data.Temperature)
            print("RH           ", Data.RH)
            halt()

def GetNephCal(Dbug, UnitNumber, ObsDate, Gain, Offset, Pstd, Tstd):
    rewind(UnitNumber)
    Line = read_line(UnitNumber, Eof)
    while Line[0] == '#':
        Line = read_line(UnitNumber, Eof)
    Pstd, Tstd = read_data(Line)
    skip(UnitNumber, 1)
    Line = read_line(UnitNumber, Eof)
    I = 1
    while Eof == 0:
        Data[I] = read_data(Line)
        I += 1
        Line = read_line(UnitNumber, Eof)
    NumRecords = I - 1
    I = 1
    while ObsDate - Data[I].Date >= 0:
        I += 1
        if I > NumRecords:
            print("GetNephCal, array overrun!")
            print("Check that last date in cal data file is AFTER current date")
            halt()
    DaystoObs = ObsDate - Data[I-1].Date
    DaysBtwCals = Data[I].Date - Data[I-1].Date
    Gain = Data[I-1].Gain + (Data[I].Gain - Data[I-1].Gain) * DaystoObs / DaysBtwCals
    Offset = Data[I-1].Offset + (Data[I].Offset - Data[I-1].Offset) * DaystoObs / DaysBtwCals
    if Dbug:
        print("ObsDate     ", ObsDate)
        print("Date-       ", Data[I-1].Date)
        print("Date+       ", Data[I].Date)
        print("DaystoObs   ", DaystoObs)
        print("DaysBtwCals ", DaysBtwCals)
        print("Gain-       ", Data[I-1].Gain)
        print("Gain+       ", Data[I].Gain)
        print("Gain        ", Gain)
        print("Offset-     ", Data[I-1].Offset)
        print("Offset+     ", Data[I].Offset)
        print("Offset      ", Offset)
        print("Pstd        ", Pstd)
        print("Tstd        ", Tstd)
        halt()

def Read_BoM_Ozone(Dbug, UnitNumber, OzTabUnit, Data, Eof):
    for I in range(1, 5):
        Jstart = (I - 1) * 8 + 1
        Data.Station, Data.LineNumber, Data.Month, Data.Year, Data.StartHourUT[Jstart-1:Jstart+7], Data.EndHourUT[Jstart-1:Jstart+7], Data.WavelengthCode[Jstart-1:Jstart+7], Data.ObsCode[Jstart-1:Jstart+7], Data.DU[Jstart-1:Jstart+7] = read_data(UnitNumber, Eof)
        if Eof != 0:
            return
        if Data.LineNumber != I:
            print("Error in reading ozone file")
            halt()
        if Dbug:
            print("Station         :", Data.Station)
            print("LineNumber      :", Data.LineNumber)
            print("Month           :", Data.Month)
            print("Year            :", Data.Year)
            print("StartHourUT     :", Data.StartHourUT[Jstart-1])
            print("EndHourUT       :", Data.EndHourUT[Jstart-1])
            print("WavelengthCode  :", Data.WavelengthCode[Jstart-1])
            print("ObsCode         :", Data.ObsCode[Jstart-1])
            print("DU              :", Data.DU[Jstart-1])
            halt()
    skip_line(UnitNumber, Eof)  # Skip comment line 5
    Ndays = days_in_month(Data.Month, Data.Year)
    NOk = 0
    for I in range(1, Ndays + 1):
        if Data.DU[I-1] > 10:
            NOk += 1
            RData[NOk-1] = float(Data.DU[I-1])
    Mean, Sdev, MinData, MaxData = stat(RData, NOk)
    if Dbug:
        print("Year, Month NDays NOk:", Data.Year, Data.Month, Ndays, NOk)
        print("Ozone mean, sdev, Min, Max:", Mean, Sdev, MinData, MaxData)
        halt()
    write_data(OzTabUnit, Data.Year, Data.Month, Ndays, NOk, Mean, Sdev, MinData, MaxData, (MaxData - MinData) / 2)

def Read_BoM_Ozone_2011(Dbug, UnitNumber, ObsDate, OzoneDU, Eof):
    StartDate = ObsDate
    StartDate.Day = 1
    rewind(UnitNumber)
    Data.Station, Data.Year, Data.Month, Data.Day, Data.StartHourUT, Data.EndHourUT, Data.WavelengthCode, Data.ObsCode, Data.DU, Data.StdErrDU, Data.Instrument, Data.SerialNumber = read_data(UnitNumber, Eof)
    DateStamp.Year = Data.Year
    DateStamp.Month = Data.Month
    DateStamp.Day = Data.Day
    if ObsDate - DateStamp < 0:
        print("Date preceeds first date in ozone file")
        print("ObsDate  :", ObsDate)
        print("DateStamp:", DateStamp)
        halt()
        return
    while Eof == 0 and DateStamp - StartDate < 0:
        Data.Station, Data.Year, Data.Month, Data.Day, Data.StartHourUT, Data.EndHourUT, Data.WavelengthCode, Data.ObsCode, Data.DU, Data.StdErrDU, Data.Instrument, Data.SerialNumber = read_data(UnitNumber, Eof)
        DateStamp.Year = Data.Year
        DateStamp.Month = Data.Month
        DateStamp.Day = Data.Day
        if Dbug:
            print("Station         :", Data.Station)
            print("Year            :", Data.Year)
            print("Month           :", Data.Month)
            print("Day             :", Data.Day)
            print("StartHourUT     :", Data.StartHourUT)
            print("EndHourUT       :", Data.EndHourUT)
            print("WavelengthCode  :", Data.WavelengthCode)
            print("ObsCode         :", Data.ObsCode)
            print("DU              :", Data.DU)
            print("StdErrDU        :", Data.StdErrDU)
            print("Inst Code       :", Data.Instrument)
            print("Ser. No.        :", Data.SerialNumber)
            halt()
    halt()
    if Eof != 0:
        print("End of file reading ozone file")
        print("EOF encountered before first day of target month")
        print("check data ranges in ozone file and target dates.")
        halt()
    if Dbug:
        print("Station         :", Data.Station)
        print("Year            :", Data.Year)
        print("Month           :", Data.Month)
        print("Day             :", Data.Day)
        print("StartHourUT     :", Data.StartHourUT)
        print("EndHourUT       :", Data.EndHourUT)
        print("WavelengthCode  :", Data.WavelengthCode)
        print("ObsCode         :", Data.ObsCode)
        print("DU              :", Data.DU)
        print("StdErrDU        :", Data.StdErrDU)
        print("Inst Code       :", Data.Instrument)
        print("Ser. No.        :", Data.SerialNumber)
        halt()
    EndDate = StartDate
    DayMonth = days_in_month(StartDate.Month, StartDate.Year)
    EndDate.Day = DayMonth
    OzoneDU[0] = 0
    OzoneDU[1:32] = -1
    N = 0
    while Eof == 0 and DateStamp - EndDate <= 0:
        N += 1
        OzoneDU[DateStamp.Day] = Data.DU
        OzoneDU[0] += Data.DU
        Data.Station, Data.Year, Data.Month, Data.Day, Data.StartHourUT, Data.EndHourUT, Data.WavelengthCode, Data.ObsCode, Data.DU, Data.StdErrDU, Data.Instrument, Data.SerialNumber = read_data(UnitNumber, Eof)
        DateStamp.Year = Data.Year
        DateStamp.Month = Data.Month
        DateStamp.Day = Data.Day
    if Eof != 0:
        print("End of file reading last month of ozone file")
        print("However, for a short final month, the ozone values")
        print("of the days present will be available")
        print("and the monthly mean will be based on the")
        print("average of the available days.")
        halt()
    OzoneDU[0] /= N
    if Dbug:
        print("Start date      :", StartDate)
        print("EndDate         :", EndDate)
        print("NumPts          :", N)
        print("MeanDU          :", OzoneDU[0])
        print("Daily           :", OzoneDU[1:DayMonth+1])
        halt()

def Read_BoM_Transmission(Dbug, UnitNumber, NumBomChannels, Wavelength_BoM, NumObs, Data, Eof):
    skip(UnitNumber, 5)  # Skip 5 header records
    NumBomChannels = read_integer(UnitNumber)
    for I in range(NumBomChannels):
        Wavelength_BoM[I] = read_float(UnitNumber)
    skip(UnitNumber, 1)  # Skip 1 header record

    N = 0
    N += 1
    Data[N].Date, Data[N].DayOfYear, Data[N].Minute, Data[N].Day, Data[N].Pressure, Data[N].SolarZD, Data[N].Airmass, Data[N].FLT, Data[N].StDev, Data[N].Tran[:NumBomChannels], Data[N].Diffuse[:NumBomChannels] = read_data(UnitNumber, Eof)

    if Dbug and Eof == 0:
        print("Date      :", Data[N].Date)
        print("DayOfYear :", Data[N].DayOfYear)
        print("Minute    :", Data[N].Minute)
        print("Hour      :", Data[N].Day)
        print("Pressure  :", Data[N].Pressure)
        print("SolarZD   :", Data[N].SolarZD)
        print("Airmass   :", Data[N].Airmass)
        print("Wavelength:", [Wavelength_BoM[I] for I in range(NumBomChannels)])
        print("Tran(1)   :", [Data[N].Tran[I] for I in range(NumBomChannels)])
        print("Diffuse(1):", [Data[N].Diffuse[I] for I in range(NumBomChannels)])
        print("Eof       :", Eof)
        halt()

    while Eof == 0:
        N += 1
        Data[N].Date, Data[N].DayOfYear, Data[N].Minute, Data[N].Day, Data[N].Pressure, Data[N].SolarZD, Data[N].Airmass, Data[N].FLT, Data[N].StDev, Data[N].Tran[:NumBomChannels], Data[N].Diffuse[:NumBomChannels] = read_data(UnitNumber, Eof)

        if Dbug and Eof == 0:
            print("Date      :", Data[N].Date)
            print("DayOfYear :", Data[N].DayOfYear)
            print("Minute    :", Data[N].Minute)
            print("Hour      :", Data[N].Day)
            print("Pressure  :", Data[N].Pressure)
            print("SolarZD   :", Data[N].SolarZD)
            print("Airmass   :", Data[N].Airmass)
            print("Wavelength:", [Wavelength_BoM[I] for I in range(NumBomChannels)])
            print("Tran(1)   :", [Data[N].Tran[I] for I in range(NumBomChannels)])
            print("Diffuse(1):", [Data[N].Diffuse[I] for I in range(NumBomChannels)])
            print("Eof       :", Eof)
            halt()

    NumObs = N - 1

def Read_BoM_Pressure(Dbug, UnitNumber, ObsDate, PData, NumPressurePoints, NumValidPresPoints, MeanPressure):
    rewind(UnitNumber)
    skip(UnitNumber, 1)  # Skip field descriptor record
    TimeZone.Hour, TimeZone.Minute = read_data(UnitNumber, Eof)

    if Eof == 0:
        TimeZone.Year = 0
        TimeZone.Month = 0
        TimeZone.Day = 0
        TimeZone.Hour = -TimeZone.Hour
        TimeZone.Minute = -TimeZone.Minute
        TimeZone.Second = 0
    else:
        print("Error in Read_BoM_Pressure")
        print("Eof while reading time zone from file")
        halt()

    Data.Id, Data.Station, Data.DateTime.Year, Data.DateTime.Month, Data.DateTime.Day, Data.DateTime.Hour, Data.DateTime.Minute, Data.MSLP, Data.Qual_MSLP, Data.Station_Level_Pressure, Data.Qual_SLP, Data.QNHP, Data.Qual_QNHP, Data.AWS_Flag, Data.EOR_Hash = read_data(UnitNumber, Eof)

    Data.DateTime.Second = 0

    DateStamp.Year = Data.DateTime.Year
    DateStamp.Month = Data.DateTime.Month
    DateStamp.Day = Data.DateTime.Day

    DateTimeStamp.Year = Data.DateTime.Year
    DateTimeStamp.Month = Data.DateTime.Month
    DateTimeStamp.Day = Data.DateTime.Day
    DateTimeStamp.Hour = Data.DateTime.Hour
    DateTimeStamp.Minute = Data.DateTime.Minute
    DateTimeStamp.Second = 0

    if ObsDate - DateStamp < 0:
        print("Date preceeds first date in pressure file")
        print("ObsDate  :", ObsDate)
        print("DateStamp:", DateStamp)
        halt()
        return

    while Eof == 0 and DateStamp - ObsDate < 0:
        Data.Id, Data.Station, Data.DateTime.Year, Data.DateTime.Month, Data.DateTime.Day, Data.DateTime.Hour, Data.DateTime.Minute, Data.MSLP, Data.Qual_MSLP, Data.Station_Level_Pressure, Data.Qual_SLP, Data.QNHP, Data.Qual_QNHP, Data.AWS_Flag, Data.EOR_Hash = read_data(UnitNumber, Eof)

        Data.DateTime.Second = 0

        DateStamp.Year = Data.DateTime.Year
        DateStamp.Month = Data.DateTime.Month
        DateStamp.Day = Data.DateTime.Day

        DateTimeStamp.Year = Data.DateTime.Year
        DateTimeStamp.Month = Data.DateTime.Month
        DateTimeStamp.Day = Data.DateTime.Day
        DateTimeStamp.Hour = Data.DateTime.Hour
        DateTimeStamp.Minute = Data.DateTime.Minute
        DateTimeStamp.Second = 0

        if Dbug:
            print("Station         :", Data.Station)
            print("Year            :", Data.DateTime.Year)
            print("Month           :", Data.DateTime.Month)
            print("Day             :", Data.DateTime.Day)
            print("Hour            :", Data.DateTime.Hour)
            print("Minute          :", Data.DateTime.Minute)
            print("Pressure        :", Data.Station_Level_Pressure)
            print("Quality flag    :", Data.Qual_SLP)
            halt()

    if Eof != 0:
        print("End of file reading BoM pressure file")
        halt()
        return

    if Dbug:
        print("Station         :", Data.Station)
        print("Year            :", Data.DateTime.Year)
        print("Month           :", Data.DateTime.Month)
        print("Day             :", Data.DateTime.Day)
        print("Hour            :", Data.DateTime.Hour)
        print("Minute          :", Data.DateTime.Minute)
        print("Pressure        :", Data.Station_Level_Pressure)
        print("Quality flag    :", Data.Qual_SLP)
        halt()

    MeanPressure = 0
    N = 0

    while Eof == 0 and DateStamp == ObsDate:
        N += 1
        PData[N].DateTime = Data.DateTime + TimeZone
        PData[N].Pressure = Data.Station_Level_Pressure
        MeanPressure += PData[N].Pressure

        Data.Id, Data.Station, Data.DateTime.Year, Data.DateTime.Month, Data.DateTime.Day, Data.DateTime.Hour, Data.DateTime.Minute, Data.MSLP, Data.Qual_MSLP, Data.Station_Level_Pressure, Data.Qual_SLP, Data.QNHP, Data.Qual_QNHP, Data.AWS_Flag, Data.EOR_Hash = read_data(UnitNumber, Eof)

        Data.DateTime.Second = 0

        DateStamp.Year = Data.DateTime.Year
        DateStamp.Month = Data.DateTime.Month
        DateStamp.Day = Data.DateTime.Day

        DateTimeStamp.Year = Data.DateTime.Year
        DateTimeStamp.Month = Data.DateTime.Month
        DateTimeStamp.Day = Data.DateTime.Day
        DateTimeStamp.Hour = Data.DateTime.Hour
        DateTimeStamp.Minute = Data.DateTime.Minute
        DateTimeStamp.Second = 0

    if Eof != 0:
        print("End of file reading BoM pressure file")
        halt()
        return

    MeanPressure /= N
    NumPressurePoints = N
    NumValidPresPoints = N

    if Dbug:
        print("Date            :", ObsDate)
        print("NumPts          :", N)
        print("MeanPressure    :", MeanPressure)
        halt()

def Get_Time_Correction(Dbug, UnitNumber, ObsDate, TimeCorrI, Eof):
    rewind(UnitNumber)
    skip(UnitNumber, 1)
    DateStamp.Year, DateStamp.Month, DateStamp.Day, TimeCorrR = read_data(UnitNumber, Eof)

    if ObsDate - DateStamp < 0:
        print("Date preceeds first date in time correction file")
        print("ObsDate  :", ObsDate)
        print("DateStamp:", DateStamp)
        halt()
        return

    while Eof == 0 and DateStamp - ObsDate < 0:
        DateStamp.Year, DateStamp.Month, DateStamp.Day, TimeCorrR = read_data(UnitNumber, Eof)

    if DateStamp != ObsDate:
        print("Datestamp should match Observed Date, but doesnt")
        print("ObsDate  :", ObsDate)
        print("DateStamp:", DateStamp)
        halt()

    TimeCorrI = round(TimeCorrR)
    print("Date, Time correction:", DateStamp, TimeCorrR, TimeCorrI)

    if Eof != 0:
        print("End of file reading time correction file")
        print("NB time correction file needs at least one valid record")
        print("after last day of last month processed.")
        halt()
        return

def Rayleigh(Dbug, Wavelength, SurfacePressure, RayleighOD):
    DensitySTP = 2.54743e19  # cm-3
    WMicron = Wavelength / 1000.0
    WInvSq = 1.0 / WMicron**2
    if Wavelength > 230.0:
        RefIndex = 1.0 + 1.0e-8 * (5791817.0 / (238.0185 - WInvSq) + 167909.0 / (57.362 - WInvSq))
    else:
        RefIndex = 1.0 + 1.0e-8 * (8060.51 + 2480990.0 / (132.274 - WInvSq) + 17455.7 / (39.32957 - WInvSq))

    NumWavTab = 36
    WavTab = [200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 240.0, 250.0, 260.0,
              270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0,
              370.0, 380.0, 390.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0,
              750.0, 800.0, 850.0, 900.0, 950.0, 1000.0]
    DepolTab = [4.545e-2, 4.384e-2, 4.221e-2, 4.113e-2, 4.004e-2, 3.895e-2, 3.785e-2, 3.675e-2,
                3.565e-2, 3.455e-2, 3.400e-2, 3.289e-2, 3.233e-2, 3.178e-2, 3.178e-2, 3.122e-2,
                3.066e-2, 3.066e-2, 3.010e-2, 3.010e-2, 3.010e-2, 2.955e-2, 2.955e-2, 2.955e-2,
                2.899e-2, 2.842e-2, 2.842e-2, 2.786e-2, 2.786e-2, 2.786e-2, 2.786e-2, 2.730e-2,
                2.730e-2, 2.730e-2, 2.730e-2, 2.730e-2]

    ColumnDensity = AvogadroConstant * SurfacePressure / (MolecularWeightDryAir * Gravity * 100.0)
    Depolarization = xlin(NumWavTab, WavTab, DepolTab, Wavelength)
    KingFactor = (6.0 + 3.0 * Depolarization) / (6.0 - 7.0 * Depolarization)
    Wcm = Wavelength * 1.0e-7
    RayleighCsa = 24.0 * math.pi**3 * (RefIndex**2 - 1.0)**2 / (Wcm**4 * DensitySTP**2 * (RefIndex**2 + 2.0)**2) * KingFactor
    RayleighOD = ColumnDensity * RayleighCsa

    if Dbug:
        print("Surface Pressure :", SurfacePressure)
        print("Column density   :", ColumnDensity)
        print("Wavelength       :", Wavelength)
        print("Ref Index        :", RefIndex)
        print("Depolarization   :", Depolarization)
        print("KingFactor       :", KingFactor)
        print("RayleighCSA      :", RayleighCsa)
        print("RayleighOD       :", RayleighOD)
        halt()

def CheckTripletCv(Dbug, NumPoints, Nstart, Airmass, CoVar, SpreadFlag, NumOK, Index):
    MinAirmass = 2
    MaxAirmass = 6
    MinPtsUnitAirmass = 2
    MaxCoVar = 1.0

    NumOK = 0
    NumPtsUnitAirmass = [0] * (MaxAirmass + 1)

    for I in range(Nstart, Nstart + NumPoints):
        if MinAirmass <= Airmass[I] <= MaxAirmass and CoVar[I] <= MaxCoVar:
            NumOK += 1
            Index[NumOK - 1] = I
            J = int(Airmass[I])
            if J > MaxAirmass:
                print("Array overflow, airmass too big!")
                print("INT(Airmass)   :", J)
                print("Max Airmass    :", MaxAirmass)
                halt()
            NumPtsUnitAirmass[J] += 1

    SpreadFlag = True
    for J in range(MinAirmass, MaxAirmass):
        if NumPtsUnitAirmass[J] < MinPtsUnitAirmass:
            SpreadFlag = False
            print("Insufficient points in airmass interval", J)

    if Dbug:
        print("NumOK    ", NumOK)
        for J in range(1, 7):
            print("Points in Airmassrange", J, NumPtsUnitAirmass[J])
        halt()

def CheckFitQuality(Dbug, NumPoints, X, Y, Nref, MaxSdevFit, SpreadFlag, FitFlag, Index=None):
    MinAirmass = 2
    MaxAirmass = 6
    MinPtsUnitAirmass = 4

    N = Nref
    Intercept, Slope, Residual, Erms, DelIntercept, DelSlope = boxfit(NumPoints, X[:, N], Y[:, N])

    if Dbug:
        print("Num points in fit=", NumPoints)
        print("Intercept         ", Intercept)
        print("Slope             ", Slope)
        print("Residual          ", Residual)
        print("Erms              ", Erms)
        print("DelIntercept      ", DelIntercept)
        print("DelSlope          ", DelSlope)
        halt()

    NumPtsUnitAirmass = [0] * (MaxAirmass + 1)
    I = 0

    for J in range(NumPoints):
        Error = Y[J, Nref] - (Intercept + Slope * X[J, Nref])
        if abs(Error) < 1.5 * Erms:
            I += 1
            X[I - 1, :] = X[J, :]
            Y[I - 1, :] = Y[J, :]
            if Index is not None:
                Index[I - 1] = J
            K = int(X[I - 1, Nref])
            NumPtsUnitAirmass[K] += 1
        else:
            print("Rejecting point, am=:", X[J, Nref], "LnV=", Y[J, Nref], "Error=", Error, "1.5sigma=", 1.5 * Erms)

    print("Points satisfying Pass 1 test=", NumPoints)
    NumPoints = I
    print("Points satisfying Pass 2 test=", NumPoints)

    SpreadFlag = True
    for J in range(MinAirmass, MaxAirmass):
        if NumPtsUnitAirmass[J] < MinPtsUnitAirmass:
            SpreadFlag = False
            print("Insufficient points in airmass interval", J)

    if Dbug:
        for J in range(1, 7):
            print("Points in Airmassrange", J, NumPtsUnitAirmass[J])
        halt()

    if SpreadFlag:
        Intercept, Slope, Residual, Erms, DelIntercept, DelSlope = boxfit(NumPoints, X[:NumPoints, N], Y[:NumPoints, N])
        if Dbug:
            print("Num points in fit=", NumPoints)
            print("Intercept         ", Intercept)
            print("Slope             ", Slope)
            print("Residual          ", Residual)
            print("Erms              ", Erms)
            print("DelIntercept      ", DelIntercept)
            print("DelSlope          ", DelSlope)
            halt()

        FitUnit = next_unit_number()
        with open("Langley2.dat", "w") as file:
            file.write(f"#Num points in fit= {NumPoints}\n")
            file.write(f"#Intercept          {Intercept}\n")
            file.write(f"#Slope              {Slope}\n")
            file.write(f"#Residual           {Residual}\n")
            file.write(f"#Erms               {Erms}\n")
            file.write(f"#DelIntercept       {DelIntercept}\n")
            file.write(f"#DelSlope           {DelSlope}\n")
            for I in range(NumPoints):
                file.write(f"{X[I, Nref]} {Y[I, Nref]}\n")

        print("Erms              ", Erms)
        print("MaxSdevFit        ", MaxSdevFit)
        if Erms < MaxSdevFit:
            FitFlag = True
        else:
            FitFlag = False
            print("FitFlag set false because standard deviation=", Erms)
    else:
        FitFlag = False

def Langley_Mean_Temperature(Dbug, NumPoints, Index, DetectorTemp, MeanDetectorTemp):
    ValidTemp = True
    MeanDetectorTemp = 0

    for I in range(NumPoints):
        MeanDetectorTemp += DetectorTemp[Index[I] - 1]
        if DetectorTemp[Index[I] - 1] < TempMin or DetectorTemp[Index[I] - 1] > TempMax:
            ValidTemp = False
            break

    if ValidTemp:
        MeanDetectorTemp /= NumPoints
    else:
        MeanDetectorTemp = -99.99

    if Dbug:
        print("NumPoints      ", NumPoints)
        print("MeanDetectorTemp", MeanDetectorTemp)
        halt()

def General_Fit(Dbug, NumPoints, Index, AirMassODAerosolRef, LnV, AirMassODRayleigh, AirMassODOzone, Intercept, DelIntercept, Slope, DelSlope):
    Weight = [1.0] * NumPoints

    X = [AirMassODAerosolRef[Index[I] - 1] for I in range(NumPoints)]
    Y = [LnV[Index[I] - 1] + AirMassODRayleigh[Index[I] - 1] + AirMassODOzone[Index[I] - 1] for I in range(NumPoints)]

    Intercept, Slope, Residual, Erms, DelIntercept, DelSlope = elfit(NumPoints, Weight, X, Y)

    if Dbug:
        print("NumPoints      ", NumPoints)
        print("X              ", X)
        print("Y              ", Y)
        print("Intercept      ", Intercept)
        print("DelIntercept   ", DelIntercept)
        print("Slope          ", Slope)
        print("DelSlope       ", DelSlope)
        print("Rms error      ", Erms)
        print("Residual       ", Residual)
        halt()

def CheckFileDosUnix(UnitNumber):
    DOSFlag = False
    Line = read_line(UnitNumber)
    L = len(Line.strip())
    CR = 13

    if ord(Line[L - 1]) == CR:
        DOSFlag = True
    else:
        DOSFlag = False

    if DOSFlag:
        print("File is in DOS format. Continue (y/n)?")
        YesNo = input()
        if YesNo.upper() == 'N':
            print("Aborting.")
            exit()

    backspace(UnitNumber)

def Combine_Stat(Ndat, Size, MeanIn, SdevIn, MinIn, MaxIn, Mean, Sdev, MinOut, MaxOut):
    if Ndat > 0:
        Sum1 = 0
        Sum2 = 0
        for N in range(Ndat):
            Sum1 += Size[N] * MeanIn[N]
            Sum2 += Size[N]

        Mean = Sum1 / Sum2

        Sum1 = 0
        Sum2 = 0
        for N in range(Ndat):
            Sum1 += Size[N] * (SdevIn[N] ** 2 + MeanIn[N] ** 2)
            Sum2 += Size[N]

        Sdev = math.sqrt(Sum1 / Sum2 - Mean ** 2)

        MinOut = MinIn[0]
        MaxOut = MaxIn[0]
        for N in range(1, Ndat):
            MinOut = min(MinOut, MinIn[N])
            MaxOut = max(MaxOut, MaxIn[N])
    else:
        print("Combine_Stat: number of data points <= 0")
        Mean = 0
        Sdev = 0

def Elfit(N, W, X, Y, A, S, R, E, EA, ES):
    SW = 0
    SWX = 0
    SWY = 0
    SWXX = 0
    SWXY = 0

    for I in range(N):
        SW += W[I]
        SWX += W[I] * X[I]
        SWY += W[I] * Y[I]
        SWXX += W[I] * X[I] * X[I]
        SWXY += W[I] * X[I] * Y[I]

    DET = SW * SWXX - SWX * SWX
    A = (SWXX * SWY - SWX * SWXY) / DET
    S = (SWXY * SW - SWX * SWY) / DET

    R = 0
    for I in range(N):
        R += W[I] * (A + S * X[I] - Y[I]) ** 2

    E = math.sqrt(R / SW)
    EA = math.sqrt((SWXX / DET) * (R / (N - 2)))
    ES = math.sqrt((SW / DET) * (R / (N - 2)))

def Boxfit(N, X, Y, A, S, R, E, EA, ES):
    S0 = 0
    SX1 = 0
    SX2 = 0
    SX1Y = 0
    SX2Y = 0

    for I in range(N):
        S0 += 1
        SX1 += 1 / X[I]
        SX2 += 1 / (X[I] * X[I])
        SX1Y += Y[I] / X[I]
        SX2Y += Y[I] / (X[I] * X[I])

    DET = S0 * SX2 - SX1 * SX1
    A = (S0 * SX2Y - SX1 * SX1Y) / DET
    S = (SX2 * SX1Y - SX1 * SX2Y) / DET

    R = 0
    T = 0
    for I in range(N):
        R += (A + S * X[I] - Y[I]) ** 2
        T += ((A + S * X[I] - Y[I]) / X[I]) ** 2

    E = math.sqrt(R / (N - 2))
    EA = math.sqrt((S0 / DET) * (T / (N - 2)))
    ES = math.sqrt((SX2 / DET) * (T / (N - 2)))

def Xlin(N, X, Y, Z):
    if X[0] < X[N - 1]:
        Xmin = X[0]
        Xmax = X[N - 1]
        Imin = 0
        Imax = N - 1
    else:
        Xmin = X[N - 1]
        Xmax = X[0]
        Imin = N - 1
        Imax = 0

    if Z <= Xmin:
        Xlin = Y[Imin]
    elif Z >= Xmax:
        Xlin = Y[Imax]
    else:
        while abs(Imax - Imin) > 1:
            I = (Imin + Imax) // 2
            if Z >= X[I]:
                Imin = I
            else:
                Imax = I

        Xlin = Y[Imin] + (Y[Imax] - Y[Imin]) * (Z - X[Imin]) / (X[Imax] - X[Imin])

    return Xlin

def GetPressure(ObsDateTime, PData, PressureOption, NumPressurePoints, NumValidPresPoints, DailyMeanPressure, DefaultSurfacePressure):
    if PressureOption == 0 or NumPressurePoints == 0 or NumValidPresPoints < MinValidPresPoints:
        GetPressure = DefaultSurfacePressure
    elif PressureOption == 1 or PressureOption == 3:
        I = 0
        TimeDiff1 = ObsDateTime - PData[I].DateTime
        I = NumPressurePoints - 1
        TimeDiff2 = ObsDateTime - PData[I].DateTime

        if TimeDiff1 < 0 or TimeDiff2 > 0:
            GetPressure = DailyMeanPressure
            print("Observation outside range of pressure tabulation: assuming daily mean")
            print("First time on file=", PData[0].DateTime)
            print("Last time  on file=", PData[NumPressurePoints - 1].DateTime)
            print("ObsTime           =", ObsDateTime)
            halt()
        else:
            IPresMatch = 0
            for I in range(NumPressurePoints):
                TimeDiff = PData[I].DateTime - ObsDateTime
                if TimeDiff > 0:
                    IPresMatch = I
                    break

            if IPresMatch == 0:
                print("Failed to match pressure time to obs time")
                print("ObsDateTime   =", ObsDateTime)
                print("PDateTime(1)  =", PData[0].DateTime)
                print("PDateTime(N)  =", PData[NumPressurePoints - 1].DateTime)
                halt()

            GetPressure = PData[IPresMatch].Pressure

            if GetPressure < 0:
                GetPressure = DailyMeanPressure

    return GetPressure

def GetAirMass(Index, ApparentZenith):
    ScaleHeight = [0.0, 6.58, 1.0, 22.0]
    EarthRadius = 6370.0

    X = math.cos(ApparentZenith * math.pi / 180)
    R = ScaleHeight[Index] / EarthRadius
    D = math.sqrt(X * X + 2 * R + R * R)

    GetAirMass = (1 + R) / D

    return GetAirMass

def VoltsAt1AU(V, R):
    VoltsAt1AU = V * R * R
    print("Measured voltage      :", V)
    print("Sun-Earth distance, au:", R)
    print("Voltage at 1au        :", VoltsAt1AU)
    halt()
    return VoltsAt1AU

def AGSNet2Aeronet(N):
    if N in [1, 2, 3, 4]:
        AGSNet2Aeronet = N + 120
    elif N in [5, 6]:
        AGSNet2Aeronet = N + 250
    else:
        AGSNet2Aeronet = N
    return AGSNet2Aeronet

def TempDep1020(TempC):
    GainCoefficient = 0.0025

    if TempMin < TempC < TempMax:
        TempDep1020 = 1.0 + GainCoefficient * (TempC - 25.0)
        TempLastValid = TempC
    else:
        TempDep1020 = 1.0

    return TempDep1020


