import math
import pandas as pd
import numpy as np
import time_module as tm
import datetime as dt
from constants_module import *

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
NumSkyA = 4
NumSkyK =4
I440 = 1
I670 = 2
I870 = 3
I1020 = 4
IWV = 10
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

def Read_Site_Configuration(Dbug, filename, ObsDate):
    
    Config = pd.read_csv(filename, skiprows=8, header=None, delimiter=r'\s+', usecols=range(0,14))
    Config.rename(columns={0:'StartDate', 1:'EndDate', 2:'TOffSetHour', 3:'TOffSetMinute', 4:'TOffSetSecond',
                           5:'CimelNumber', 6:'CimelModel', 7:'Version', 8:'Card', 9:'Baro', 10:'Wind', 11:'Neph',
                           12:'CimCalFile',13:'SkyCalDate'}, inplace=True)
    
    Config['TOffset'] = [dt.timedelta(hours=h,minutes=m,seconds=s) for h,m,s in zip(Config.TOffSetHour,Config.TOffSetMinute,Config.TOffSetSecond)]
    Config['StartDate'] = Config['StartDate'].apply(dateint2dateclass)
    Config['EndDate'] = Config['EndDate'].apply(dateint2dateclass)
    #Config['SkyCalDate'] = Config['SkyCalDate'].apply(dateint2dateclass)  # input file contains 'dummy'
   
    with open(filename,'r') as f:
        next(f)
        next(f)
        Name = f.readline().rstrip()
        BoMNumber = int(f.readline().rstrip())
        Id, lat, lon= f.readline().rstrip().split()
        Height = float(f.readline().rstrip().split()[0])
        TimeZone = float(f.readline().rstrip().split()[0])
    Config.attrs = {'Name':Name, 'BoMNumber':BoMNumber, 'Id':Id, 'Latitude':float(lat), 'Longitude':float(lon), 'Height':Height, 'TimeZone':TimeZone}
    pre_obs = [i for i,x in enumerate(Config['StartDate']) if (x<ObsDate)]  # 
    #print(pre_obs)
    ConfigOut = Config.iloc[pre_obs[-1]]
#    if Dbug:
#        print("Config.Name          =", Config.Name)
#        print("Config.Id            =", Config.Id)
#        print("Config.Latitude      =", Config.Latitude)
#        print("Config.Longitude     =", Config.Longitude)
#        print("Config.Height        =", Config.Height)
#        print("Config.TimeZone      =", Config.TimeZone)
#        #halt()
#
#    if Dbug:
#        print("Config.Date(0).Year  =", Config.StartDate[0].year)
#        print("Config.Date(0).Month =", Config.StartDate[0].month)
#        print("Config.Date(0).Day   =", Config.StartDate[0].day)
#        print("Config.Date(1).Year  =", Config.EndDate[0].year)
#        print("Config.Date(1).Month =", Config.EndDate[0].month)
#        print("Config.Date(1).Day   =", Config.EndDate[0].day)
##            print("Config.TOffset.Hour  =", Config.TOffset.Hour)
##            print("Config.TOffset.Minute=", Config.TOffset.Minute)
##            print("Config.TOffset.Second=", Config.TOffset.Second)
#        print("Config.CimelNumber   =", Config.CimelNumber[0])
#        print("Config.CimelModel    =", Config.CimelModel[0])
#        print("Config.Version       =", Config.Version[0])
#        print("Config.Card          =", Config.Card[0])
#        print("Config.Baro          =", Config.Baro[0])
#        print("Config.Wind          =", Config.Wind[0])
#        print("Config.Neph          =", Config.Neph[0])
#        print("Config.CimCalFile    =", Config.CimCalFile[0])
#        print("Config.SkyCalDate    =", Config.SkyCalDate[0])
#            #halt()
#
#        if ObsDate - Config.Date[0] < 0:
#            if Dbug:
#                print("Date(tab)  =", Config.Date[0])
#                print("Date(entry)=", ObsDate)
#                print("entry-tab=", ObsDate - Config.Date[0])
#            print("Error: specified time is earlier than first entry on configuration file")
#            #halt()
#            return
#
#        while ObsDate - Config.Date[0] > 0 and ObsDate - Config.Date[1] > 0:
#            if Dbug:
#                print("Read_Config, in loop")
#                print("Start Date(tab)  =", Config.Date[0])
#                print("End Date  (tab)  =", Config.Date[1])
#                print("Test Date        =", ObsDate)
#                print("test-start=", ObsDate - Config.Date[0])
#                print("test-end  =", ObsDate - Config.Date[1])
#                #halt()
#
#
#            if Eof != 0:
#                print("Error: specified time is later than last entry on configuration file")
#                #halt())
#                return
#
#    if Dbug:
#        print("Config.Id            =", Config.Id)
#        print("Config.Latitude      =", Config.Latitude)
#        print("Config.Longitude     =", Config.Longitude)
#        print("Config.Date(0).Year  =", Config.Date[0].Year)
#        print("Config.Date(0).Month =", Config.Date[0].Month)
#        print("Config.Date(0).Day   =", Config.Date[0].Day)
#        print("Config.Date(1).Year  =", Config.Date[1].Year)
#        print("Config.Date(1).Month =", Config.Date[1].Month)
#        print("Config.Date(1).Day   =", Config.Date[1].Day)
#        print("Config.TOffset.Hour  =", Config.TOffset.Hour)
#        print("Config.TOffset.Minute=", Config.TOffset.Minute)
#        print("Config.TOffset.Second=", Config.TOffset.Second)
#        print("Config.CimelNumber   =", Config.CimelNumber)
#        print("Config.CimelModel    =", Config.CimelModel)
#        print("Config.Version       =", Config.Version)
#        print("Config.Card          =", Config.Card)
#        print("Config.Baro          =", Config.Baro)
#        print("Config.Wind          =", Config.Wind)
#        print("Config.Neph          =", Config.Neph)
#        print("Config.CimCalFile    =", Config.CimCalFile)
#        print("Config.SkyCalDate    =", Config.SkyCalDate)
#        #halt()


    return ConfigOut

def dateint2dateclass(d):
    d=str(d)
    d=dt.date(int(d[0:4]),int(d[4:6]),int(d[6:8]))
    return d

def read_triple_sun_record(filename, Model):
    if Model in [1,3]:
        order = [4,3,2,1,7,10,6,5]
    elif Model==2:
        order = [4,5,2,1,6,3,10,7]
    elif Model==4:
        order = [4,3,2,1,6,10,7,5]
    elif Model in [5,6]:
        order = [4,9,3,2,1,7,8,10,6,5]
    elif Model==7:
        order = [4,3,2,1,5,10,6,7]
    elif Model==8:
        order = [4,3,2,1,5,7,10,6,8]
  
    s = pd.read_csv(filename, header=None)
    s.rename(columns={32:'Temperature'}, inplace=True)
    s['Date'] = s[0].apply(lambda x: dt.datetime.strptime(x,'%d:%m:%y'))
    s['Time'] = s[1].apply(lambda x: dt.datetime.strptime(x,'%H:%M:%S').time())
    for ii,ch in enumerate(order):
        name='Ch{}'.format(ch)
        s[name] = [[x,y,z] for x,y,z in zip(s[ii+2], s[ii+12], s[ii+22])] # not work for n=/=10
    
    return s

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

def GetPressure(ObsDateTime, PData, PressureOption, NumPressurePoints, NumValidPresPoints, DailyMeanPressure, DefaultSurfacePressure):
    if PressureOption == 0 or NumPressurePoints == 0 or NumValidPresPoints < MinValidPresPoints:
        GetPressure = DefaultSurfacePressure
    elif PressureOption == 1 or PressureOption == 3:
        I = 0
        TimeDiff1 = ObsDateTime - PData.iloc[I].DateTime
        I = NumPressurePoints - 1
        TimeDiff2 = ObsDateTime - PData.iloc[I].DateTime
        if TimeDiff1.total_seconds() < 0 or TimeDiff2.total_seconds() > 0:
            GetPressure = DailyMeanPressure
            print("Observation outside range of pressure tabulation: assuming daily mean")
            print("First time on file=", PData.iloc[0].DateTime)
            print("Last time  on file=", PData.iloc[NumPressurePoints - 1].DateTime)
            print("ObsTime           =", ObsDateTime)
        else:
            IPresMatch = 0
            for I in range(NumPressurePoints):
                TimeDiff = PData.iloc[I].DateTime - ObsDateTime
                if TimeDiff.total_seconds() > 0:
                    IPresMatch = I
                    break
            if IPresMatch == 0:
                print("Failed to match pressure time to obs time")
                print("ObsDateTime   =", ObsDateTime)
                print("PDateTime(1)  =", PData.iloc[0].DateTime)
                print("PDateTime(N)  =", PData.iloc[NumPressurePoints - 1].DateTime)

            GetPressure = PData.iloc[IPresMatch].Pressure

            if GetPressure < 0:
                GetPressure = DailyMeanPressure

    return GetPressure

def Rayleigh(Dbug, Wavelength, SurfacePressure):
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

    ColumnDensity = AVOGADRO_CONSTANT * SurfacePressure / (MOLECULAR_WEIGHT_DRY_AIR * GRAVITY * 100.0)
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
    
    return RayleighOD

def xlin(N, X, Y, Z):
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

def GetAirMass(Index, ApparentZenith):
    ScaleHeight = [0.0, 6.58, 1.0, 22.0]
    EarthRadius = 6370.0

    X = math.cos(ApparentZenith * math.pi / 180)
    R = ScaleHeight[Index] / EarthRadius
    D = math.sqrt(X * X + 2 * R + R * R)

    GetAirMass = (1 + R) / D

    return GetAirMass

def CheckTripletCv(Dbug, NumPoints, Nstart, Airmass, CoVar):
    MinAirmass = 2
    MaxAirmass = 6
    MinPtsUnitAirmass = 2
    MaxCoVar = 1.0
    Index = np.full(len(Airmass),None)
    NumOK = 0
    NumPtsUnitAirmass = [0] * (MaxAirmass + 1)

    for I in range(Nstart, Nstart + NumPoints):
        #print('am[i] = ',Airmass[I])
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

    return SpreadFlag,NumOK,Index

