import math
pi = math.pi
rad2deg = 180.0 / pi
deg2rad = pi /180.0
avogadro = 6.022529e+23  # mol^{-1}
M_dryair = 2.8964e-02  # kg/mol
earthrad = 6370.0 #km
gravity = 9.80665  # ms^{-2}, standard value


def getairmass(idx, appsolzen):
    #Finds airmass components given apparent solar zenith angle in degrees.
    #     Based on Forgan 1988, Appl Opt 27, 2546-2548.
    #     Components are:
    #     1. Molecular atmosphere, scale height 6.58 km
    #     2. Aerosol   atmosphere, scale height 1.0  km
    #     3. Ozone layer,          scale height 22.0 km
    
    scaleheight = [0.0, 6.58, 1.0, 22]
    x = math.cos(appsolzen * deg2rad)
    r = scaleheight[idx]/earthrad
    d = math.sqrt(x*x + 2*r + r*r)
    airmass = (1+r)/d
    return airmass


def rayleigh(Wavelength, SurfacePressure):
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

    ColumnDensity = avogadro * SurfacePressure / (M_dryair * gravity * 100.0)
    Depolarization = xlin(NumWavTab, WavTab, DepolTab, Wavelength)
    KingFactor = (6.0 + 3.0 * Depolarization) / (6.0 - 7.0 * Depolarization)
    Wcm = Wavelength * 1.0e-7
    RayleighCsa = 24.0 * pi**3 * (RefIndex**2 - 1.0)**2 / (Wcm**4 * DensitySTP**2 * (RefIndex**2 + 2.0)**2) * KingFactor
    RayleighOD = ColumnDensity * RayleighCsa

    #print("Surface Pressure :", SurfacePressure)
    #print("Column density   :", ColumnDensity)
    #print("Wavelength       :", Wavelength)
    #print("Ref Index        :", RefIndex)
    #print("Depolarization   :", Depolarization)
    #print("KingFactor       :", KingFactor)
    #print("RayleighCSA      :", RayleighCsa)
    #print("RayleighOD       :", RayleighOD)
    
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


# x[model][channel]
ozonecoef = [
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0447, 0.0002, 0.0327, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0004, 0.0322, 0.0087, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0447, 0.0002, 0.0327, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0447, 0.0002, 0.0327, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0002, 0.0327, 0.0327, 0.0000, 0.0000, 0.0000],
    [0.0041, 0.0462, 0.0000, 0.0000, 0.0002, 0.0327, 0.0327, 0.0327, 0.0000, 0.0000]
]

wavelength = [
    [440.0, 670.0, 870.0, 1020.0, 0.0, 0.0, 0.0, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 870.1, 870.2, 870.3, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 340.0, 380.0, 500.0, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 380.0, 500.0, 778.0, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 340.0, 380.0, 500.0, 1021.0, 1246.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 340.0, 380.0, 500.0, 1021.0, 1640.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 412.0, 500.0, 550.0, 0.0, 0.0, 936.0],
    [440.0, 670.0, 870.0, 1020.0, 412.0, 500.0, 531.0, 550.0, 0.0, 936.0]
]
