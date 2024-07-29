import numpy as np
import math

            
def CheckTripletCv(Dbug, numPoints, nstart, Airmass, CoVar):
 # Checks a sequence of sun measurements for suitability
 # for Langley analysis, using the criteria that there must
 # be at least two data points
 # in each of airmass ranges 0-1, 1-2, 2-3, 3-4, 4-5 and 5-6
 # with triplet CV of less than 1% (aod error .002-.005).  

    MinAirmass = 2
    MaxAirmass = 6
    MinPtsUnitAirmass = 2
    MaxCoVar = 1.0
    index =[]
    numOk = 0
    numPtsUnitAirmass = [0] * (MaxAirmass + 1)

    for i in range(nstart, nstart + numPoints):
        #print('am[i] = ',Airmass[i])
        if MinAirmass <= Airmass[i] <= MaxAirmass and CoVar[i] <= MaxCoVar:
            numOk += 1
            index.append(i)
            j = int(Airmass[i])
            if j > MaxAirmass:
                print("Array overflow, airmass too big!")
                print("inT(Airmass)   :", j)
                print("Max Airmass    :", MaxAirmass)
                halt()
            numPtsUnitAirmass[j] += 1

    SpreadFlag = True
    for j in range(MinAirmass, MaxAirmass):
        if numPtsUnitAirmass[j] < MinPtsUnitAirmass:
            SpreadFlag = False
            print("insufficient points in airmass interval", j)

    if Dbug:
        print("numOk    ", numOk)
        for j in range(1, 7):
            print("Points in Airmassrange", j, numPtsUnitAirmass[j])

    return SpreadFlag,numOk,index




def CheckFitQuality(Dbug, numPoints, X, Y, nref, MaxSdevFit):
 # Checks a sequence of sun measurements for suitability
 # for Langley analysis, using the following criteria:
 # 1. There must be at least two data points
 #    in each of airmass ranges 0-1, 1-2, 2-3, 3-4, 4-5 and 5-6
 #    with triplet CV of less than 1% (aod error .002-.005).
 # 2. if this test is passed, a preliminary linear regression
 #    of lnV against airmass is performed. Then the points
 #    selected on the basis of step 1 are reexamined.
 #    Any point lying more than 1.5 standard deviations from
 #    the regression line is rejected.
 # 3. The number of points in each airmass interval is recorded,
 #    and the period is deemed satisfactory if the criterion
 #    in step 1 are met.
 # 4. The standard deviation of the fit must be less than 0.005,
 #    corresponding to an error of about 0.002 in aod for m=2.    
    
    MinAirmass = 2
    MaxAirmass = 6
    MinPtsUnitAirmass = 4
    index=[]
    n = nref
    intercept, Slope, Residual, Erms, Delintercept, DelSlope = boxfit(numPoints, X[:, n], Y[:, n])

    if Dbug:
        print("num points in fit=", numPoints)
        print("intercept         ", intercept)
        print("Slope             ", Slope)
        print("Residual          ", Residual)
        print("Erms              ", Erms)
        print("Delintercept      ", Delintercept)
        print("DelSlope          ", DelSlope)

    numPtsUnitAirmass = [0] * (MaxAirmass + 1)
    i = 0

    for j in range(numPoints):
        Error = Y[j, nref] - (intercept + Slope * X[j, nref])
        if abs(Error) < 1.5 * Erms:
            i += 1
            X[i - 1, :] = X[j, :]
            Y[i - 1, :] = Y[j, :]
            index.append(j)
            k = int(X[i - 1, nref])
            numPtsUnitAirmass[k] += 1
        else:
            print("Rejecting point, am=:", X[j, nref], "LnV=", Y[j, nref], "Error=", Error, "1.5sigma=", 1.5 * Erms)

    print("Points satisfying Pass 1 test=", numPoints)
    numPoints = i
    print("Points satisfying Pass 2 test=", numPoints)

    SpreadFlag = True
    for j in range(MinAirmass, MaxAirmass):
        if numPtsUnitAirmass[j] < MinPtsUnitAirmass:
            SpreadFlag = False
            print("insufficient points in airmass interval", j)

    if Dbug:
        for j in range(1, 7):
            print("Points in Airmassrange", j, numPtsUnitAirmass[j])

    if SpreadFlag:
        intercept, Slope, Residual, Erms, Delintercept, DelSlope = boxfit(numPoints, X[:numPoints, n], Y[:numPoints, n])
        if Dbug:
            print("num points in fit=", numPoints)
            print("intercept         ", intercept)
            print("Slope             ", Slope)
            print("Residual          ", Residual)
            print("Erms              ", Erms)
            print("Delintercept      ", Delintercept)
            print("DelSlope          ", DelSlope)

        print("Erms              ", Erms)
        print("MaxSdevFit        ", MaxSdevFit)
        if Erms < MaxSdevFit:
            FitFlag = True
        else:
            FitFlag = False
            print("FitFlag set false because standard deviation=", Erms)
    else:
        FitFlag = False
        
    return  SpreadFlag, FitFlag, index


def boxfit(n, X, Y):    
#      Subroutine Boxfit fits a linear model to n data points
#      stored in arrays X and Y. it is unweigted in the context of
#      fitting Langley plots where lnV is plotted against airmass m.
#      A standard least-squares fit automatically weights the fit
#      increasingly toward high airmass. The present algorithm was
#      determined to minimise deviations in optical depth from the mean.

#      Reference: Hermann, Box, Reagan and Evans, "Alternative approach
#      to the analysis of solar photometer data", Applied Optics, V20, 
#      pp. 2925--2928, 1981.

#      The slope of the line is returned in S and the intercept in A.
#      The residual sum of squares is returned in R.
#      The RMS deviation is returned in E.
#      The standard deviation of the intercept is returned in EA.
#      The standard deviation of the slope     is returned in ES.    
    
    S0 = 0
    SX1 = 0
    SX2 = 0
    SX1Y = 0
    SX2Y = 0

    for i in range(n):
        S0 += 1
        SX1 += 1 / X[i]
        SX2 += 1 / (X[i] * X[i])
        SX1Y += Y[i] / X[i]
        SX2Y += Y[i] / (X[i] * X[i])

    DET = S0 * SX2 - SX1 * SX1
    A = (S0 * SX2Y - SX1 * SX1Y) / DET
    S = (SX2 * SX1Y - SX1 * SX2Y) / DET

    R = 0
    T = 0
    for i in range(n):
        R += (A + S * X[i] - Y[i]) ** 2
        T += ((A + S * X[i] - Y[i]) / X[i]) ** 2

    E = math.sqrt(R / (n - 2))
    EA = math.sqrt((S0 / DET) * (T / (n - 2)))
    ES = math.sqrt((SX2 / DET) * (T / (n - 2)))  
    
    return A, S, R, E, EA, ES

def elfit(N, W, X, Y):
    # Subroutine ELFIT fits a linear model to N data points stored in
    # arrays X and Y.
    # The data is weighted by weights in array W.

    # The slope of the line is returned in S and the intercept in A.
    # The residual weighted sum of squares is returned in R.
    # The RMS deviation is returned in E.
    # The RMS error associated with the intercept is returned in EA.
    # The RMS error associated with the slope     is returned in ES.    
    
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
    
    return A, S , R, E #, EA, ES
