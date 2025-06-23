import math
import random

def correlation_coefficient(X, Y):
    """
    Computes the correlation coefficient of vectors X and Y.
    """
    NX = len(X)
    NY = len(Y)

    if NX <= 0 or NY <= 0 or NX != NY:
        print("Error in correlation coefficient subroutine!")
        print(f"Size of X vector        = {NX}")
        print(f"Size of Y vector        = {NY}")
        print("Correlation coefficient set to zero.")
        input()
        return 0.0

    mean_X = sum(X) / NX
    mean_Y = sum(Y) / NY
    residual_X = [x - mean_X for x in X]
    residual_Y = [y - mean_Y for y in Y]
    cov_XX = sum(x ** 2 for x in residual_X)
    cov_XY = sum(x * y for x, y in zip(residual_X, residual_Y))
    cov_YY = sum(y ** 2 for y in residual_Y)

    if cov_XX <= 0 or cov_YY <= 0:
        print("Error in correlation coefficient subroutine!")
        print(f"Autocovariance of X     = {cov_XX}")
        print(f"Autocovariance of Y     = {cov_YY}")
        print("Correlation coefficient set to zero.")
        input()
        return 0.0

    return cov_XY / math.sqrt(cov_XX * cov_YY)

def sample_uniform(mean, std_deviation):
    """
    Draws a sample from a uniform distribution with the given mean and standard deviation.
    """
    x = random.uniform(0, 1)
    return mean + std_deviation * (2 * x - 1)

def sample_uniform_relative(mean, std_deviation):
    """
    Draws a sample from a uniform distribution with the given mean and relative standard deviation.
    """
    x = random.uniform(0, 1)
    return mean * (1 + std_deviation * (2 * x - 1))

def linear_fit(X, Y):
    """
    Fits a linear model to the data points in X and Y.
    Returns the slope, intercept, and residual sum of squares.
    """
    N = len(X)
    sx, sy, sxx, sxy = sum(X), sum(Y), sum(x ** 2 for x in X), sum(x * y for x, y in zip(X, Y))
    det = N * sxx - sx ** 2
    slope = (sxy * N - sx * sy) / det
    intercept = (sxx * sy - sx * sxy) / det
    residual_sum_squares = sum((intercept + slope * x - y) ** 2 for x, y in zip(X, Y))
    return slope, intercept, residual_sum_squares

def weighted_linear_fit(self, W, X, Y):
    """
    Fits a linear model to the weighted data points in W, X, and Y.
    Returns the slope, intercept, residual weighted sum of squares,
    RMS errors for slope and intercept, and RMS deviation.
    """
    N = len(X)
    sw = sum(W)
    swx = sum(w * x for w, x in zip(W, X))
    swy = sum(w * y for w, y in zip(W, Y))
    swxx = sum(w * x ** 2 for w, x in zip(W, X))
    swxy = sum(w * x * y for w, x, y in zip(W, X, Y))
    det = sw * swxx - swx ** 2
    slope = (swxy * sw - swx * swy) / det
    intercept = (swxx * swy - swx * swxy) / det
    residual_weighted_sum_squares = sum(w * (intercept + slope * x - y) ** 2 for w, x, y in zip(W, X, Y))
    rms_deviation = math.sqrt(residual_weighted_sum_squares / sw)
    rms_error_intercept = math.sqrt((swxx / det) * (residual_weighted_sum_squares / (N - 2)))
    rms_error_slope = math.sqrt((sw / det) * (residual_weighted_sum_squares / (N - 2)))
    return slope, intercept, residual_weighted_sum_squares, rms_error_slope, rms_error_intercept, rms_deviation

def coefficient_variation(X):
    """
    Computes the coefficient of variation of samples in vector X.
    """
    mean = sum(X) / len(X)
    variance = sum((x - mean) ** 2 for x in X) / len(X)
    return math.sqrt(variance) / mean

def stat(data):
    """
    Returns the mean, standard deviation, minimum, and maximum of the input data.
    """
    ndat = len(data)
    if ndat > 0:
        mean = sum(data) / ndat
        sum_squared_deviations = sum((x - mean) ** 2 for x in data)
        sdev = math.sqrt(sum_squared_deviations / (ndat - 1)) if ndat > 1 else 0
        min_data = min(data)
        max_data = max(data)
        return mean, sdev, min_data, max_data
    else:
        print("Statflos: number of data points <= 0")
        return 0, 0, None, None
