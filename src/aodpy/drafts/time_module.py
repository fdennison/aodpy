from datetime import datetime, timedelta
from typing import NamedTuple

class Date(NamedTuple):
    year: int
    month: int
    day: int

class Time(NamedTuple):
    hour: int
    minute: int
    second: int

class DateTime(NamedTuple):
    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: int

def days_in_month(month: int, year: int) -> int:
    if month == 2:
        if (year % 4 == 0 and year % 100 != 0) or year % 400 == 0:
            return 29
        else:
            return 28
    else:
        return [31, None, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month - 1]

def day_of_year(day: int, month: int, year: int) -> int:
    day_of_year = sum(days_in_month(m, year) for m in range(1, month))
    return day_of_year + day

def dayofyear(date: Date) -> int:
    return day_of_year(date.day, date.month, date.year)

def daybreak(year: int, numday: int) -> tuple[int, int]:
    month = 1
    day = 1
    while numday > days_in_month(month, year):
        numday -= days_in_month(month, year)
        month += 1
    return numday, month

def num_days_year(year: int) -> int:
    return sum(days_in_month(month, year) for month in range(1, 13))

def daycount(date1: tuple[int, int, int], date2: tuple[int, int, int]) -> int:
    year1, month1, day1 = date1
    year2, month2, day2 = date2
    date1 = datetime(year1, month1, day1)
    date2 = datetime(year2, month2, day2)
    return (date2 - date1).days

def add_time(ymdhms: list[int], delta_sec: int) -> None:
    dt = datetime(*ymdhms[:3]) + timedelta(seconds=ymdhms[3] * 3600 + ymdhms[4] * 60 + ymdhms[5] + delta_sec)
    ymdhms[:] = [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]

def del_time(time1: list[int], time2: list[int]) -> int:
    dt1 = datetime(*time1[:3]) + timedelta(hours=time1[3], minutes=time1[4], seconds=time1[5])
    dt2 = datetime(*time2[:3]) + timedelta(hours=time2[3], minutes=time2[4], seconds=time2[5])
    return int((dt2 - dt1).total_seconds())

def hour2hms(hour: float) -> Time:
    sign = 1 if hour >= 0 else -1
    hour = abs(hour)
    h = int(hour)
    m = int((hour - h) * 60)
    s = int((hour - h - m / 60) * 3600)
    return Time(sign * h, sign * m, sign * s)

def hms2hour(hms: tuple[int, int, int]) -> float:
    h, m, s = hms
    sign = -1 if h < 0 or m < 0 or s < 0 else 1
    return sign * (abs(h) + abs(m) / 60 + abs(s) / 3600)

def timediff(time2: Time, time1: Time) -> int:
    return (time2.second - time1.second) + (time2.minute - time1.minute) * 60 + (time2.hour - time1.hour) * 3600

def datediff(date2: Date, date1: Date) -> int:
    start_date = min(date1, date2, key=lambda d: (d.year, d.month, d.day))
    end_date = max(date1, date2, key=lambda d: (d.year, d.month, d.day))
    sign = 1 if date2 >= date1 else -1

    if end_date.year == start_date.year:
        return sign * (dayofyear(end_date) - dayofyear(start_date))

    days = num_days_year(start_date.year) - dayofyear(start_date)
    for year in range(start_date.year + 1, end_date.year):
        days += num_days_year(year)
    days += dayofyear(end_date)
    return sign * days

def test_time_equal(time2: Time, time1: Time) -> bool:
    return (
        time2.hour == time1.hour
        and time2.minute == time1.minute
        and time2.second == time1.second
    )

def test_time_not_equal(time2: Time, time1: Time) -> bool:
    return (
        time2.hour != time1.hour
        or time2.minute != time1.minute
        or time2.second != time1.second
    )

def test_date_equal(date2: Date, date1: Date) -> bool:
    return (
        date2.year == date1.year
        and date2.month == date1.month
        and date2.day == date1.day
    )

def test_date_not_equal(date2: Date, date1: Date) -> bool:
    return (
        date2.year != date1.year
        or date2.month != date1.month
        or date2.day != date1.day
    )

def test_datetime_equal(datetime2: DateTime, datetime1: DateTime) -> bool:
    return (
        datetime2.year == datetime1.year
        and datetime2.month == datetime1.month
        and datetime2.day == datetime1.day
        and datetime2.hour == datetime1.hour
        and datetime2.minute == datetime1.minute
        and datetime2.second == datetime1.second
    )

def test_datetime_not_equal(datetime2: DateTime, datetime1: DateTime) -> bool:
    return (
        datetime2.year != datetime1.year
        or datetime2.month != datetime1.month
        or datetime2.day != datetime1.day
        or datetime2.hour != datetime1.hour
        or datetime2.minute != datetime1.minute
        or datetime2.second != datetime1.second
    )

def datetimediff(datetime2: DateTime, datetime1: DateTime) -> int:
    dt1 = datetime(datetime1.year, datetime1.month, datetime1.day, datetime1.hour, datetime1.minute, datetime1.second)
    dt2 = datetime(datetime2.year, datetime2.month, datetime2.day, datetime2.hour, datetime2.minute, datetime2.second)
    return int((dt2 - dt1).total_seconds())

def add_datetime(dt1: DateTime, dt2: DateTime) -> DateTime:
    dt1 = datetime(dt1.year, dt1.month, dt1.day, dt1.hour, dt1.minute, dt1.second)
    dt2 = datetime(dt2.year, dt2.month, dt2.day, dt2.hour, dt2.minute, dt2.second)
    result = dt1 + timedelta(
        days=dt2.day,
        hours=dt2.hour,
        minutes=dt2.minute,
        seconds=dt2.second,
    )
    return DateTime(result.year, result.month, result.day, result.hour, result.minute, result.second)
