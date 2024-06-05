import math

def Elapsed_Seconds(TimeIni, TimeEnd):
    Hour, Minute, Second = map(int, TimeIni[:6].split())
    MilliSecond = int(TimeIni[7:10])
    Ini = MilliSecond + 1000 * (Second + 60 * (Minute + 60 * Hour))

    Hour, Minute, Second = map(int, TimeEnd[:6].split())
    MilliSecond = int(TimeEnd[7:10])
    End = MilliSecond + 1000 * (Second + 60 * (Minute + 60 * Hour))

    Elapsed_Seconds = (End - Ini) / 1000
    return Elapsed_Seconds

def Bell(Pips):
    Beep = chr(7)
    for _ in range(Pips):
        print(Beep, end='')

def Halt(Message=None):
    if Message:
        print(Message)
    else:
        print('Halt!')
    input()

def String_to_Integer(String):
    Length = len(String.strip())
    String_to_Integer = 0
    for I in range(Length):
        Digit = ord(String[I]) - 48
        if 0 <= Digit <= 9:
            String_to_Integer += 10 ** (Length - I - 1) * Digit
        else:
            print(f"String_to_Integer: argument {String[I]} out of range")
            String_to_Integer = -1
            return String_to_Integer
    return String_to_Integer

def Integer_to_String(N, PadLength):
    Length = 8
    Integer_to_String = ' ' * Length
    I = 0
    X = N
    X += 0.1  # To prevent rounding problems
    while X >= 1.0:
        I += 1
        X /= 10.0
    if I > Length:
        print("Error in Integer_to_String!")
        print(f"More than {Length} digits in integer arg!")
        Halt()
    for J in range(I):
        X *= 10
        Digit = int(X)
        Integer_to_String = Integer_to_String[:J] + chr(Digit + 48) + Integer_to_String[J+1:]
        X -= Digit
    while I < PadLength:
        I += 1
        Integer_to_String = '0' + Integer_to_String
    return Integer_to_String.strip()

def Next_Unit_Number():
    UnitNumber = 1
    while True:
        try:
            with open(f"unit{UnitNumber}", "r"):
                UnitNumber += 1
        except FileNotFoundError:
            return UnitNumber

def Skip(DiscUnit, NumLines):
    I = 1
    Error = 0
    while I <= NumLines and Error == 0:
        try:
            next(DiscUnit)
            I += 1
        except StopIteration:
            Error = 1
    if Error != 0:
        print("Error reading in subroutine Skip!")
        input()

def Skip_Comments(DiscUnit, Comment='#'):
    LocalComment = Comment
    TestChar = next(DiscUnit, '').strip()
    while TestChar == LocalComment:
        TestChar = next(DiscUnit, '').strip()
    if TestChar:
        DiscUnit = iter([TestChar] + list(DiscUnit))
    else:
        FileName = getattr(DiscUnit, 'name', 'Unknown')
        print("Error in Skip_Comments!")
        print(f"File name = {FileName}")
        print(f"Unit      = {DiscUnit}")
        input()
    return DiscUnit

def Scan_File(FileName, BlockNumber):
    with open(FileName, 'r') as FileUnit:
        SkipLines = 0
        for _ in range(BlockNumber):
            Line = next(FileUnit, '').strip()
            while Line == '#':
                SkipLines += 1
                Line = next(FileUnit, '').strip()
            SkipLines += 1
            DataLines = 0
            while Line != '#':
                DataLines += 1
                SkipLines += 1
                Line = next(FileUnit, '').strip()
        SkipLines -= DataLines + 1
    return SkipLines, DataLines
