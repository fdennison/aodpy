class SatelliteAxes:
    def __init__(self, r=0.0, s=0.0, t=(0.0, 0.0, 0.0), u=(0.0, 0.0, 0.0),
                 v=(0.0, 0.0, 0.0), w=(0.0, 0.0, 0.0)):
        self.r = r
        self.s = s
        self.t = t
        self.u = u
        self.v = v
        self.w = w


sat = SatelliteAxes()

sat.u = (1.0,2.0,3.0)
sat.r = 2.0

sat.u = tuple(x / sat.r for x in sat.u)
[print(x) for x in dir(sat) if not x.startswith('__')]

