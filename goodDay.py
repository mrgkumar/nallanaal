import swisseph as swe
from math import floor
import skyfield.api
import skyfield
from datetime import datetime, date, timezone
import numpy as np
import almanac
from almanac import Almanac, find_discrete
from params import NAKSHATRAS, VARAM, RASI, THITHI, MONTHS, DateRange

tithi_list = [t for t in THITHI]
rasi_list = [r for r in RASI]
varam_list = [v for v in VARAM]
month_list = [m for m in MONTHS]
nakshatras_deg_span = 360. / 27

deg_span = 12
rasi_deg_span = 30
yoga_deg_span = 6


class Obj:
    start = 0.
    stop = 0.


class GoodDay:
    @staticmethod
    def get_good_Tharai(inp: NAKSHATRAS):
        nakshatrasEnum = [k for k in NAKSHATRAS]
        nakshatras = [int(k) for k in NAKSHATRAS]
        startIdx = np.squeeze(np.argwhere(np.array(nakshatras) == int(inp)))
        origList = nakshatras[startIdx:] + nakshatras[:startIdx]
        nakshatrasEnum = nakshatrasEnum[startIdx:] + nakshatrasEnum[:startIdx]
        tmp = list(np.array(nakshatras[startIdx:] + list(np.array(nakshatras[:startIdx]) + 27)) - int(inp))
        modTmp = [(t % 9) + 1 for t in tmp]
        print(modTmp)
        print(origList)
        print(nakshatrasEnum)


def deg_min_sec(dec, aya=None):
    if aya is not None:
        dec = dec - aya
    deg = int(floor(dec))
    minutes = dec % 1.0 * 60
    secs = minutes % 1.0 * 60
    return f'{deg}-{int(floor(minutes))}-{secs}'


a = Almanac()
a.compute(datetime.now(timezone.utc))
print(a.tithi)
for t in tithi_list:
    print(f'requesting for {t.__str__()}')
    a.get_next_thiti(t)


# a.get_dates_on_day(VARAM.TUE, DateRange(date(2021, 9, 1), date(2025, 12, 31)))
# from flatlib.datetime import Datetime
# from flatlib.geopos import GeoPos
# import flatlib.const as const
# from flatlib.chart import Chart

# now = Datetime('2009/07/15', '05:30')
# pos = GeoPos(13.0827, 80.2707)  # chennai
now = swe.julday(2019, 8, 26, .5, swe.GREG_CAL)

ts = skyfield.api.load.timescale()
planets = skyfield.api.load('de421.bsp')

# for k in range(0,39):
#     swe.set_sid_mode(k)
#     ayanamsa = swe.get_ayanamsa(now)
#     print(swe.get_ayanamsa_name(k),ayanamsa)
swe.set_sid_mode(swe.SIDM_LAHIRI)
# chart = Chart(now, pos, hsys=const.HOUSES_EQUAL)
# sun = chart.getObject(const.SUN)
# moon = chart.getObject(const.MOON)
chennai = skyfield.api.Topos('13.0827 N', '80.2707 E')
t0 = ts.utc(2021, 9, 1)
t1 = ts.utc(2025, 12, 31)
t, y = find_discrete(t0, t1, almanac.sunrise_sunset(planets, chennai))

now = t[0].tt
ayanamsa = swe.get_ayanamsa(now)

sun = swe.calc(now, swe.SUN, swe.FLG_SIDEREAL)
moon = swe.calc(now, swe.MOON, swe.FLG_SIDEREAL)
mars = swe.calc(now, swe.MARS, swe.FLG_SIDEREAL)
mercury = swe.calc(now, swe.MERCURY, swe.FLG_SIDEREAL)
jupiter = swe.calc(now, swe.JUPITER, swe.FLG_SIDEREAL)
venus = swe.calc(now, swe.VENUS, swe.FLG_SIDEREAL)
saturn = swe.calc(now, swe.SATURN, swe.FLG_SIDEREAL)
rahu = swe.calc(now, swe.MEAN_NODE, swe.FLG_SIDEREAL)
ketu = swe.calc(now, swe.TRUE_NODE, swe.FLG_SIDEREAL)
print('ayanamsa', ayanamsa, deg_min_sec(ayanamsa))

# print('sun', deg_min_sec(sun[0]), int(floor(sun[0] / 30)) + 1)
# print('moon', deg_min_sec(moon[0]), int(floor(moon[0] / 30)) + 1)
# print('mars', deg_min_sec(mars[0]), int(floor(mars[0] / 30)) + 1)
# print('mercury', deg_min_sec(mercury[0]), int(floor(mercury[0] / 30)) + 1)
