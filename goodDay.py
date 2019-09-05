from enum import IntEnum
import swisseph as swe
from math import floor
import skyfield.api
import skyfield
from skyfield import almanac
from datetime import datetime, time, date, timezone, timedelta
import numpy as np
import dateutil.relativedelta as relativedelta
import dateutil.rrule as rrule
from scipy.optimize import minimize_scalar


class Nakshatras(IntEnum):
    __order__ = "ASHWINI BHARANI KRITTIKAI ROHINI MIRGASIRISHAM THIRUVATHIRAI PUNARPOOSAM POOSAM AAYILYAM MAGAM " \
                "POORAM UTHIRAM HASTHAM CHITIRAI SWATHI VISAKAM ANUSHAM KETAI MOOLAM POORADAM UTHTHARADAM THIRUVOONAM" \
                " AVITAM SATHAYAM POORATTATHI UTHTHRATTATHI REVATHI"
    ASHWINI = 1
    BHARANI = 2
    KRITTIKAI = 3
    ROHINI = 4
    MIRGASIRISHAM = 5
    THIRUVATHIRAI = 6
    PUNARPOOSAM = 7
    POOSAM = 8
    AAYILYAM = 9
    MAGAM = 10
    POORAM = 11
    UTHIRAM = 12
    HASTHAM = 13
    CHITIRAI = 14
    SWATHI = 15
    VISAKAM = 16
    ANUSHAM = 17
    KETAI = 18
    MOOLAM = 19
    POORADAM = 20
    UTHTHARADAM = 21
    THIRUVOONAM = 22
    AVITAM = 23
    SATHAYAM = 24
    POORATTATHI = 25
    UTHTHRATTATHI = 26
    REVATHI = 27


class Varam(IntEnum):
    __order__ = "SUN MON TUE WED THU FRI SAT"
    SUN = 1
    MON = 2
    TUE = 3
    WED = 4
    THU = 5
    FRI = 6
    SAT = 7

    @staticmethod
    def get_relative_day(day):
        table = {Varam.SUN: relativedelta.SU,
                 Varam.MON: relativedelta.MO,
                 Varam.TUE: relativedelta.TU,
                 Varam.WED: relativedelta.WE,
                 Varam.THU: relativedelta.TH,
                 Varam.FRI: relativedelta.FR,
                 Varam.SAT: relativedelta.SA
                 }
        return table[day]


class Rasi(IntEnum):
    __order__ = "MESHAM RISHABAM MITHUNAM KADAKAM SIMHAM KANNI THULAM VRICHIKAM DHANUSH MAKARAM KUMBAM MEENAM"
    MESHAM = 1
    RISHABAM = 2
    MITHUNAM = 3
    KADAKAM = 4
    SIMHAM = 5
    KANNI = 6
    THULAM = 7
    VRICHIKAM = 8
    DHANUSH = 9
    MAKARAM = 10
    KUMBAM = 11
    MEENAM = 12


class Thithi(IntEnum):
    SUK_PRATHAMAI = 1
    SUK_DVITHIYAI = 2
    SUK_THRITHIYAI = 3
    SUK_CHATHURTHI = 4
    SUK_PANCHAMI = 5
    SUK_SHASTI = 6
    SUK_SAPTHAMI = 7
    SUK_ASHTAMI = 8
    SUK_NAVAMI = 9
    SUK_DASAMI = 10
    SUK_YEKADASI = 11
    SUK_DVADASI = 12
    SUK_THRIYODASI = 13
    SUK_CHATHURDASI = 14
    POURNAMI = 15
    KRI_PRATHAMAI = 16
    KRI_DVITHIYAI = 17
    KRI_THRITHIYAI = 18
    KRI_CHATHURTHI = 19
    KRI_PANCHAMI = 20
    KRI_SHASTI = 21
    KRI_SAPTHAMI = 22
    KRI_ASHTAMI = 23
    KRI_NAVAMI = 24
    KRI_DASAMI = 25
    KRI_YEKADASI = 26
    KRI_DVADASI = 27
    KRI_THRIYODASI = 28
    KRI_CHATHURDASI = 29
    AMMAVASAI = 30


class Months(IntEnum):
    __order__ = "MESHAM RISHABAM MITHUNAM KADAKAM SIMHAM KANNI THULAM VRICHIKAM DHANUSH MAKARAM KUMBAM MEENAM"
    MESHAM = 1
    RISHABAM = 2
    MITHUNAM = 3
    KADAKAM = 4
    SIMHAM = 5
    KANNI = 6
    THULAM = 7
    VRICHIKAM = 8
    DHANUSH = 9
    MAKARAM = 10
    KUMBAM = 11
    MEENAM = 12


class AmruthathiYogams(IntEnum):
    SIDHA = 1
    AMRUTHA = 2
    MARANA = 3


class TimeRange:
    def __init__(self, start_hr, start_min, stop_hr, stop_min):
        self.start = time(start_hr, start_min)
        self.stop = time(stop_hr, stop_min)


class DateRange:
    def __init__(self, start_date: date, stop_date: date):
        self.start = start_date
        self.stop = stop_date


class Directions(IntEnum):
    NORTH = 1
    EAST = 2
    SOUTH = 3
    WEST = 4
    NORTH_EAST = 5
    SOUTH_EAST = 6
    NORTH_WEST = 7
    SOUTH_WEST = 8


marana_yoga_table = {
    Varam.SUN: [Nakshatras.MAGAM, Nakshatras.VISAKAM, Nakshatras.ANUSHAM, Nakshatras.KETAI, Nakshatras.AVITAM],
    Varam.MON: [Nakshatras.KRITTIKAI, Nakshatras.MAGAM, Nakshatras.VISAKAM, Nakshatras.UTHTHARADAM,
                Nakshatras.POORATTATHI],
    Varam.TUE: [Nakshatras.THIRUVATHIRAI, Nakshatras.VISAKAM, Nakshatras.SATHAYAM, Nakshatras.POORATTATHI],
    Varam.WED: [Nakshatras.ASHWINI, Nakshatras.HASTHAM, Nakshatras.MOOLAM, Nakshatras.AVITAM, Nakshatras.REVATHI],
    Varam.THU: [Nakshatras.KRITTIKAI, Nakshatras.ROHINI, Nakshatras.MIRGASIRISHAM, Nakshatras.THIRUVATHIRAI,
                Nakshatras.UTHIRAM, Nakshatras.SATHAYAM],
    Varam.FRI: [Nakshatras.ROHINI, Nakshatras.POOSAM, Nakshatras.AAYILYAM, Nakshatras.MAGAM, Nakshatras.KETAI,
                Nakshatras.THIRUVOONAM],
    Varam.SAT: [Nakshatras.AAYILYAM, Nakshatras.UTHIRAM, Nakshatras.HASTHAM, Nakshatras.CHITIRAI,
                Nakshatras.POORATTATHI, Nakshatras.REVATHI]
}

soolai_table = {
    Varam.SUN: [Directions.WEST, Directions.NORTH_WEST],
    Varam.MON: [Directions.EAST, Directions.SOUTH_WEST],
    Varam.TUE: [Directions.NORTH, Directions.NORTH_WEST],
    Varam.WED: [Directions.NORTH, Directions.NORTH_EAST],
    Varam.THU: [Directions.SOUTH, Directions.SOUTH_EAST],
    Varam.FRI: [Directions.WEST, Directions.SOUTH_WEST],
    Varam.SAT: [Directions.EAST, Directions.SOUTH_WEST]
}

avoid_table = {
    Varam.SUN: ([Thithi.KRI_PANCHAMI, Thithi.SUK_PANCHAMI], Nakshatras.HASTHAM),
    Varam.MON: ([Thithi.SUK_SAPTHAMI, Thithi.KRI_SAPTHAMI], Nakshatras.THIRUVOONAM),
    Varam.TUE: ([Thithi.SUK_SAPTHAMI, Thithi.KRI_SAPTHAMI], Nakshatras.ASHWINI),
    Varam.WED: ([Thithi.KRI_ASHTAMI, Thithi.SUK_ASHTAMI], Nakshatras.ANUSHAM),
    Varam.THU: ([Thithi.KRI_THRITHIYAI, Thithi.SUK_THRITHIYAI], Nakshatras.POOSAM),
    Varam.FRI: ([Thithi.KRI_NAVAMI, Thithi.SUK_NAVAMI], Nakshatras.REVATHI),
    Varam.SAT: ([Thithi.KRI_YEKADASI, Thithi.SUK_YEKADASI], Nakshatras.REVATHI)
}

kari_naal = {
    Months.MESHAM: [6, 15],
    Months.RISHABAM: [7, 16, 17],
    Months.MITHUNAM: [1, 6],
    Months.KADAKAM: [2, 10, 20],
    Months.SIMHAM: [2, 9, 28],
    Months.KANNI: [16, 29],
    Months.THULAM: [6, 20],
    Months.VRICHIKAM: [1, 10, 17],
    Months.DHANUSH: [6, 9, 11],
    Months.MAKARAM: [1, 2, 3, 11, 17],
    Months.KUMBAM: [15, 16, 17],
    Months.MEENAM: [6, 15, 19]
}

dhaniya_naal = {
    Months.MESHAM: [3, 20],
    Months.RISHABAM: [9, 22],
    Months.MITHUNAM: [8, 22],
    Months.KADAKAM: [7, 20],
    Months.SIMHAM: [7, 18],
    Months.KANNI: [9, 26],
    Months.THULAM: [8, 10],
    Months.VRICHIKAM: [8, 14],
    Months.DHANUSH: [8, 26],
    Months.MAKARAM: [8, 15],
    Months.KUMBAM: [15, 24],
    Months.MEENAM: [18, 24]
}

###                 set id, curr stars,                                             compatible_set
tara_palan_table = [(1, [Nakshatras.ASHWINI, Nakshatras.MAGAM, Nakshatras.MOOLAM], [2, 4, 6, 8, 9]),
                    (2, [Nakshatras.BHARANI, Nakshatras.POORAM, Nakshatras.POORADAM], [1, 3, 5, 7, 9]),
                    (3, [Nakshatras.KRITTIKAI, Nakshatras.UTHIRAM, Nakshatras.UTHTHARADAM], [1, 2, 4, 6, 8]),
                    (4, [Nakshatras.ROHINI, Nakshatras.HASTHAM, Nakshatras.THIRUVOONAM], [2, 3, 5, 7, 9]),
                    (5, [Nakshatras.MIRGASIRISHAM, Nakshatras.CHITIRAI, Nakshatras.AVITAM], [1, 3, 4, 6, 8]),
                    (6, [Nakshatras.THIRUVATHIRAI, Nakshatras.SWATHI, Nakshatras.SATHAYAM], [2, 4, 5, 7, 9]),
                    (7, [Nakshatras.PUNARPOOSAM, Nakshatras.VISAKAM, Nakshatras.POORADAM], [1, 3, 5, 6, 8]),
                    (8, [Nakshatras.POOSAM, Nakshatras.ANUSHAM, Nakshatras.UTHTHRATTATHI], [2, 4, 6, 7, 9]),
                    (9, [Nakshatras.AAYILYAM, Nakshatras.KETAI, Nakshatras.REVATHI], [1, 3, 5, 7, 8])]
nakshaktra_ignore_list = [Nakshatras.BHARANI, Nakshatras.KRITTIKAI, Nakshatras.THIRUVATHIRAI, Nakshatras.AAYILYAM,
                          Nakshatras.POORADAM, Nakshatras.POORAM, Nakshatras.POORATTATHI, Nakshatras.KETAI,
                          Nakshatras.VISAKAM, Nakshatras.CHITIRAI, Nakshatras.SWATHI, Nakshatras.MAGAM]
good_horai = {Varam.SUN: [TimeRange(7, 30, 10, 0), TimeRange(14, 0, 16, 30), TimeRange(21, 0, 23, 59)],
              Varam.MON: [TimeRange(6, 0, 7, 0), TimeRange(12, 0, 14, 0), TimeRange(18, 0, 21, 0),
                          TimeRange(22, 0, 23, 0)],
              Varam.TUE: [TimeRange(10, 30, 11, 0), TimeRange(12, 0, 13, 0), TimeRange(16, 30, 18, 0),
                          TimeRange(19, 0, 20, 0)],
              Varam.WED: [TimeRange(9, 0, 10, 0), TimeRange(13, 30, 15, 0), TimeRange(16, 0, 17, 0),
                          TimeRange(21, 0, 22, 0), TimeRange(23, 0, 23, 59)],
              Varam.THU: [TimeRange(9, 0, 10, 30), TimeRange(13, 0, 13, 30), TimeRange(16, 30, 18, 0),
                          TimeRange(18, 0, 19, 0), TimeRange(20, 0, 21, 0)],
              Varam.FRI: [TimeRange(6, 0, 9, 0), TimeRange(13, 0, 13, 30), TimeRange(17, 0, 18, 0),
                          TimeRange(20, 0, 21, 0), TimeRange(22, 30, 23, 0)],
              Varam.SAT: [TimeRange(7, 0, 7, 30), TimeRange(10, 30, 13, 0), TimeRange(17, 0, 19, 30),
                          TimeRange(21, 0, 22, 0)]}

tithi_list = [t for t in Thithi]
rasi_list = [r for r in Rasi]
varam_list = [v for v in Varam]
month_list = [m for m in Months]
nakshatras_deg_span = 360. / 27

deg_span = 12
rasi_deg_span = 30
yoga_deg_span = 6


class Obj:
    start = 0.
    stop = 0.


class GoodDay:
    @staticmethod
    def get_good_Tharai(inp: Nakshatras):
        nakshatrasEnum = [k for k in Nakshatras]
        nakshatras = [int(k) for k in Nakshatras]
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


class Almanac:
    def __init__(self, latitude='13.0827 N', longitude='80.2707 E'):
        self.latitude = latitude
        self.longitude = longitude
        self.time_scale, self.planets = Almanac.init_skyfield()
        self.tithi = None
        self.vara = None
        self.nakshaktra = None
        self.yoga = None
        self.karanam = None
        swe.set_sid_mode(swe.SIDM_LAHIRI)

    @staticmethod
    def init_skyfield():
        ts = skyfield.api.load.timescale()
        planets = skyfield.api.load('de421.bsp')
        return ts, planets

    def compute(self, in_date: datetime):
        in_date = in_date.astimezone(timezone.utc)
        _, in_date_jul = swe.utc_to_jd(in_date.year, in_date.month, in_date.day, in_date.hour, in_date.minute,
                                       in_date.second, swe.GREG_CAL)
        # print(in_date_jul)

        # in_date_jul = swe.julday(in_date.year, in_date.month, in_date.day,
        #                          in_date.hour + in_date.minute / 60 + in_date.second / 3600)
        sun_pos = swe.calc(in_date_jul, swe.SUN, swe.FLG_SIDEREAL)[0]
        moon_pos = swe.calc(in_date_jul, swe.MOON, swe.FLG_SIDEREAL)[0]
        self.tithi, self.pos_deg = self._calc_tithi(sun_pos, moon_pos)

    def _calc_tithi(self, sun_pos: float, moon_pos: float) -> Thithi:
        if moon_pos < sun_pos:
            moon_pos += 360.
        delta = moon_pos - sun_pos
        idx = int(floor(delta / 12))
        return tithi_list[idx], delta

    def get_dates_on_day(self, day: Varam, range: DateRange):
        rr = rrule.rrule(rrule.WEEKLY, byweekday=Varam.get_relative_day(day),
                         dtstart=datetime.fromordinal(range.start.toordinal()), cache=True)
        print(rr.between(after=datetime.fromordinal(range.start.toordinal()),
                         before=datetime.fromordinal(range.stop.toordinal()), inc=True))

        # GoodDay.get_good_Tharai(Nakshatras.THIRUVATHIRAI)

    def _get_date_time_from_float(self, x: float):
        d_int = int(floor(x))
        d = date.fromordinal(d_int)
        t_float = x - d_int
        hr_float = t_float * 24
        hr = int(floor(hr_float))
        minute = int(floor((hr_float - hr) * 60))
        t = time(hr, minute)
        dt = datetime(d.year, d.month, d.day, hr, minute, tzinfo=timezone.utc)
        return dt

    def get_next_thiti(self, target: Thithi):
        low_deg = (int(target - 1)) * 12.
        target_deg = low_deg
        now = datetime.now(timezone.utc)
        self.compute(now)
        first_target = False
        if target_deg < self.pos_deg:
            orig_target = target_deg
            target_deg = 360
            first_target = True

        def calc(x: float):
            dt = self._get_date_time_from_float(x)
            self.compute(dt)
            # print(f'curr_pos {self.pos_deg} target_deg {target_deg}')
            return np.abs(self.pos_deg - target_deg)

        if first_target:
            result = minimize_scalar(calc, bounds=(now.toordinal(), now.toordinal() + 30), method='Bounded',
                                     options={'disp': False})
            now = self._get_date_time_from_float(result.x)
            target_deg = orig_target
        result = minimize_scalar(calc, bounds=(now.toordinal(), now.toordinal() + 30), method='Bounded',
                                 options={'disp': False})

        dt = self._get_date_time_from_float(result.x)

        self.compute(dt)
        print(dt, dt.astimezone(timezone(timedelta(hours=5, minutes=30))), self.tithi, self.pos_deg,
              int(self.tithi) * 12 - self.pos_deg, 'left')


a = Almanac()
a.compute(datetime.now(timezone.utc))
print(a.tithi)
a.get_next_thiti(Thithi.SUK_PANCHAMI)
a.get_dates_on_day(Varam.TUE, DateRange(date(2019, 9, 1), date(2019, 12, 31)))
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
t0 = ts.utc(2019, 8, 25)
t1 = ts.utc(2019, 8, 26)
t, y = almanac.find_discrete(t0, t1, almanac.sunrise_sunset(planets, chennai))

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

print('sun', deg_min_sec(sun[0]), int(floor(sun[0] / 30)) + 1)
print('moon', deg_min_sec(moon[0]), int(floor(moon[0] / 30)) + 1)
print('mars', deg_min_sec(mars[0]), int(floor(mars[0] / 30)) + 1)
print('mercury', deg_min_sec(mercury[0]), int(floor(mercury[0] / 30)) + 1)
from astral import Astral, Location, AstralGeocoder
