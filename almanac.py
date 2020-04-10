from datetime import datetime, timezone, date, time, timedelta
from math import floor
from typing import Tuple
import numpy as np
import skyfield.api
from skyfield.jpllib import SpiceKernel
from skyfield.timelib import Timescale, Time
from skyfield.units import Angle
from skyfield.toposlib import Topos
from skyfield.almanac import find_discrete, sunrise_sunset
import swisseph as swe
from dateutil import rrule as rrule
from scipy.optimize import minimize_scalar, fmin
from params import THITHI, VARAM, DateRange, RASI, MONTHS, NAKSHATRAS, Varam
from params import Thithi, Nakshaktram, Month
from params import tara_palan_table

tithi_list = Thithi._list
rasi_list = [r for r in RASI]
varam_list = [v for v in VARAM]
month_list = [m for m in MONTHS]


class Almanac:
    tara_palan_dict = dict()

    def __init__(self, latitude='12.9399348 N', longitude='77.7337622 E'):
        self.latitude = latitude
        self.longitude = longitude
        self.topo = Topos(latitude=latitude, longitude=longitude)
        self.time_scale, self.planets = Almanac.init_skyfield()
        self.vara = None
        self.yoga = None
        self.karanam = None
        self.earth = self.planets['earth']
        self.moon = self.planets['moon']
        self.sun = self.planets['sun']
        self.time = self.time_scale.now()
        self.ayanamsa = None
        self.moon_pos = None
        self.sun_pos = None
        self.thithi = None
        self.nakshaktra = None
        self.day = None

        self.compute_almanac()
        Almanac._compute_tara_palan_dict()

        swe.set_sid_mode(swe.SIDM_LAHIRI)

    def _compute_tara_palan_dict():
        tmp_dict = dict()
        for idx, nak, comp in tara_palan_table:
            tmp_dict[idx] = nak, comp
        for key, (nak, comp) in tmp_dict.items():
            for n in nak:
                for c in comp:
                    target_nak, _ = tmp_dict[c]
                    if n in Almanac.tara_palan_dict.keys():
                        Almanac.tara_palan_dict[n] += target_nak
                    else:
                        Almanac.tara_palan_dict[n] = [*target_nak]

    def compute_almanac(self):
        self.ayanamsa = Almanac._compute_ayanamsa(self.time)
        self.moon_pos = self._compute_moon_pos(self.time)
        self.sun_pos = self._compute_sun_pos(self.time)
        self.thithi = self._compute_tithi_on(self.time)
        self.nakshaktra = self._compute_nakshaktram(self.time)
        if len(self.time.shape) > 0:
            func = lambda x: Varam._dict[x.astimezone().weekday()]
            vfunc = np.vectorize(func, otypes=[VARAM])
            self.day = vfunc(self.time.utc_datetime())
        else:
            self.day = Varam._dict[self.time.utc_datetime().weekday()]
        self.time.utc_datetime()

    @staticmethod
    def init_skyfield() -> Tuple[Timescale, SpiceKernel]:
        ts = skyfield.api.load.timescale()
        planets = skyfield.api.load('de438s.bsp')  # ('de421.bsp')
        return ts, planets

    def compute(self, in_date: datetime):
        in_date = in_date.astimezone(timezone.utc)
        _, in_date_jul = swe.utc_to_jd(in_date.year, in_date.month, in_date.day, in_date.hour, in_date.minute,
                                       in_date.second, swe.GREG_CAL)
        # print(in_date_jul)

        # in_date_jul = swe.julday(in_date.year, in_date.month, in_date.day,
        #                          in_date.hour + in_date.minute / 60 + in_date.second / 3600)
        sun_pos, _ = swe.calc(in_date_jul, swe.SUN, swe.FLG_SIDEREAL)
        moon_pos, _ = swe.calc(in_date_jul, swe.MOON, swe.FLG_SIDEREAL)
        sun_pos = sun_pos[0]
        moon_pos = moon_pos[0]
        self.tithi, self.pos_deg = self._calc_tithi(sun_pos, moon_pos)

    def _calc_tithi(self, sun_pos: float, moon_pos: float) -> Tuple[THITHI, float]:
        if moon_pos < sun_pos:
            moon_pos += 360.
        delta = moon_pos - sun_pos
        idx = int(floor(delta / 12))
        return tithi_list[idx], delta

    def get_dates_on_day(self, day: VARAM, range: DateRange):
        rr = rrule.rrule(rrule.WEEKLY, byweekday=VARAM.get_relative_day(day),
                         dtstart=datetime.fromordinal(range.start.toordinal()), cache=True)
        print(rr.between(after=datetime.fromordinal(range.start.toordinal()),
                         before=datetime.fromordinal(range.stop.toordinal()), inc=True))

    def get_nakshaktra_on_day(self):
        # GoodDay.get_good_Tharai(Nakshatras.THIRUVATHIRAI)
        pass

    def _get_date_time_from_float(self, x: float):
        d_int = int(floor(x))
        d = date.fromordinal(d_int)
        t_float = x - d_int
        hr_float = t_float * 24
        hr = int(floor(hr_float))
        minute_float = (hr_float - hr) * 60
        minute = int(floor(minute_float))
        second = (minute_float - minute) * 60
        t = time(hr, minute, int(second))
        dt = datetime(d.year, d.month, d.day, hr, minute, int(second), tzinfo=timezone.utc)
        return dt

    def _get_date_time_as_float(self, x: datetime) -> float:
        return float(x.toordinal()) + (x.hour + (x.minute + x.second / 60) / 60) / 24

    def get_next_thiti(self, target: THITHI):
        low_deg = (int(target - 1)) * 12.
        target_deg = low_deg
        now = datetime.now(timezone.utc)
        self.compute(now)
        first_target = False
        if target_deg < (self.tithi - 1) * 12:
            orig_target = target_deg
            target_deg = 360
            first_target = True

        def calc(x: float):
            dt = self._get_date_time_from_float(x)
            self.compute(dt)
            # print(f'{x} --- curr_pos {self.pos_deg} target_deg {target_deg}')
            return (target_deg - self.pos_deg) ** 2

        if first_target:
            result = minimize_scalar(calc, bounds=(now.toordinal(), now.toordinal() + 30), method='Bounded',
                                     options={'disp': False})
            result = fmin(calc, result.x, ftol=1 / 24 / 60 / 60, disp=False)
            now = self._get_date_time_from_float(result[0])
            target_deg = orig_target

        result = minimize_scalar(calc, bounds=(float(now.toordinal()), float(now.toordinal()) + 29.7), method='Bounded',
                                 options={'xatol': 1 / 24 / 60 / 60})
        result = fmin(calc, result.x, disp=False)
        dt = self._get_date_time_from_float(result[0])

        self.compute(dt)
        while self.tithi != target:
            if self.tithi < target:
                delta = timedelta(minutes=1)
            else:
                delta = timedelta(minutes=-1)
            dt = dt + delta
            self.compute(dt)
        print(dt, dt.astimezone(timezone(timedelta(hours=5, minutes=30))), self.tithi, self.pos_deg,
              100 * (int(self.tithi) * 12 - self.pos_deg) / 12, '% left')

    def _compute_tithi_on(self, time: Time) -> Thithi:
        moon_pos = self._compute_moon_pos(time, False)
        sun_pos = self._compute_sun_pos(time, False)
        if len(time.shape) == 0:
            if moon_pos.degrees < sun_pos.degrees:
                moon_pos.degrees = moon_pos.degrees + 360.
        else:
            mask = moon_pos.degrees < sun_pos.degrees
            moon_pos.degrees[mask] = moon_pos.degrees[mask] + 360.
        delta = moon_pos.degrees - sun_pos.degrees
        return Thithi(Angle(degrees=delta))

    @staticmethod
    def _compute_ayanamsa(time: Time) -> Angle:
        utc_time = time.utc_datetime()

        def foo(utc_time_obj):
            a = -6.92416 + 16.90709 * (utc_time_obj.year * 1e-3) - 0.757371 * (utc_time_obj.year ** 2) * 1e-6
            b = ((utc_time_obj.month - 1) + utc_time_obj.day / 30) * 1.1574074 * 1e-3
            return a + b

        if len(time.shape) == 0:
            return Angle(degrees=foo(utc_time))
        else:
            vfunc = np.vectorize(foo)
            return Angle(degrees=vfunc(utc_time))

    def _compute_moon_pos(self, time: Time, correct_ayanamsa: bool = True) -> Angle:
        e = self.earth.at(time)
        _, moon_pos, _ = e.observe(self.moon).apparent().ecliptic_latlon('date')
        moon_pos.degrees = moon_pos.degrees % 360
        if correct_ayanamsa:
            ayanamsa = Almanac._compute_ayanamsa(time)
            moon_pos = Angle(degrees=(moon_pos.degrees - ayanamsa.degrees) % 360)
        return moon_pos

    def _compute_sun_pos(self, time: Time, correct_ayanamsa: bool = True) -> Angle:
        e = self.earth.at(time)
        _, sun_pos, _ = e.observe(self.sun).apparent().ecliptic_latlon('date')
        sun_pos.degrees = sun_pos.degrees % 360
        if correct_ayanamsa:
            ayanamsa = self._compute_ayanamsa(time)
            sun_pos.degrees = (sun_pos.degrees - ayanamsa.degrees) % 360
        return sun_pos

    def _compute_nakshaktram(self, time: Time) -> Nakshaktram:
        moon_pos = self._compute_moon_pos(time, True)
        return Nakshaktram(moon_pos)

    def _compute_month(self, time: Time) -> Month:
        sun_pos = self._compute_sun_pos(time, True)
        return Month(sun_pos)

    def _compute_sun_rise_sun_set(self, start: Time, stop: Time):
        times, rise_set = find_discrete(start, stop, sunrise_sunset(self.planets, self.topo))
        return times[rise_set], times[~rise_set]
