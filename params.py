from enum import IntEnum

from dateutil import relativedelta as relativedelta
from datetime import datetime, time, date, timezone, timedelta
from skyfield.units import Angle
from skyfield.timelib import Time
from math import floor
import numpy as np


class NAKSHATRAS(IntEnum):
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


class Nakshaktram:
    _list = [n for n in NAKSHATRAS]
    _span_deg = 13 + 20 / 60

    def __init__(self, moon_lon: Angle):
        """
        :param moon_lon: ayanamsa corrected
        """
        self.angle = moon_lon
        if type(moon_lon.degrees) is np.float64:
            self.name = Nakshaktram._list[int(floor(moon_lon.degrees / Nakshaktram._span_deg))]
        elif type(moon_lon.degrees) is np.ndarray:
            self.name = np.asarray(
                [Nakshaktram._list[int(x)] for x in np.floor(moon_lon.degrees / Nakshaktram._span_deg)], NAKSHATRAS)
        else:
            raise ValueError(f"unknown data type for moon_lon.degrees {type(moon_lon.degrees)}")
        self.rasi = Rasi(moon_lon)


class VARAM(IntEnum):
    __order__ = "SUN MON TUE WED THU FRI SAT"
    SUN = relativedelta.SU.weekday
    MON = relativedelta.MO.weekday
    TUE = relativedelta.TU.weekday
    WED = relativedelta.WE.weekday
    THU = relativedelta.TH.weekday
    FRI = relativedelta.FR.weekday
    SAT = relativedelta.SA.weekday

    @staticmethod
    def get_relative_day(day):
        table = {VARAM.SUN: relativedelta.SU,
                 VARAM.MON: relativedelta.MO,
                 VARAM.TUE: relativedelta.TU,
                 VARAM.WED: relativedelta.WE,
                 VARAM.THU: relativedelta.TH,
                 VARAM.FRI: relativedelta.FR,
                 VARAM.SAT: relativedelta.SA
                 }
        return table[day]

class Varam:
    _dict = {int(v): v for v in VARAM}

class RASI(IntEnum):
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


class Rasi:
    _list = [r for r in RASI]
    _span_deg = 30.0

    def __init__(self, moon_pos: Angle):
        if type(moon_pos.degrees) is np.float64:
            self.name = Rasi._list[int(floor(moon_pos.degrees / Rasi._span_deg))]
        elif type(moon_pos.degrees) is np.ndarray:
            self.name = np.array([Rasi._list[int(r)] for r in np.floor(moon_pos.degrees / Rasi._span_deg)], RASI)
        else:
            raise ValueError(f"unknown data type for moon_pos.degrees {type(moon_pos.degrees)}")


class THITHI(IntEnum):
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


class Thithi:
    _list = [t for t in THITHI]
    _span = 360 / len(_list)

    def __init__(self, angle: Angle):
        self.angle = angle
        self.angle.degrees = angle.degrees % 360
        if type(angle.degrees) is np.float64:
            self.name = Thithi._list[int(floor(angle.degrees / 12))]
        elif type(angle.degrees) is np.ndarray:
            self.name = np.array([Thithi._list[int(t)] for t in np.floor(angle.degrees / 12.)], THITHI)
        else:
            ValueError(f"Unknown type for angle.degrees {type(angle.degrees)}")


class MONTHS(IntEnum):
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


class Month:
    _list = [m for m in MONTHS]
    _span_deg = 30.0

    def __init__(self, sun_pos: Angle):
        self.name = Month._list[int(floor(sun_pos.degrees / Month._span_deg))]


class AMRUTHATHI_YOGAMS(IntEnum):
    SIDHA = 1
    AMRUTHA = 2
    MARANA = 3


class DIRECTIONS(IntEnum):
    NORTH = 1
    EAST = 2
    SOUTH = 3
    WEST = 4
    NORTH_EAST = 5
    SOUTH_EAST = 6
    NORTH_WEST = 7
    SOUTH_WEST = 8


class TimeRange:
    def __init__(self, start_hr, start_min, stop_hr, stop_min):
        self.start = time(start_hr, start_min)
        self.stop = time(stop_hr, stop_min)


class DateRange:
    def __init__(self, start_date: date, stop_date: date):
        self.start = start_date
        self.stop = stop_date


marana_yoga_table = {
    VARAM.SUN: [NAKSHATRAS.MAGAM, NAKSHATRAS.VISAKAM, NAKSHATRAS.ANUSHAM, NAKSHATRAS.KETAI, NAKSHATRAS.AVITAM],
    VARAM.MON: [NAKSHATRAS.KRITTIKAI, NAKSHATRAS.MAGAM, NAKSHATRAS.VISAKAM, NAKSHATRAS.UTHTHARADAM,
                NAKSHATRAS.POORATTATHI],
    VARAM.TUE: [NAKSHATRAS.THIRUVATHIRAI, NAKSHATRAS.VISAKAM, NAKSHATRAS.SATHAYAM, NAKSHATRAS.POORATTATHI],
    VARAM.WED: [NAKSHATRAS.ASHWINI, NAKSHATRAS.HASTHAM, NAKSHATRAS.MOOLAM, NAKSHATRAS.AVITAM, NAKSHATRAS.REVATHI],
    VARAM.THU: [NAKSHATRAS.KRITTIKAI, NAKSHATRAS.ROHINI, NAKSHATRAS.MIRGASIRISHAM, NAKSHATRAS.THIRUVATHIRAI,
                NAKSHATRAS.UTHIRAM, NAKSHATRAS.SATHAYAM],
    VARAM.FRI: [NAKSHATRAS.ROHINI, NAKSHATRAS.POOSAM, NAKSHATRAS.AAYILYAM, NAKSHATRAS.MAGAM, NAKSHATRAS.KETAI,
                NAKSHATRAS.THIRUVOONAM],
    VARAM.SAT: [NAKSHATRAS.AAYILYAM, NAKSHATRAS.UTHIRAM, NAKSHATRAS.HASTHAM, NAKSHATRAS.CHITIRAI,
                NAKSHATRAS.POORATTATHI, NAKSHATRAS.REVATHI]
}

soolai_table = {
    VARAM.SUN: [DIRECTIONS.WEST, DIRECTIONS.NORTH_WEST],
    VARAM.MON: [DIRECTIONS.EAST, DIRECTIONS.SOUTH_WEST],
    VARAM.TUE: [DIRECTIONS.NORTH, DIRECTIONS.NORTH_WEST],
    VARAM.WED: [DIRECTIONS.NORTH, DIRECTIONS.NORTH_EAST],
    VARAM.THU: [DIRECTIONS.SOUTH, DIRECTIONS.SOUTH_EAST],
    VARAM.FRI: [DIRECTIONS.WEST, DIRECTIONS.SOUTH_WEST],
    VARAM.SAT: [DIRECTIONS.EAST, DIRECTIONS.SOUTH_WEST]
}

kari_naal = {
    MONTHS.MESHAM: [6, 15],
    MONTHS.RISHABAM: [7, 16, 17],
    MONTHS.MITHUNAM: [1, 6],
    MONTHS.KADAKAM: [2, 10, 20],
    MONTHS.SIMHAM: [2, 9, 28],
    MONTHS.KANNI: [16, 29],
    MONTHS.THULAM: [6, 20],
    MONTHS.VRICHIKAM: [1, 10, 17],
    MONTHS.DHANUSH: [6, 9, 11],
    MONTHS.MAKARAM: [1, 2, 3, 11, 17],
    MONTHS.KUMBAM: [15, 16, 17],
    MONTHS.MEENAM: [6, 15, 19]
}

dhaniya_naal = {
    MONTHS.MESHAM: [3, 20],
    MONTHS.RISHABAM: [9, 22],
    MONTHS.MITHUNAM: [8, 22],
    MONTHS.KADAKAM: [7, 20],
    MONTHS.SIMHAM: [7, 18],
    MONTHS.KANNI: [9, 26],
    MONTHS.THULAM: [8, 10],
    MONTHS.VRICHIKAM: [8, 14],
    MONTHS.DHANUSH: [8, 26],
    MONTHS.MAKARAM: [8, 15],
    MONTHS.KUMBAM: [15, 24],
    MONTHS.MEENAM: [18, 24]
}

###                 set id, curr stars,                                             compatible_set
tara_palan_table = [(1, [NAKSHATRAS.ASHWINI, NAKSHATRAS.MAGAM, NAKSHATRAS.MOOLAM], [2, 4, 6, 8, 9]),
                    (2, [NAKSHATRAS.BHARANI, NAKSHATRAS.POORAM, NAKSHATRAS.POORADAM], [1, 3, 5, 7, 9]),
                    (3, [NAKSHATRAS.KRITTIKAI, NAKSHATRAS.UTHIRAM, NAKSHATRAS.UTHTHARADAM], [1, 2, 4, 6, 8]),
                    (4, [NAKSHATRAS.ROHINI, NAKSHATRAS.HASTHAM, NAKSHATRAS.THIRUVOONAM], [2, 3, 5, 7, 9]),
                    (5, [NAKSHATRAS.MIRGASIRISHAM, NAKSHATRAS.CHITIRAI, NAKSHATRAS.AVITAM], [1, 3, 4, 6, 8]),
                    (6, [NAKSHATRAS.THIRUVATHIRAI, NAKSHATRAS.SWATHI, NAKSHATRAS.SATHAYAM], [2, 4, 5, 7, 9]),
                    (7, [NAKSHATRAS.PUNARPOOSAM, NAKSHATRAS.VISAKAM, NAKSHATRAS.POORATTATHI], [1, 3, 5, 6, 8]),
                    (8, [NAKSHATRAS.POOSAM, NAKSHATRAS.ANUSHAM, NAKSHATRAS.UTHTHRATTATHI], [2, 4, 6, 7, 9]),
                    (9, [NAKSHATRAS.AAYILYAM, NAKSHATRAS.KETAI, NAKSHATRAS.REVATHI], [1, 3, 5, 7, 8])]

nakshaktra_ignore_list = [NAKSHATRAS.BHARANI, NAKSHATRAS.KRITTIKAI, NAKSHATRAS.THIRUVATHIRAI, NAKSHATRAS.AAYILYAM,
                          NAKSHATRAS.POORADAM, NAKSHATRAS.POORAM, NAKSHATRAS.POORATTATHI, NAKSHATRAS.KETAI,
                          NAKSHATRAS.VISAKAM, NAKSHATRAS.CHITIRAI, NAKSHATRAS.SWATHI, NAKSHATRAS.MAGAM]
good_horai = {VARAM.SUN: [TimeRange(7, 30, 10, 0), TimeRange(14, 0, 16, 30), TimeRange(21, 0, 23, 59)],
              VARAM.MON: [TimeRange(6, 0, 7, 0), TimeRange(12, 0, 14, 0), TimeRange(18, 0, 21, 0),
                          TimeRange(22, 0, 23, 0)],
              VARAM.TUE: [TimeRange(10, 30, 11, 0), TimeRange(12, 0, 13, 0), TimeRange(16, 30, 18, 0),
                          TimeRange(19, 0, 20, 0)],
              VARAM.WED: [TimeRange(9, 0, 10, 0), TimeRange(13, 30, 15, 0), TimeRange(16, 0, 17, 0),
                          TimeRange(21, 0, 22, 0), TimeRange(23, 0, 23, 59)],
              VARAM.THU: [TimeRange(9, 0, 10, 30), TimeRange(13, 0, 13, 30), TimeRange(16, 30, 18, 0),
                          TimeRange(18, 0, 19, 0), TimeRange(20, 0, 21, 0)],
              VARAM.FRI: [TimeRange(6, 0, 9, 0), TimeRange(13, 0, 13, 30), TimeRange(17, 0, 18, 0),
                          TimeRange(20, 0, 21, 0), TimeRange(22, 30, 23, 0)],
              VARAM.SAT: [TimeRange(7, 0, 7, 30), TimeRange(10, 30, 13, 0), TimeRange(17, 0, 19, 30),
                          TimeRange(21, 0, 22, 0)]}

# avoid_thithis = [THITHI.KRI_PRATHAMAI, THITHI.SUK_PRATHAMAI,
#                  THITHI.SUK_NAVAMI, THITHI.KRI_NAVAMI,
#                  THITHI.SUK_ASHTAMI, THITHI.KRI_ASHTAMI]

avoid_nakshaktras_on_day = {VARAM.SUN: [NAKSHATRAS.BHARANI, NAKSHATRAS.KRITTIKAI, NAKSHATRAS.MIRGASIRISHAM,
                                        NAKSHATRAS.MAGAM, NAKSHATRAS.VISAKAM, NAKSHATRAS.ANUSHAM, NAKSHATRAS.KETAI,
                                        NAKSHATRAS.POORATTATHI],
                            VARAM.MON: [NAKSHATRAS.CHITIRAI, NAKSHATRAS.KRITTIKAI, NAKSHATRAS.MAGAM, NAKSHATRAS.VISAKAM,
                                        NAKSHATRAS.ANUSHAM, NAKSHATRAS.POORAM, NAKSHATRAS.POORATTATHI],
                            VARAM.TUE: [NAKSHATRAS.UTHTHARADAM, NAKSHATRAS.THIRUVATHIRAI, NAKSHATRAS.KETAI,
                                        NAKSHATRAS.THIRUVOONAM, NAKSHATRAS.AVITAM, NAKSHATRAS.SATHAYAM],
                            VARAM.WED: [NAKSHATRAS.AVITAM, NAKSHATRAS.ASHWINI, NAKSHATRAS.BHARANI, NAKSHATRAS.KRITTIKAI,
                                        NAKSHATRAS.MOOLAM, NAKSHATRAS.THIRUVOONAM],
                            VARAM.THU: [NAKSHATRAS.KETAI, NAKSHATRAS.MIRGASIRISHAM, NAKSHATRAS.PUNARPOOSAM,
                                        NAKSHATRAS.POOSAM, NAKSHATRAS.POORADAM, NAKSHATRAS.REVATHI],
                            VARAM.FRI: [NAKSHATRAS.POORADAM, NAKSHATRAS.ROHINI, NAKSHATRAS.MIRGASIRISHAM,
                                        NAKSHATRAS.POOSAM, NAKSHATRAS.VISAKAM, NAKSHATRAS.HASTHAM, NAKSHATRAS.ANUSHAM,
                                        NAKSHATRAS.AVITAM],
                            VARAM.SAT: [NAKSHATRAS.REVATHI, NAKSHATRAS.PUNARPOOSAM, NAKSHATRAS.POOSAM,
                                        NAKSHATRAS.UTHIRAM, NAKSHATRAS.HASTHAM]
                            }

avoid_table = {
    VARAM.SUN: ([THITHI.KRI_PANCHAMI, THITHI.SUK_PANCHAMI], NAKSHATRAS.HASTHAM),
    VARAM.MON: ([THITHI.SUK_SAPTHAMI, THITHI.KRI_SAPTHAMI], NAKSHATRAS.THIRUVOONAM),
    VARAM.TUE: ([THITHI.SUK_SAPTHAMI, THITHI.KRI_SAPTHAMI], NAKSHATRAS.ASHWINI),
    VARAM.WED: ([THITHI.KRI_ASHTAMI, THITHI.SUK_ASHTAMI], NAKSHATRAS.ANUSHAM),
    VARAM.THU: ([THITHI.KRI_THRITHIYAI, THITHI.SUK_THRITHIYAI], NAKSHATRAS.POOSAM),
    VARAM.FRI: ([THITHI.KRI_NAVAMI, THITHI.SUK_NAVAMI], NAKSHATRAS.REVATHI),
    VARAM.SAT: ([THITHI.KRI_YEKADASI, THITHI.SUK_YEKADASI], NAKSHATRAS.REVATHI)
}

good_thithi_on_days = {VARAM.SUN: [THITHI.KRI_ASHTAMI, THITHI.SUK_ASHTAMI],
                       VARAM.MON: [THITHI.KRI_NAVAMI, THITHI.SUK_NAVAMI],
                       VARAM.TUE: [THITHI.SUK_SHASTI, THITHI.KRI_SHASTI],
                       VARAM.WED: [THITHI.KRI_THRITHIYAI, THITHI.SUK_THRITHIYAI],
                       VARAM.THU: [THITHI.KRI_YEKADASI, THITHI.SUK_YEKADASI],
                       VARAM.FRI: [THITHI.KRI_THRIYODASI, THITHI.SUK_THRIYODASI],
                       VARAM.SAT: [THITHI.KRI_CHATHURDASI, THITHI.SUK_CHATHURDASI]}

avoid_thithi_on_days = {VARAM.SUN: [THITHI.KRI_CHATHURDASI, THITHI.SUK_CHATHURDASI],
                        VARAM.TUE: [THITHI.SUK_SAPTHAMI, THITHI.KRI_SAPTHAMI],
                        VARAM.MON: [THITHI.SUK_SHASTI, THITHI.KRI_SHASTI],
                        VARAM.WED: [THITHI.KRI_DVITHIYAI, THITHI.SUK_DVITHIYAI],
                        VARAM.THU: [THITHI.KRI_ASHTAMI, THITHI.SUK_ASHTAMI],
                        VARAM.FRI: [THITHI.KRI_NAVAMI, THITHI.SUK_NAVAMI],
                        VARAM.SAT: [THITHI.KRI_SAPTHAMI, THITHI.SUK_SAPTHAMI]}

good_thithi_list = [THITHI.SUK_DASAMI, THITHI.SUK_YEKADASI, THITHI.SUK_DVADASI, THITHI.SUK_THRIYODASI,
                    THITHI.KRI_DVITHIYAI, THITHI.KRI_THRITHIYAI, THITHI.KRI_CHATHURTHI, THITHI.KRI_PANCHAMI,
                    THITHI.KRI_SHASTI, THITHI.KRI_SAPTHAMI]

avoid_thitis_list = [THITHI.AMMAVASAI, THITHI.POURNAMI, THITHI.SUK_ASHTAMI, THITHI.KRI_ASHTAMI, THITHI.SUK_NAVAMI,
                     THITHI.KRI_NAVAMI, THITHI.SUK_CHATHURDASI, THITHI.KRI_CHATHURDASI, THITHI.SUK_PRATHAMAI,
                     THITHI.KRI_PRATHAMAI]
