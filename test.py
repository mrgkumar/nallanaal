from almanac import Almanac
from params import nakshaktra_ignore_list, avoid_nakshaktras_on_day, avoid_thithi_on_days, avoid_thitis_list
from params import avoid_table, marana_yoga_table
from params import NAKSHATRAS
import numpy as np

user_nakshatra = NAKSHATRAS.MAGAM
num_days = 90
mins_step = 30
a = Almanac()
t0 = a.time_scale.utc(2020, 4, 1)
t1 = a.time_scale.utc(2020, 4, 1, 0, 24*60*366)
month_times, months = a._get_month_start_end(t0, t1)
# sun_rise, sun_set = a._compute_sun_rise_sun_set(t0, t1)
tframe = a.time_scale.utc(2020, 4, 1, 0, range(0, 24 * 60 * num_days, mins_step))

a.time = tframe
a.compute_almanac()
mask = np.ones(tframe.shape, dtype=np.bool)
print(f'started with {mask.sum()}')
for t in avoid_thitis_list:
    mask[(t == a.thithi.name)] = False
print(f'after removing thithis {avoid_thitis_list} -- {mask.sum()}')
for n in nakshaktra_ignore_list:
    mask[n == a.nakshaktra.name] = False
print(f'after removing nakshaktras {nakshaktra_ignore_list} -- {mask.sum()}')
for day, thithis in avoid_thithi_on_days.items():
    cur_days = (a.vara == day)
    for t in thithis:
        mask[(t == a.thithi.name) & cur_days] = False
print(f'after removing thithi on days {avoid_thithi_on_days} -- {mask.sum()}')
for day, nakshaktras in avoid_nakshaktras_on_day.items():
    cur_days = (a.vara == day)
    for n in nakshaktras:
        mask[(n == a.nakshaktra.name) & cur_days] = False
print(f'after removing nakshaktras on days {avoid_nakshaktras_on_day} -- {mask.sum()}')
for day, nakshaktras in marana_yoga_table.items():
    cur_days = (a.vara == day)
    for n in nakshaktras:
        mask[(n == a.nakshaktra.name) & cur_days] = False
print(f'after removing marana yogas -- {mask.sum()}')

for day, (thithis, nakshaktra) in avoid_table.items():
    cur_days = (a.vara == day)
    for t in thithis:
        mask[(t == a.thithi.name) & cur_days & (a.nakshaktra.name == nakshaktra)] = False
print(f'after removing avoid table thihi, nakshaktra on days {avoid_table} -- {mask.sum()}')
func = lambda x: x in Almanac.tara_palan_dict[user_nakshatra]
vfunc = np.vectorize(func)
mask &= vfunc(a.nakshaktra.name)
print(f'after filtering for tara palan {user_nakshatra} -- {mask.sum()} ')
idxs = np.argwhere(mask)
for k in range(idxs.shape[0]):
    cur_id = idxs[k,0]
    print(a.time.utc_datetime()[cur_id].astimezone(), a.vara[cur_id], a.nakshaktra.name[cur_id], a.thithi.name[cur_id])

pass


