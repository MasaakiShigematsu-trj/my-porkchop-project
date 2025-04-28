import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pykep as pk
from pykep.planet import jpl_lp
from astropy.time import Time

# ----------------------
# 1. 出発日・到着日（ユリウス日）
# ----------------------
# 出発日と到着日を定義
# 出発日・到着日
launch_epoch = pk.epoch_from_string('2030-12-16 00:00:00')
arrival_epoch = pk.epoch_from_string('2031-09-23 00:00:00')

TOF_sec = (arrival_epoch.mjd2000 - launch_epoch.mjd2000) * 86400  # [秒]



# ----------------------
# 2. 惑星データの取得
# ----------------------
earth = jpl_lp('earth')
mars = jpl_lp('mars')

r1, v1 = earth.eph(launch_epoch)
r2, v2 = mars.eph(arrival_epoch)
mu_sun = pk.MU_SUN  # [km^3/s^2]

# ----------------------
# 3. ランベール問題を解く
# ----------------------
lambert_solutions = pk.lambert_problem(r1, r2, TOF_sec, mu_sun)
v1_departure = lambert_solutions.get_v1()[0]
v2_arrival = lambert_solutions.get_v2()[0]

# ----------------------
# 4. 軌道を積分する（宇宙機）
# ----------------------
num_points = 1000
times = np.linspace(0, TOF_sec, num_points)
r_sc = []

for t in times:
    M = np.linalg.norm(v1_departure)**2/2 - mu_sun/np.linalg.norm(r1)
    a = -mu_sun/(2*M)
    r, _ = pk.propagate_lagrangian(r1, v1_departure, t, mu_sun)
    r_sc.append(r)

r_sc = np.array(r_sc)

# ----------------------
# 5. 地球と火星の軌道を描画
# ----------------------
# 周回軌道をサンプリング
theta = np.linspace(0, 2*np.pi, num_points)

# 地球軌道（円軌道近似）
r_earth_orbit = []
for angle in theta:
    r, _ = earth.eph(launch_epoch + (angle/(2*np.pi))*365.25)
    r_earth_orbit.append(r)
r_earth_orbit = np.array(r_earth_orbit)

# 火星軌道（円軌道近似）
r_mars_orbit = []
for angle in theta:
    r, _ = mars.eph(launch_epoch + (angle/(2*np.pi))*686.98)
    r_mars_orbit.append(r)
r_mars_orbit = np.array(r_mars_orbit)

# ----------------------
# 6. 3Dプロット
# ----------------------
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 太陽描画
ax.scatter(0, 0, 0, color='yellow', s=100, label='Sun')

# 地球軌道
ax.plot(r_earth_orbit[:,0], r_earth_orbit[:,1], r_earth_orbit[:,2], 'b', label='Earth Orbit')
# 火星軌道
ax.plot(r_mars_orbit[:,0], r_mars_orbit[:,1], r_mars_orbit[:,2], 'r', label='Mars Orbit')
# 宇宙機遷移軌道
ax.plot(r_sc[:,0], r_sc[:,1], r_sc[:,2], 'g', label='Spacecraft Trajectory')

# 地球出発点
ax.scatter(r1[0], r1[1], r1[2], color='blue', s=50, label='Earth Departure')
# 火星到着点
ax.scatter(r2[0], r2[1], r2[2], color='red', s=50, label='Mars Arrival')

ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
ax.set_zlabel('Z [km]')
ax.set_title('Earth-Mars Transfer (Lambert Solution)')
ax.legend()
ax.grid(True)
ax.set_box_aspect([1,1,1])  # aspect ratioを1:1:1に

# 保存（オプション）
plt.savefig("trj1.png", dpi=300, bbox_inches='tight')

# 表示
plt.show(block=True)
