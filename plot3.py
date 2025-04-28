import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.dates as mdates

# ファイルから結果を読み込み
data = np.loadtxt("result.txt")

# データを各軸に分解
dep_JD = data[:, 0]
arr_JD = data[:, 1]
tof_days = data[:, 2]
vinf1 = data[:, 3]
vinf2 = data[:, 4]
C3 = data[:, 5]

dep_unique = np.unique(dep_JD)
arr_unique = np.unique(arr_JD)

Dep_JD_grid, Arr_JD_grid = np.meshgrid(dep_unique, arr_unique)

C3_grid = np.full_like(Dep_JD_grid, np.nan, dtype=np.float64)

for i in range(len(data)):
    dep_idx = np.where(dep_unique == dep_JD[i])[0][0]
    arr_idx = np.where(arr_unique == arr_JD[i])[0][0]
    C3_grid[arr_idx, dep_idx] = C3[i]

# C3が90以上のところは無効にする
C3_grid[C3_grid > 90] = np.nan

# JDをdatetimeに変換
dep_dates = Time(dep_unique, format='jd').to_datetime()
arr_dates = Time(arr_unique, format='jd').to_datetime()

Dep_grid_dates, Arr_grid_dates = np.meshgrid(dep_dates, arr_dates)

# 等高線のレベルを指定
#levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]
levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

# プロット
plt.figure(figsize=(12, 9))

# 塗りつぶし等高線
cp = plt.contourf(Dep_grid_dates, Arr_grid_dates, C3_grid, levels=levels, cmap='jet', extend='both')

# # 等高線（線だけ）
# cs = plt.contour(Dep_grid_dates, Arr_grid_dates, C3_grid, levels=levels, colors='black', linewidths=0.8)

# # ラベルつける
# plt.clabel(cs, inline=True, fontsize=8, fmt="%.0f")

# カラーバー
plt.colorbar(cp, label="C3 (km²/s²)")

# 軸設定
plt.xlabel("Departure Date")
plt.ylabel("Arrival Date")
plt.title("Porkchop Plot: Earth to Mars (C3 <= 90)")

# 日付フォーマットにする
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.gca().yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

# x軸・y軸に自動で日付の間隔をいい感じに調整
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=2))
plt.gca().yaxis.set_major_locator(mdates.MonthLocator(interval=2))

# x軸の目盛りを斜めに回転
plt.xticks(rotation=90)

plt.grid(True)

# 保存
plt.savefig("porkchop_c3_below90_dates.png", dpi=300, bbox_inches='tight')

# 表示
plt.show(block=True)
