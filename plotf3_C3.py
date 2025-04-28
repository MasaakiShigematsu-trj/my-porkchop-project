import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.dates as mdates

# ファイルから結果を読み込み
data = np.loadtxt("result_f3.txt")

# データを各軸に分解
dep_JD = data[:, 0]
arr_JD = data[:, 1]
tof_days = data[:, 2]
vinf1 = data[:, 3]
vinf2 = data[:, 4]
C3 = data[:, 5]
DVtotal = data[:, 6]

dep_unique = np.unique(dep_JD)
arr_unique = np.unique(arr_JD)

Dep_JD_grid, Arr_JD_grid = np.meshgrid(dep_unique, arr_unique)

C3_grid = np.full_like(Dep_JD_grid, np.nan, dtype=np.float64)
TOF_grid = np.full_like(Dep_JD_grid, np.nan, dtype=np.float64)

for i in range(len(data)):
    dep_idx = np.where(dep_unique == dep_JD[i])[0][0]
    arr_idx = np.where(arr_unique == arr_JD[i])[0][0]
    C3_grid[arr_idx, dep_idx] = C3[i]

for i in range(len(data)):
    dep_idx = np.where(dep_unique == dep_JD[i])[0][0]
    arr_idx = np.where(arr_unique == arr_JD[i])[0][0]
    TOF_grid[arr_idx, dep_idx] = tof_days[i]

# JDをdatetimeに変換
dep_dates = Time(dep_unique, format='jd').to_datetime()
arr_dates = Time(arr_unique, format='jd').to_datetime()

Dep_grid_dates, Arr_grid_dates = np.meshgrid(dep_dates, arr_dates)


print(np.isnan(TOF_grid).sum())  # NaNの数を確認

non_nan_count = np.sum(~np.isnan(TOF_grid))
print("NaNではない値の数:", non_nan_count)

# 等高線のレベルを指定
#levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
#levels2 = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
levels2 = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

# プロット
plt.figure(figsize=(10, 8))

# 塗りつぶし等高線
cp = plt.contourf(Dep_grid_dates, Arr_grid_dates, C3_grid, levels=levels, cmap='jet', extend='both')

# 等高線（線だけ）
cs = plt.contour(Dep_grid_dates, Arr_grid_dates, C3_grid, levels=levels, colors='black', linewidths=0.8)
TOFs = plt.contour(Dep_grid_dates, Arr_grid_dates, TOF_grid, levels=levels2, colors='black', linewidths=0.5)

# ラベルつける
plt.clabel(cs, inline=True, fontsize=8, fmt="%.0f")
plt.clabel(TOFs, inline=True, fontsize=8, fmt="%.0f")

# カラーバー
plt.colorbar(cp, label="C3 [km²/s²]")

# 軸ラベルとか
plt.xlabel("Departure Date")
plt.ylabel("Arrival Date")
plt.title("Porkchop Plot: Earth to Mars (C3)")

# 日付フォーマットにする
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.gca().yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

# x軸・y軸の日付の間隔を月単位で調整
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
plt.gca().yaxis.set_major_locator(mdates.MonthLocator(interval=1))

# x軸とy軸の範囲を調整（見たい範囲を指定）
plt.xlim([min(dep_dates), max(dep_dates)])  # 出発日の範囲
plt.ylim([min(arr_dates), max(arr_dates)])  # 到着日の範囲

# x軸の目盛りを斜めに回転
plt.xticks(rotation=80)


plt.grid(True)

# 保存（オプション）
plt.savefig("porkchop_f3_c3.png", dpi=300, bbox_inches='tight')

# 表示
plt.show(block=True)