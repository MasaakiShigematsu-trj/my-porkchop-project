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

# 出発日と到着日を一意の値に
dep_unique = np.unique(dep_JD)
arr_unique = np.unique(arr_JD)

# メッシュグリッドを作成
Dep_JD_grid, Arr_JD_grid = np.meshgrid(dep_unique, arr_unique)

# グリッドを初期化
vinf1_grid = np.full_like(Dep_JD_grid, np.nan, dtype=np.float64)
vinf2_grid = np.full_like(Dep_JD_grid, np.nan, dtype=np.float64)
DV_grid = np.full_like(Dep_JD_grid, np.nan, dtype=np.float64)
TOF_grid = np.full_like(Dep_JD_grid, np.nan, dtype=np.float64)

# 各データを対応するグリッドに格納
for i in range(len(data)):
    dep_idx = np.where(dep_unique == dep_JD[i])[0][0]
    arr_idx = np.where(arr_unique == arr_JD[i])[0][0]
    vinf1_grid[arr_idx, dep_idx] = vinf1[i]
    vinf2_grid[arr_idx, dep_idx] = vinf2[i]
    DV_grid[arr_idx, dep_idx] = DVtotal[i]
    TOF_grid[arr_idx, dep_idx] = tof_days[i]

# JDをdatetimeに変換
dep_dates = Time(dep_unique, format='jd').to_datetime()
arr_dates = Time(arr_unique, format='jd').to_datetime()

Dep_grid_dates, Arr_grid_dates = np.meshgrid(dep_dates, arr_dates)

# 等高線のレベルを指定
levels_vinf1 = [2, 3, 4, 5, 6, 7, 8, 9, 10]
levels_vinf2 = [2, 3, 4, 5, 6, 7, 8, 9, 10]
levels_DV = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
levels_TOF = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

# 等高線の描画関数
def plot_contours(grid, levels, title, filename, colorbar_label, cmap='jet'):
    plt.figure(figsize=(10, 8))

    # 塗りつぶし等高線
    cp = plt.contourf(Dep_grid_dates, Arr_grid_dates, grid, levels=levels, cmap=cmap, extend='both')

    # 等高線（線だけ）
    cs = plt.contour(Dep_grid_dates, Arr_grid_dates, grid, levels=levels, colors='black', linewidths=0.8)
    TOFs = plt.contour(Dep_grid_dates, Arr_grid_dates, TOF_grid, levels=levels_TOF, colors='black', linewidths=0.5)

    # ラベルつける
    plt.clabel(cs, inline=True, fontsize=8, fmt="%.0f")
    plt.clabel(TOFs, inline=True, fontsize=8, fmt="%.0f")

    # カラーバー
    plt.colorbar(cp, label=colorbar_label)

    # 軸ラベルとか
    plt.xlabel("Departure Date")
    plt.ylabel("Arrival Date")
    plt.title(title)

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

    # 保存
    plt.savefig(filename, dpi=300, bbox_inches='tight')

    # 表示
    # plt.show(block=True)
    plt.show()

# 各等高線を描画して保存
plot_contours(vinf1_grid, levels_vinf1, "Porkchop Plot : Vinf_0", "porkchop_f3_vinf1.png", "Vinf_0 [km/s]")
plot_contours(vinf2_grid, levels_vinf2, "Porkchop Plot : Vinf_f", "porkchop_f3_vinf2.png", "Vinf_f [km/s]")
plot_contours(DV_grid, levels_DV, "Porkchop Plot : V_total ( = Vinf_0 + Vinf_f)", "porkchop_f3_vtotal.png", "V_total [km/s]")
