import numpy as np
import matplotlib.pyplot as plt

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

# 等高線のレベルを指定
levels_vinf1 = [2, 3, 4, 5, 6, 7, 8, 9, 10]
levels_vinf2 = [2, 3, 4, 5, 6, 7, 8, 9, 10]
levels_DV = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
levels_TOF = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

# 等高線の描画関数
def plot_contours(grid, levels, title, filename, colorbar_label, cmap='jet'):
    plt.figure(figsize=(10, 8))

    # 塗りつぶし等高線
    cp = plt.contourf(Dep_JD_grid, Arr_JD_grid, grid, levels=levels, cmap=cmap, extend='both')

    # 等高線（線だけ）
    cs = plt.contour(Dep_JD_grid, Arr_JD_grid, grid, levels=levels, colors='black', linewidths=0.8)
    TOFs = plt.contour(Dep_JD_grid, Arr_JD_grid, TOF_grid, levels=levels_TOF, colors='black', linewidths=0.5)

    # ラベルつける
    plt.clabel(cs, inline=True, fontsize=8, fmt="%.0f")
    plt.clabel(TOFs, inline=True, fontsize=8, fmt="%.0f")

    # カラーバー
    plt.colorbar(cp, label=colorbar_label)

    # 軸ラベルとか
    plt.xlabel("Departure JD")
    plt.ylabel("Arrival JD")
    plt.title(title)
    plt.grid(True)

    # 保存
    plt.savefig(filename, dpi=300, bbox_inches='tight')

    # 表示
    # plt.show(block=True)
    plt.show()

# 各等高線を描画して保存
plot_contours(vinf1_grid, levels_vinf1, "Porkchop Plot: vinf1", "porkchop_f3_vinf1.png", "vinf1 (km/s)")
plot_contours(vinf2_grid, levels_vinf2, "Porkchop Plot: vinf2", "porkchop_f3_vinf2.png", "vinf2 (km/s)")
plot_contours(DV_grid, levels_DV, "Porkchop Plot: DVtotal", "porkchop_f3_vtotal.png", "DVtotal (km/s)")
