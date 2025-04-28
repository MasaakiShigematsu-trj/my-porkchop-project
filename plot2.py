import numpy as np
import matplotlib.pyplot as plt

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

# 等高線のレベルを指定
# levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

# プロット
plt.figure(figsize=(10, 8))

# 塗りつぶし等高線
cp = plt.contourf(Dep_JD_grid, Arr_JD_grid, C3_grid, levels=levels, cmap='jet', extend='both')

# 等高線（線だけ）
cs = plt.contour(Dep_JD_grid, Arr_JD_grid, C3_grid, levels=levels, colors='black', linewidths=0.8)

# ラベルつける
plt.clabel(cs, inline=True, fontsize=8, fmt="%.0f")

# カラーバー
plt.colorbar(cp, label="C3 (km²/s²)")

# 軸ラベルとか
plt.xlabel("Departure JD")
plt.ylabel("Arrival JD")
plt.title("Porkchop Plot: Earth to Mars (C3=20,25,30,35,40)")
plt.grid(True)

# 保存（オプション）
plt.savefig("porkchop_c3_fixed_levels.png", dpi=300, bbox_inches='tight')

# 表示
plt.show(block=True)
