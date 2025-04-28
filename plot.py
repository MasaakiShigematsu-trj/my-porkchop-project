import numpy as np
import matplotlib.pyplot as plt

# ファイルから結果を読み込み
data = np.loadtxt("result.txt")

# データを各軸に分解
dep_JD = data[:, 0]  # 出発時刻 (JD)
arr_JD = data[:, 1]  # 到着時刻 (JD)
tof_days = data[:, 2]  # 飛行時間（日）
vinf1 = data[:, 3]  # 出発時の相対速度
vinf2 = data[:, 4]  # 到着時の相対速度
C3 = data[:, 5]      # 出発時のC3

# グリッド化する
# dep_JDとarr_JDのユニークな値を取得して並べる
dep_unique = np.unique(dep_JD)
arr_unique = np.unique(arr_JD)

# dep_unique, arr_uniqueを使ってグリッドを作る
Dep_JD_grid, Arr_JD_grid = np.meshgrid(dep_unique, arr_unique)

# C3をグリッドに割り当てる
C3_grid = np.full_like(Dep_JD_grid, np.nan, dtype=np.float64)

# 各点に対応するC3を格納
for i in range(len(data)):
    dep_idx = np.where(dep_unique == dep_JD[i])[0][0]
    arr_idx = np.where(arr_unique == arr_JD[i])[0][0]
    C3_grid[arr_idx, dep_idx] = C3[i]

# プロット
plt.figure(figsize=(10, 8))
cp = plt.contourf(Dep_JD_grid, Arr_JD_grid, C3_grid, levels=30)
plt.colorbar(cp, label="C3 (km²/s²)")
plt.xlabel("Departure JD")
plt.ylabel("Arrival JD")
plt.title("Porkchop Plot: Earth to Mars (C3 contour)")
plt.grid(True)
#plt.show()
plt.show(block=True)