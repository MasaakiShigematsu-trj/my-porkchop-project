import numpy as np
from astropy.time import Time

# データを読み込む（nanを正しく読み取る）
data = np.loadtxt('result_f2.txt')

# 3列目（インデックス2）のデータと7列目（インデックス6）のデータを取り出す
col3 = data[:, 2]
col7 = data[:, 6]

# 3列目の値が250以下かつ、7列目の値が最小となる行を探す
valid_rows = (col3 <= 250) & (~np.isnan(col7))  # 3列目が250以下かつ7列目がNaNでない行

# 条件を満たす行から7列目の最小値を探す
valid_data = data[valid_rows]  # 条件を満たすデータを抽出
col7_valid = valid_data[:, 6]  # 7列目のデータ

# 最小値のインデックスを探す
min_index_valid = np.argmin(col7_valid)

# 最小値を持つ行を取得
min_row = valid_data[min_index_valid]

print("3列目が250以下で7列目が最小の行：")
print(min_row)

# ユリウス日を西暦に変換
julian_date = min_row[0]
gregorian_date1 = Time(julian_date, format='jd').iso

print("出発日：")
print(gregorian_date1)

# ユリウス日を西暦に変換
julian_date = min_row[1]
gregorian_date2 = Time(julian_date, format='jd').iso

print("到着日：")
print(gregorian_date2)