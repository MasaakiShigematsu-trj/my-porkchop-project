import numpy as np
from astropy.time import Time

# データを読み込む（nanを正しく読み取る）
data = np.loadtxt('result_f2.txt')

# 7列目（インデックス6）のデータを取り出す
col7 = data[:, 6]

# NaNを除外して最小値を探す
valid_indices = ~np.isnan(col7)  # NaNでない場所だけTrue
col7_valid = col7[valid_indices]

# 最小値のインデックスを探す
min_index_valid = np.argmin(col7_valid)

# 元のdataに対応するインデックスに直す
valid_rows = np.where(valid_indices)[0]
min_row_index = valid_rows[min_index_valid]

# 最小値を持つ行を取得
min_row = data[min_row_index]

print("7列目が最小の行：")
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