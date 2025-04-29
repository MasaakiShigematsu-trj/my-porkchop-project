import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 地球のエフェメリス情報を読み込み
ephemeris_E = []

with open("family2_E.txt", "r") as file:
    buffer = []  # 一時的に行を保持するリスト
    for line in file:
        # 行をバッファに追加
        buffer.append(line.strip())

        # バッファに4行たまったら処理を開始
        if len(buffer) == 4:
            for i, buf_line in enumerate(buffer):
                if "= A.D." in buf_line:  # `= A.D.`を含む場合
                    # `=A.D.` の前の値を取得
                    JD_value = buf_line.split()[0]
                    
                    # 次の行の3つの値を取得
                    r_line = buffer[i + 1].split() if i + 1 < len(buffer) else []
                    v_line = buffer[i + 2].split() if i + 2 < len(buffer) else []
                    
                    if len(r_line) >= 3 and len(v_line) >= 3:
                        values = [float(JD_value)] + list(map(float, r_line[:3])) + list(map(float, v_line[:3]))
                        ephemeris_E.append(values)
            
            # 処理後、バッファをクリア
            buffer = []


# 火星のエフェメリス情報を読み込み
ephemeris_M = []

with open("family2_M.txt", "r") as file:
    buffer = []  # 一時的に行を保持するリスト
    for line in file:
        # 行をバッファに追加
        buffer.append(line.strip())

        # バッファに4行たまったら処理を開始
        if len(buffer) == 4:
            for i, buf_line in enumerate(buffer):
                if "= A.D." in buf_line:  # `= A.D.`を含む場合
                    # `=A.D.` の前の値を取得
                    JD_value = buf_line.split()[0]
                    
                    # 次の行の3つの値を取得
                    r_line = buffer[i + 1].split() if i + 1 < len(buffer) else []
                    v_line = buffer[i + 2].split() if i + 2 < len(buffer) else []
                    
                    if len(r_line) >= 3 and len(v_line) >= 3:
                        values = [float(JD_value)] + list(map(float, r_line[:3])) + list(map(float, v_line[:3]))
                        ephemeris_M.append(values)
                        
            
            # 処理後、バッファをクリア
            buffer = []


# ephemeris_Eの2, 3, 4列目をx, y, z座標として抽出
x_Ekm = [row[1] for row in ephemeris_E]
y_Ekm = [row[2] for row in ephemeris_E]
z_Ekm = [row[3] for row in ephemeris_E]

x_Mkm = [row[1] for row in ephemeris_M]
y_Mkm = [row[2] for row in ephemeris_M]
z_Mkm = [row[3] for row in ephemeris_M]

AU = 1.495978707e8
x_E = [xi / AU for xi in x_Ekm]
y_E = [yi / AU for yi in y_Ekm]
z_E = [zi / AU for zi in z_Ekm]

x_M = [xi / AU for xi in x_Mkm]
y_M = [yi / AU for yi in y_Mkm]
z_M = [zi / AU for zi in z_Mkm]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_E, y_E, z_E, color='lightblue', label='Earth Orbit')
ax.plot(x_M, y_M, z_M, color='orange', label='Mars Orbit')

ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')

# 軸の表示範囲を調整
ax.set_xlim(-2.0, 2.0)
ax.set_ylim(-2.0, 2.0)
ax.set_zlim(-0.5, 0.5)

# 立方体の枠で表示
ax.set_box_aspect([1, 1, 1])

# ★ Z軸の目盛り間隔を0.25 AUに設定 ★
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.zaxis.set_major_locator(ticker.MultipleLocator(0.5))

ax.legend()
plt.tight_layout()
plt.show()


# 2次元プロット
plt.figure()
plt.plot(x_E, y_E, color='lightblue', label='Earth Orbit (2D)')
plt.plot(x_M, y_M, color='orange', label='Mars Orbit (2D)')
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.legend()
plt.grid()
plt.axis('equal')
plt.show()