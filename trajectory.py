import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def lambert_v(GM, X1, X2, tof, string):
    R1 = X1[:3]  # X1（6次元）の最初の3つの要素を取り出す
    r1 = np.linalg.norm(R1)  # ノルムを計算する


    R2 = X2[:3]  # X2（6次元）の最初の3つの要素を取り出す
    r2 = np.linalg.norm(R2)  # ノルムを計算する


    cross12 = np.cross(R1, R2)  # R1とR2の外積を計算する
    theta = np.arccos(np.dot(R1, R2) / (r1*r2))  # 内積を計算し、アークコサインを取る


    if string == 'pro':
        if cross12[2] < 0:
            theta = 2 * np.pi - theta  # 第3象限または第4象限へ
    elif string == 'retro':
        if cross12[2] > 0:
            theta = 2 * np.pi - theta  # 第3象限または第4象限へ
   
    A = np.sin(theta) * np.sqrt(r1*r2/(1 - np.cos(theta)))
    
    ## 必要な関数の定義
    def C(z):
        if z == 0:
            return 1/2
        elif z > 0:
            return (1 - np.cos(np.sqrt(z))) / z
        elif z < 0:
            return (np.cosh(np.sqrt(-z)) - 1) / (-z)
    
    def S(z):
        if z == 0:
            return 1/6
        elif z > 0:
            return (np.sqrt(z) - np.sin(np.sqrt(z))) / z**(3/2)
        elif z < 0:
            return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (-z)**(3/2)

    def y(z):
        return r1 + r2 + A*(z*S(z) - 1)/np.sqrt(C(z))
    
    def F(z):
        return S(z)*(y(z)/C(z))**(3/2) + A*np.sqrt(y(z)) - np.sqrt(GM)*tof
    
    def dF(z):
        if z == 0:
            return np.sqrt(2)/40*y(0)**(3/2) + A/8*(np.sqrt(y(0)) + A*np.sqrt(1/(2*y(0))))
        else:
            return (y(z)/C(z))**(3/2) * ((C(z) - 3*S(z)/2/C(z))/2/z + 3*S(z)**2/4/C(z)) + A/8*(3*S(z)/C(z)*np.sqrt(y(z)) + A*np.sqrt(C(z)/y(z)))

    ## ニュートン法によるzの計算

    ## まずはy(z)が正の値になるまでzを増加させる
    z = 0  # 初期値
    #print(y(-100))
    # while y(z) < 0:
    #    z = z + 0.1
       #print(z, y(z))
    #print("initial value 0")
    #print(z, y(z))

    ## 次にF(z)が正の値になるまでzを増加させる
    while F(z) < 0:
       z = z + 0.1
       #print(z, F(z))
    #print("initial value 1")
    #print(z, F(z))

    # z = 0.1

    imax = 1000
    for i in range(imax):
        z = z - F(z)/dF(z)  # ニュートン法の更新式
        if y(z) < 0:
            print("y(z) is negative")
        elif abs(F(z)/dF(z)) < 1e-13:  # 収束条件
            break
        #else:
            #print(i, z, y(z))
            #print(abs(F(z)/dF(z)))        
    
    #print(i, z)
    #print(abs(F(z)/dF(z)))
    

    ## ラグランジュの係数の計算
    f = 1 - y(z)/r1
    g = A*np.sqrt(y(z)/GM)
    gdot = 1 - y(z)/r2

    ## 速度ベクトルの計算
    V1 = (1/g)*(R2 - f*R1)
    V2 = (1/g)*(gdot*R2 - R1)

    Vp1 = X1[3:6]  # X1（6次元）の4番目から6番目の要素を取り出す
    Vp2 = X2[3:6]  # X2（6次元）の4番目から6番目の要素を取り出す

    Vinf1 = V1 - Vp1  # V1からVp1を引く
    Vinf2 = V2 - Vp2  # V2からVp2を引く
    Vinf1_norm = np.linalg.norm(Vinf1)  # Vinf1のノルムを計算する  
    Vinf2_norm = np.linalg.norm(Vinf2)  # Vinf2のノルムを計算する
    C3 = Vinf1_norm**2
    DV_total = Vinf1_norm + Vinf2_norm

    # print("Vinf1:", Vinf1_norm)  # Vinf1を表示する
    # print("Vinf2:", Vinf2_norm)  # Vinf2を表示する

    return V1, V2, tof, Vinf1_norm, Vinf2_norm, C3, DV_total  # Vinf1とVinf2を返す

# 地球のエフェメリス情報を読み込み
ephemeris_E = []

with open("eph_E.txt", "r") as file:
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

with open("eph_M.txt", "r") as file:
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

# 地球のエフェメリス情報を読み込み
ephemeris_E_full = []

with open("eph_E_full.txt", "r") as file:
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
                        ephemeris_E_full.append(values)
            
            # 処理後、バッファをクリア
            buffer = []

# 地球のエフェメリス情報を読み込み
ephemeris_M_full = []

with open("eph_M_full.txt", "r") as file:
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
                        ephemeris_M_full.append(values)
            
            # 処理後、バッファをクリア
            buffer = []


# ephemeris_Eの2, 3, 4列目をx, y, z座標として抽出
x_Ekm = [row[1] for row in ephemeris_E]
y_Ekm = [row[2] for row in ephemeris_E]
z_Ekm = [row[3] for row in ephemeris_E]

x_Mkm = [row[1] for row in ephemeris_M]
y_Mkm = [row[2] for row in ephemeris_M]
z_Mkm = [row[3] for row in ephemeris_M]

x_E_fullkm = [row[1] for row in ephemeris_E_full]
y_E_fullkm = [row[2] for row in ephemeris_E_full]
z_E_fullkm = [row[3] for row in ephemeris_E_full]

x_M_fullkm = [row[1] for row in ephemeris_M_full]
y_M_fullkm = [row[2] for row in ephemeris_M_full]
z_M_fullkm = [row[3] for row in ephemeris_M_full]

AU = 1.495978707e8
x_E = [xi / AU for xi in x_Ekm]
y_E = [yi / AU for yi in y_Ekm]
z_E = [zi / AU for zi in z_Ekm]

x_M = [xi / AU for xi in x_Mkm]
y_M = [yi / AU for yi in y_Mkm]
z_M = [zi / AU for zi in z_Mkm]

x_E_full = [xi / AU for xi in x_E_fullkm]
y_E_full = [yi / AU for yi in y_E_fullkm]
z_E_full = [zi / AU for zi in z_E_fullkm]

x_M_full = [xi / AU for xi in x_M_fullkm]
y_M_full = [yi / AU for yi in y_M_fullkm]
z_M_full = [zi / AU for zi in z_M_fullkm]


i_dep = 0  # 地球の最初の行を取得
i_arr = len(ephemeris_M) - 1  # 火星の最終行を取得


GM_sun = 1.32712440041279419e11 # 太陽の重力定数 [km^3/s^2]
tof = (ephemeris_M[i_arr][0] - ephemeris_E[i_dep][0])*24*60*60  # [s]
X1 = np.array(ephemeris_E[i_dep][1:7])  # 地球の初期状態ベクトル
X2 = np.array(ephemeris_M[i_arr][1:7])  # 火星の初期状態ベクトル
string = 'pro'

l = lambert_v(GM_sun, X1, X2, tof, string)

# 太陽との二体問題の微分方程式
def two_body_problem(t, X, GM):
    r = np.linalg.norm(X[:3])  # 位置ベクトルのノルム
    dXdt = np.zeros(6)
    dXdt[:3] = X[3:]  # 速度ベクトル
    dXdt[3:] = -GM * X[:3] / r**3  # 加速度ベクトル
    return dXdt

# 初期値
X0 = np.zeros(6)
X0[:3] = X1[:3]  # 地球出発の位置ベクトル
X0[3:6] = l[0]   # 地球出発の速度ベクトル

# 時間範囲
t_span = (0, tof)  # 出発から到着まで
t_eval = np.linspace(0, tof, 1000)  # 時間の評価点

# 微分方程式を解く
sol = solve_ivp(two_body_problem, t_span, X0, args=(GM_sun,), t_eval=t_eval, rtol=1e-9, atol=1e-9)

# 解の位置ベクトルを取得
trajectory_x = sol.y[0] / AU
trajectory_y = sol.y[1] / AU
trajectory_z = sol.y[2] / AU


# 3次元プロット
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_E, y_E, z_E, color='lightblue', label='Earth Orbit')
ax.plot(x_M, y_M, z_M, color='orange', label='Mars Orbit')
ax.plot(x_E_full, y_E_full, z_E_full, color='lightblue', alpha=0.5, linestyle='--', label='Earth Orbit (Full)')
ax.plot(x_M_full, y_M_full, z_M_full, color='orange', alpha=0.5, linestyle='--', label='Mars Orbit (Full)')
# 軌道をプロット
ax.plot(trajectory_x, trajectory_y, trajectory_z, color='green', label='Transfer Orbit')
ax.scatter(trajectory_x[0], trajectory_y[0], trajectory_z[0], color='blue', label='Departure (Earth)')
ax.scatter(trajectory_x[-1], trajectory_y[-1], trajectory_z[-1], color='brown', label='Arrival (Mars)')
ax.scatter(0, 0, 0, color='red', label='Sun', s=200)  # 太陽の位置をプロット（2倍の大きさ）


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
plt.plot(x_E_full, y_E_full, color='lightblue', alpha=0.5, linestyle='--', label='Earth Orbit (Full)')
plt.plot(x_M_full, y_M_full, color='orange', alpha=0.5, linestyle='--', label='Mars Orbit (Full)')
plt.plot(trajectory_x, trajectory_y, color='green', label='Transfer Orbit (2D)')
plt.scatter(trajectory_x[0], trajectory_y[0], color='blue', label='Departure (Earth)')
plt.scatter(trajectory_x[-1], trajectory_y[-1], color='brown', label='Arrival (Mars)')
plt.scatter(0, 0, color='red', label='Sun', s=200)  # 太陽の位置をプロット（2倍の大きさ）
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.legend()
plt.grid()
plt.axis('equal')
plt.show()