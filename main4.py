## [2](1)専用

#######################################
#  普遍変数を用いてランベルト問題を解く  #
#######################################
import pykep as pk
import numpy as np
from numpy import array

## 入力引数
## GM : 中心天体の重力定数
## X1, X2 : 出発天体と目標天体の６次元状態ベクトル [km, km/s]
## tof : 航行時間 [s]
## string : 'pro'ならば'prograde', 'retro'ならば'retrograde'


def lambert(GM, X1, X2, tof, string):
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

    return Vinf1_norm, Vinf2_norm, C3, DV_total  # Vinf1とVinf2を返す

GM_earth = 398600.4354360959 # 地球の重力定数 [km^3/s^2]
GM_sun = 1.32712440041279419e11 # 太陽の重力定数 [km^3/s^2]

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


print(len(ephemeris_E))
print(len(ephemeris_M))



result = []
for i_dep in range(1): #(len(ephemeris_E)):
  print("hello3")
  print(i_dep)
  for i_arr in range(len(ephemeris_M)):
    tof = (ephemeris_M[i_arr][0] - ephemeris_E[i_dep][0])*24*60*60  # [s]
    if tof > 50*24*60*60:
      try:
        X1 = np.array(ephemeris_E[i_dep][1:7])  # 地球の初期状態ベクトル
        X2 = np.array(ephemeris_M[i_arr][1:7])  # 火星の初期状態ベクトル
        string = 'pro'
        l = lambert(GM_sun, X1, X2, tof, string)
        #l = pk.lambert_problem(ephemeris_E[i_dep][1:4], ephemeris_M[i_arr][1:4], tof, GM_sun)
      except ValueError:
        continue
      else:
        tof_days = tof / (24*60*60)
        if l[2] < 100:
          result.append([ephemeris_E[i_dep][0], ephemeris_M[i_arr][0], tof_days, l[0], l[1], l[2], l[3]])
        #result.append([ephemeris_E[i_dep][0], ephemeris_M[i_arr][0], tof_days, l.get_v1()[0], l.get_v2()[0]])
        
with open("result0.txt", "w") as file:
  for res in result:
    file.write(" ".join(map(str, res)) + "\n")