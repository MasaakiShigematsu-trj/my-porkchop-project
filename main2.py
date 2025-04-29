## 楕円遷移に限定するのであれば、zの初期値は0でよいかも？（わざわざ-100から始める必要ない？）

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

def lambert(GM, R1, R2, tof, string):
    #R1 = X1[:3]  # X1（6次元）の最初の3つの要素を取り出す
    r1 = np.linalg.norm(R1)  # ノルムを計算する


    #R2 = X2[:3]  # X2（6次元）の最初の3つの要素を取り出す
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
    z = -100  # 初期値
    print(y(-100))
    while y(z) < 0:
       z = z + 0.1
       #print(z, y(z))
    print("initial value 0")
    print(z, y(z))

    ## 次にF(z)が正の値になるまでzを増加させる
    while F(z) < 0:
       z = z + 0.1
       #print(z, F(z))
    print("initial value 1")
    print(z, F(z))

    # z = 0.1

    imax = 1000
    for i in range(imax):
        z = z - F(z)/dF(z)  # ニュートン法の更新式
        if y(z) < 0:
            print("y(z) is negative")
        elif abs(F(z)/dF(z)) < 1e-13:  # 収束条件
            break
        else:
            print(i, z, y(z))
            print(abs(F(z)/dF(z)))        
    
    print(i, z)
    print(abs(F(z)/dF(z)))
    

    ## ラグランジュの係数の計算
    f = 1 - y(z)/r1
    g = A*np.sqrt(y(z)/GM)
    gdot = 1 - y(z)/r2

    ## 速度ベクトルの計算
    V1 = (1/g)*(R2 - f*R1)
    V2 = (1/g)*(gdot*R2 - R1)

    return V1, V2

GM_earth = 398600.4354360959 # 地球の重力定数 [km^3/s^2]
GM_sun = 1.32712440041279419e11 # 太陽の重力定数 [km^3/s^2]
#X1 = np.array([0.01, -6778.0, 0.0, 0.0, 0.0, 0.0])
#X2 = np.array([0.01, 42164.0, 0.0, 0.0, 0.0, 0.0])
#R1 = np.array([1.621843172290461E+07, 1.465106266682866E+08, -4.394324564635754E+02])
JD1 = 2462851.500000000 # A.D. 2030-Dec-16 00:00:00.0000 TDB 
R1 = np.array([1.621843172290461E+07, 1.465106266682866E+08, -4.394324564635754E+02])

#R2 = np.array([1.463122336310571E+08, 1.635698047472437E+08, -1.411212307718471E+05])
#R2 = np.array([-2.096570358345060E+08, 1.335054354669643E+08, 7.952008874344416E+06])
#R2 = np.array([1.107482404379682E+08, -1.775947265582940E+08, -6.435435377618231E+06])
JD2 = 2463132.500000000 # A.D. 2031-Sep-23 00:00:00.0000 TDB 
R2 = np.array([1.107482404379682E+08, -1.775947265582940E+08, -6.435435377618231E+06])

#tof = (2461931.500000000 - 2461771.500000000)*24*60*60  # [s]
#tof = (2462136.500000000 - 2461771.500000000)*24*60*60  # [s]
#tof = (2463132.500000000 - 2462851.500000000)*24*60*60  # [s]
tof = (JD2 - JD1)*24*60*60  # [s]

string = 'pro'

l = lambert(GM_sun, R1, R2, tof, string)

print("my lambert")
print(l[0], l[1])


l_pykep = pk.lambert_problem(R1, R2, tof, GM_sun)
V1_pykep = l_pykep.get_v1()[0]
V2_pykep = l_pykep.get_v2()[0]
print("pykep lambert")
print(V1_pykep, V2_pykep)