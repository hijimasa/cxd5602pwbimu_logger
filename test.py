import os
import subprocess
import time
import struct
import numpy as np

GRAVITY_AMOUNT = 9.80665
EARTH_ROTATION_SPEED_AMOUNT = 7.2921159e-5

MESUREMENT_FREQUENCY = 960
GYRO_NOISE_DENSITY = 1.0e-3 * np.pi / 180
ACCEL_NOISE_DENSITY = 14.0e-6 * GRAVITY_AMOUNT
GYRO_NOISE_AMOUNT = GYRO_NOISE_DENSITY * (MESUREMENT_FREQUENCY**0.5) * 40
ACCEL_NOISE_AMOUNT = ACCEL_NOISE_DENSITY * (MESUREMENT_FREQUENCY**0.5) * 40
ACCEL_BIAS_DRIFT = 4.43e-6 * GRAVITY_AMOUNT
GYRO_BIAS_DRIFT = 0.39 * np.pi / 180

def hex_to_float(hex_str, endian='big'):
    """
    16進数文字列（例: "41cb8300"）を IEEE754 の float に変換する。
    endian が 'big' の場合はネットワークバイトオーダー（ビッグエンディアン）で変換。
    """
    # 余計な空白があれば除去
    hex_str = hex_str.strip()
    # 16進数文字列をバイト列に変換
    b = bytes.fromhex(hex_str)
    if endian == 'big':
        return struct.unpack('!f', b)[0]
    elif endian == 'little':
        return struct.unpack('<f', b)[0]
    else:
        raise ValueError("endian は 'big' か 'little' のどちらかにしてください")

def diff_quaternion(q, omega):
    w, x, y, z = q
    omega_x, omega_y, omega_z = omega
    dqdt = [ 
        0.5 * (-x * omega_x - y * omega_y - z * omega_z),
        0.5 * ( w * omega_x + y * omega_z - z * omega_y),
        0.5 * ( w * omega_y - x * omega_z + z * omega_x),
        0.5 * ( w * omega_z + x * omega_y - y * omega_x)
    ]
    return dqdt

def runge_kutta_update(q, omega, h):
    k1 = diff_quaternion(q, omega)
    q1 = [q_i + (h/2)*k1_i for q_i, k1_i in zip(q, k1)]
    
    k2 = diff_quaternion(q1, omega)
    q2 = [q_i + (h/2)*k2_i for q_i, k2_i in zip(q, k2)]
    
    k3 = diff_quaternion(q2, omega)
    q3 = [q_i + h*k3_i for q_i, k3_i in zip(q, k3)]
    
    k4 = diff_quaternion(q3, omega)
    
    q_next = [q_i + (h/6)*(k1_i + 2*k2_i + 2*k3_i + k4_i) 
              for q_i, k1_i, k2_i, k3_i, k4_i in zip(q, k1, k2, k3, k4)]
    
    # 正規化
    norm = sum([comp**2 for comp in q_next]) ** 0.5
    q_next = [comp / norm for comp in q_next]
    
    return q_next

def cross_product(a, b):
    return np.cross(a, b)

def update_vector_rk4(v, omega, h):
    # RK4 ステップ
    k1 = cross_product(omega, v)
    k2 = cross_product(omega, v + (h/2)*k1)
    k3 = cross_product(omega, v + (h/2)*k2)
    k4 = cross_product(omega, v + h*k3)
    
    # RK4 による更新
    v_next = v + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    return v_next

biased_velocity = 0.0
zero_velocity_counter = 0
# ゼロ速度補正関数
def zero_velocity_correction(velocity, dt):
    global biased_velocity
    global zero_velocity_counter
    biased_velocity += ACCEL_BIAS_DRIFT * dt
    if abs(velocity[0]) < biased_velocity and abs(velocity[1]) < biased_velocity and abs(velocity[2]) < biased_velocity:
        zero_velocity_counter += 1
    else:
        zero_velocity_counter = 0
    if zero_velocity_counter > MESUREMENT_FREQUENCY:
        biased_velocity = 0.0
        zero_velocity_counter = 0
        velocity = [0.0, 0.0, 0.0]
    return velocity

# FIFOのパスを指定
FIFO_PATH = "/tmp/sony_imu_fifo"

# FIFOが存在しなければ作成する
if not os.path.exists(FIFO_PATH):
    os.mkfifo(FIFO_PATH)

try:
    # screenセッションをFIFOをログファイルとして起動
    # オプションの順序に注意（-dmS の後に -L -Logfile など）
    subprocess.run(["screen", "-dmS", "sony_imu", "-L", "-Logfile", FIFO_PATH, "/dev/ttyUSB0", "115200"])

    # 少し待機して接続を確立
    time.sleep(2)

    # コマンド送信
    subprocess.run(["screen", "-S", "sony_imu", "-X", "stuff", f"pwbimu_logger -s {MESUREMENT_FREQUENCY} -a 8 -g 2000 -f 1 -o uart\n"])

    # FIFOからカンマ区切りの8個の4バイトデータを読み出して表示
    with open(FIFO_PATH, "r") as f:
        while True:
            is_updated = False

            line = f.readline().strip()
            if line:
                # lineをカンマで分割
                raw_data = line.split(",")
                if len(raw_data) != 18:
                    continue
                # 最初は32bitの整数値、あとの7つはfloat値に変換して表示
                try:
                    data = [int(raw_data[0], 16)] + [hex_to_float(x) for x in raw_data[1:]]
                    print(data)
                #すべてのエラーを待ち受け
                except Exception as e:
                    #print(e)
                    continue
                is_updated = True
            

except KeyboardInterrupt:
    pass
finally:
    # screenセッションの終了（セッション番号を取得してkillする例）
    session_info = subprocess.check_output(["screen", "-ls"]).decode("utf-8")
    for line in session_info.splitlines():
        if "sony_imu" in line:
            session_id = line.split(".")[0].strip()
            subprocess.run(["screen", "-S", f"{session_id}.sony_imu", "-X", "quit"])
            break
    # FIFOの削除
    os.remove(FIFO_PATH)
