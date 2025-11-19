import sys
import math

def main():
    if len(sys.argv) != 2:
        print("Usage: python calc_pi_list.py N")
        return

    try:
        N = int(sys.argv[1])
        if N <= 0:
            raise ValueError
    except ValueError:
        print("N must be a positive integer.")
        return

    PI = math.pi
    step = PI / N

    # 生成从 0 到 PI 的 N+1 个点
    for i in range(N):
        value = step * i
        print(f"{value:.6f}")

if __name__ == "__main__":
    main()
