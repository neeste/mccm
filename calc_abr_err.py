import sys

data = """
 0.50  20  19.2 11.4
 0.50  40  39.2 10.2
 0.50  60  59.2  9.3
 0.50  80  79.2  8.7
 1.00  20  20.0  7.5
 1.00  40  40.0  6.4
 1.00  60  60.0  5.8
 1.00  80  80.0  5.3
 2.00  20  20.4  5.1
 2.00  40  40.4  4.0
 2.00  60  60.4  3.5
 2.00  80  80.4  3.1
 4.00  20  20.8  3.3
 4.00  40  40.8  2.7
 4.00  60  60.8  2.1
 4.00  80  80.8  1.9
"""

mse = 0.0
count = 0
b=12.9; c=5; d=0.413

for line in data.strip().split('\n'):
    parts = line.split()
    if len(parts) == 4:
        fr = float(parts[0])
        lv = float(parts[1])
        lat = float(parts[3])
        i = lv / 100.0
        abr = b * (c**(-i)) * (fr**(-d))
        err = lat - abr
        mse += err**2
        count += 1
        print(f"fr={fr}, lv={lv}, lat={lat}, abr={abr:.3f}, err={err:.3f}")

if count > 0:
    print(f"MSE: {mse/count:.3f}")

