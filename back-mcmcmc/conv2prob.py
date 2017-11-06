from math import exp
import sys


data = {}
target = sys.argv[2]
L = int(sys.argv[3])
for line in open(sys.argv[1]):
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    tag, f = line.split()
    data[tag] = float(f)

max_f = max(data.values())
denominator = sum(exp(f - max_f) for f in data.values())
if target in data:
    prob = exp(data[target] - max_f) / denominator
    print(prob, prob*2**(L*L))
else:
    print('nan')
