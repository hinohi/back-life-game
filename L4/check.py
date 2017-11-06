from pathlib import Path
from math import exp


def get_num(path):
    data = {}
    for line in path.open():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        tag, f = line.split()
        data[tag] = float(f)

    if '0' in data:
        max_f = max(data.values())
        denominator = sum(exp(f - max_f) for f in data.values())
        prob = exp(data['0'] - max_f) / denominator
        return prob*2**16
    else:
        return None


def get_all():
    d = {}
    for line in open('all.txt'):
        spin, n = line.strip().split()
        d[spin] = int(n)
    return d


def main():
    d = get_all()
    data_path = Path('../back-mcmcmc/data/')
    e0 = 0
    e1 = 0
    max_error = 0.0
    for hist in data_path.glob('L0004_seed*_p0004.hist'):
        n = get_num(hist)
        spin = ''
        for line in open(data_path / (hist.name[:-5] + '.spin')):
            spin += line.strip()
        if n is None:
            if spin in d:
                e1 += 1
        else:
            if spin not in d:
                e0 += 1
            else:
                error = abs(d[spin] - n) / d[spin]
                max_error = max(max_error, error)
    print('ないものをあるとした: %i' % e0)
    print('あるものを見逃した: %i' % e1)
    print('最大誤差: %e' % max_error)


if __name__ == '__main__':
    main()
