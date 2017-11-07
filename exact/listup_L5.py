def bit_count32(i):
    i = i - ((i >> 1) & 0x55555555)
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333)
    i = (i + (i >> 4)) & 0x0f0f0f0f
    i = i + (i >> 8)
    i = i + (i >> 16)
    return i & 0x3f


def bit_count16(i):
    i = i - ((i >> 1) & 0x5555)
    i = (i & 0x3333) + ((i >> 2) & 0x3333)
    i = (i + (i >> 4)) & 0x0f0f
    i = i + (i >> 8)
    return i & 0x1f


def main():
    units = ['11001', '11100', '01110', '00111', '10011']
    filters = []
    for a in units:
        for b in units:
            s = '0b'
            for aa in a:
                if aa == '1':
                    s += b
                else:
                    s += '00000'
            filters.append(eval(s))
    d = {}
    for a in range(2**25):
        b = 0
        for f in filters:
            b <<= 1
            if bit_count32(a & f) in [3, 4]:
                b += 1
        d[b] = d.get(b, 0) + 1
    for b in sorted(d):
        print(bin(b+0x2000000)[3:], d[b])


if __name__ == '__main__':
    main()
