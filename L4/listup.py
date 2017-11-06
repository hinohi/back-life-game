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
    filters = [
        0b1101_1101_0000_1101,
        0b1110_1110_0000_1110,
        0b0111_0111_0000_0111,
        0b1011_1011_0000_1011,
        0b1101_1101_1101_0000,
        0b1110_1110_1110_0000,
        0b0111_0111_0111_0000,
        0b1011_1011_1011_0000,
        0b0000_1101_1101_1101,
        0b0000_1110_1110_1110,
        0b0000_0111_0111_0111,
        0b0000_1011_1011_1011,
        0b1101_0000_1101_1101,
        0b1110_0000_1110_1110,
        0b0111_0000_0111_0111,
        0b1011_0000_1011_1011,
    ]
    d = {}
    for a in range(2**16):
        b = 0
        for f in filters:
            b <<= 1
            if bit_count16(a & f) in [3, 4]:
                b += 1
        d.setdefault(b, []).append(a)
    for b in sorted(d):
        print(bin(b+0x010000)[3:], len(d[b]))


if __name__ == '__main__':
    main()
