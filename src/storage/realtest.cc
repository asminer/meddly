
#include <cstdio>


union word32 {
    unsigned u;
    float f;
};

union word64 {
    unsigned long u;
    double d;
};

void showbits(float f)
{
    printf("%17.8e:\n\t", f);

    word32 w;
    w.f = f;

    bool sign;
    unsigned expo=0;
    unsigned frac=0;
    unsigned b=32;
    do {
        b--;
        unsigned mask = (0x000001 << b);
        bool bit = (w.u & mask);
        fputc(bit+'0', stdout);
        if (31 == b) {
            fputc(' ', stdout);
            sign = bit;
            continue;
        }
        if (23 <  b) {
            expo *= 2;
            expo += bit;
            continue;
        }
        if (23 == b) {
            expo *= 2;
            expo += bit;
            fputc(' ', stdout);
            continue;
        }
        frac *= 2;
        frac += bit;
    } while (b);
    printf(" = %c Exp %u Frac %x\n", (sign ? '-' : '+'), expo, frac);
}

float buildFloat(bool sign, unsigned expo, unsigned frac)
{
    word32 w;
    w.u = (sign ? 0x80000000 : 0) |
          ( (expo & 0xFF) << 23 ) |
          ( frac & 0x7FFFFF);
    return w.f;
}

int main()
{
    showbits(1.0f);
    showbits(1.9999999f);
    showbits(2.0f);

    showbits(3.2f);

    showbits(buildFloat(0, 0x1f, 0x555555));
    showbits(buildFloat(1, 0x1f, 0x555555));

    showbits(buildFloat(0, 0, 0));
    showbits(buildFloat(1, 0, 0));
    showbits(buildFloat(0, 0xff, 0));
    showbits(buildFloat(1, 0xff, 0));
    return 0;
}
