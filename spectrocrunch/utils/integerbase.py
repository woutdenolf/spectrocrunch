# -*- coding: utf-8 -*-


class integerbase:
    def __init__(
        self,
        digs=[
            "0",
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
        ],
    ):
        self.digs = digs
        self.base = len(self.digs)

    def int2base(self, x):
        if x < 0:
            sign = -1
        elif x == 0:
            return self.digs[0]
        else:
            sign = 1
        x *= sign
        digits = []
        while x:
            digits.append(self.digs[x % self.base])
            x //= self.base
        if sign < 0:
            digits.append("-")
        digits.reverse()
        return "".join(digits)

    def base2int(self, x):
        y = x[::-1]
        val = 0
        for i in range(len(y)):
            val += self.digs.index(y[i]) * self.base**i
        return val
