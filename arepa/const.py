"""
"""

MSOL = 1.988e+33
PC = 3.0856e+18


class PARTS:
    GAS = 0
    DM = 1
    STAR = 4
    BH = 5


class CONV_ILL_TO_SOL:
    MASS = 1.989e43 / MSOL
    DIST = 3.085678e24 / (1e6*PC)
