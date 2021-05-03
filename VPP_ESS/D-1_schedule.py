import numpy as np
from scipy.optimize import linprog


def opt_ess(P_pv, SMP, rhoEss, dT, T, recPrice, pcsCap, battCap, soc0, socMin, socMax, rhoPv, etaC, etaD):
    c = np.zeros(2 * T)
    for t in range(T):
        c[t] = -SMP[t] * dT - rhoPv * recPrice * dT
        c[t + T] = -SMP[t] * dT - rhoEss[t] * recPrice * dT

    A_ub = np.zeros([5 * T, 2 * T])
    b_ub = np.zeros(5 * T)
    for t in range(T):
        A_ub[t, :t + 1] = -etaC / battCap * dT
        A_ub[t, T:t + T + 1] = -1 / etaD / battCap * dT
        b_ub[t] = socMax - soc0
        A_ub[t + T, :t + 1] = etaC / battCap * dT
        A_ub[t + T, T:t + T + 1] = 1 / etaD / battCap * dT
        b_ub[t + T] = -socMin + soc0
        A_ub[t + 2 * T, t] = 1
        A_ub[t + 2 * T, t + T] = 1
        b_ub[t + 2 * T] = pcsCap
        A_ub[t + 3 * T, t] = -1
        A_ub[t + 3 * T, t + T] = -1
        b_ub[t + 3 * T] = pcsCap
        A_ub[t + 4 * T, t] = -1
        A_ub[t + 4 * T, t + T] = -1
        b_ub[t + 4 * T] = P_pv[t]

    A_eq = np.zeros([T + 1, 2 * T])
    for t in range(10):
        A_eq[t, t] = 1
    for t in range(16, 24):
        A_eq[t, t] = 1
    for t in range(10, 16):
        A_eq[t, t + T] = 1
    for t in range(T):
        A_eq[T, t] = etaC / battCap * dT
        A_eq[T, t + T] = 1 / etaD / battCap * dT
    b_eq = np.zeros(T + 1)

    bounds = []
    for t in range(T):
        bounds.append((-pcsCap, 0))
    for t in range(T):
        bounds.append((0, pcsCap))

    x = linprog(c, A_ub, b_ub, A_eq, b_eq, bounds, 'simplex').x
    P_chg = x[:T]
    P_dchg = x[T:2 * T]
    P_ess = P_chg + P_dchg

    soc = np.zeros(T)
    for t in range(T):
        soc[t] = soc0 - (etaC * P_chg[:t + 1]).sum() / battCap - (1 / etaD * P_dchg[:t + 1]).sum() / battCap

    revSmp = SMP * (P_pv + P_ess) * dT
    revRec = (rhoPv * recPrice * (P_pv + P_chg) + rhoEss * recPrice * P_dchg) * dT
    rev = revSmp.sum() + revRec.sum()
    return P_ess, soc, rev


if __name__ == "__main__":
    pvCap = 500
    P_pv = np.array([0, 0, 0, 0, 0, 0, 31, 263, 423, 430, 414, 385, 357, 411, 402, 424, 360, 109, 1, 0, 0, 0, 0, 0])
    SMP = np.array([83.71, 81.69, 81.52, 80.13, 80.03, 80.14, 81.86, 82.45, 83.9, 87.45, 0, 0, 0, 0, 0,
                    0, 86.45, 86.45, 84.7, 84.2, 84.24, 83.5, 83.5, 82.42])
    r_ess = 4
    rhoEss = np.ones(24) * r_ess
    rhoEss[10:16] = 0

    [P_ess, soc, rev] = opt_ess(P_pv=P_pv,
                                SMP=SMP,
                                rhoEss=rhoEss,
                                dT=1,
                                T=int(24 / 1),
                                recPrice=40,
                                pcsCap=500,
                                battCap=2000,
                                soc0=0.2,
                                socMin=0.2,
                                socMax=0.8,
                                rhoPv=1.1,
                                etaC=0.98,
                                etaD=0.98)

    from matplotlib import pyplot as plt

    plt.plot(list(range(1, 25)), P_ess)
    plt.plot(list(range(1, 25)), P_pv)
    plt.show()
    plt.plot(soc)
    plt.show()
    plt.plot(SMP)
    plt.show()
