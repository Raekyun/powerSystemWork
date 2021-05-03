import pandas as pd
import numpy as np
from scipy.optimize import linprog
from matplotlib import pyplot as plt


def dailyEssOpt(dailySmp, pcsCap, battCap):
    c = -dailySmp
    a_eq = np.ones([1, 24])
    b_eq = battCap
    x_bounds = []
    for k in range(24):
        x_bounds.append((0, pcsCap))
    res = linprog(c=c, A_eq=a_eq, b_eq=b_eq, bounds=x_bounds)
    return res.x


def preprocessEtaSolarNetData(df, netIndex):
    net = df.iloc[netIndex, :]
    columns = []
    for k in range(24):
        columns.append(f'{k}')
    net.columns = columns
    index = pd.date_range(net.index[0], net.index[-1], freq='H')
    index = pd.date_range(index[0], index[-1] + index[-1].freq * 23, freq='H')
    values = net.values
    values = values.flatten()
    columns = ['net']
    net = pd.DataFrame(values, index, columns)
    return net


def preprocessEtaSolarEssDData(df, dchgIndex):
    essD = df.iloc[dchgIndex, :]
    columns = []
    for k in range(24):
        columns.append(f'{k}')
    essD.columns = columns
    index = pd.date_range(essD.index[0], essD.index[-1], freq='H')
    index = pd.date_range(index[0], index[-1] + index[-1].freq * 23, freq='H')
    values = essD.values
    values = values.flatten()
    columns = ['essD']
    essD = pd.DataFrame(values, index, columns)
    return essD


def preprocessEtaSolarEssCData(df, chgIndex):
    essC = df.iloc[chgIndex, :]
    columns = []
    for k in range(24):
        columns.append(f'{k}')
    essC.columns = columns
    index = pd.date_range(essC.index[0], essC.index[-1], freq='H')
    index = pd.date_range(index[0], index[-1] + index[-1].freq * 23, freq='H')
    values = essC.values
    values = values.flatten()
    columns = ['essC']
    essC = pd.DataFrame(values, index, columns)
    return essC


def preprocessEpsisData(df):
    columns = []
    for k in range(24):
        columns.append(f'{k}')
    df.columns = columns
    df = df.sort_index()
    index = pd.date_range(df.index[0], df.index[-1], freq='H')
    index = pd.date_range(index[0], index[-1] + index[-1].freq * 23, freq='H')
    values = df.values
    values = values.flatten()
    columns = ['smp']
    smp = pd.DataFrame(values, index, columns)
    return smp


def socStabilizer(soc0, socB, pv):
    if soc0 < socB:
        return
    t = 16
    dchg = []
    eTot = (soc0 - socB) * battCap * etaD / 100
    print(eTot)
    print(0.7*pvCap - pv[t])
    pMax = min(0.7*pvCap - pv[t], pcsCap)
    while eTot >= pMax:
        dchgTemp = pMax
        print(dchgTemp)
        dchg.append(dchgTemp)
        soc0 = soc0 - dchgTemp/etaD/battCap * 100
        eTot = (soc0 - socB) * battCap * etaD / 100
        t += 1
        pMax = min(0.7*pvCap - pv[t], pcsCap)
    dchgLast = eTot
    dchg.append(dchgLast)
    t += 1
    soc0 = soc0 - dchgLast/etaD/battCap * 100
    return dchg, t, soc0


if __name__ == "__main__":
    data = pd.read_excel('청명 20201001~31 데이터.xlsx', index_col=0, usecols="A,E:AB", engine='openpyxl')
    smp = pd.read_excel('시간별SMP 20201001~31.xlsx', index_col=0, usecols="A:Y", engine='openpyxl')
    dailySoc = [80.00, 70.90, 58.90, 75.70, 80.00, 80.00, 80.00, 79.90, 80.00, 80.00, 75.40,
                80.00, 80.00, 80.00, 80.00, 59.50, 80.00, 80.00, 80.00, 80.00, 26.10, 37.40,
                80.00, 80.00, 80.00, 79.90, 80.00, 80.00, 80.00, 80.00, 80.00]
    socB = 40
    pvCap = 2747.52
    pcsCap = 2000
    battCap = 7488.24
    etaC = 0.95
    etaD = 0.916

    netIndex = list(range(31))
    energyTable = preprocessEtaSolarNetData(data, netIndex)
    dchgIndex = list(range(31, 93, 2))
    chgIndex = list(range(32, 94, 2))
    energyTable['essD'] = preprocessEtaSolarEssDData(data, dchgIndex)
    energyTable['essC'] = preprocessEtaSolarEssCData(data, chgIndex)
    energyTable['ess'] = energyTable['essD'] - energyTable['essC']
    energyTable['pv'] = energyTable['net'] - energyTable['ess']

    energyTableDailySum = energyTable.resample('1D').sum()
    energyTableDailySum['soc'] = dailySoc

    energyTable['smp'] = preprocessEpsisData(smp)

    [dchg, t0, soc0] = socStabilizer(80, socB, energyTable['pv']['2020-10-01'].values/1000)

    # energyTable['ess_uniform'] = [pcsCap if 17 <= h < 17 + battCap / pcsCap else 0 for h in energyTable.index.hour]
    # energyTable['revenue_uniform'] = energyTable['smp'] * energyTable['ess_uniform']
    #
    # analysisPeriod = int(energyTable.index.shape[0] / 24)
    # ess = np.array([])
    # for k in range(analysisPeriod):
    #     dailySmp = energyTable['smp_for_calc'].iloc[k * 24:(k + 1) * 24]
    #     dailyEss = dailyEssOpt(dailySmp, pcsCap, battCap)
    #     ess = np.append(ess, dailyEss)
    #
    # energyTable['ess_opt'] = ess
    # energyTable['revenue_opt'] = energyTable['smp'] * energyTable['ess_opt']
    #
    # print(f'Total revenue with uniform ESS: {energyTable["revenue_uniform"].sum()}')
    # print(f'Total revenue with optimized ESS: {energyTable["revenue_opt"].sum()}')
    #
    # date = '2020-01-05'
    # right_ax_ylim = energyTable.loc[:, 'smp'][date].max() + 5
    # ax = energyTable[['ess_uniform', 'ess_opt', 'smp']].loc[date].plot(secondary_y='smp', ylim=[0, 1500])
    # ax.right_ax.set_ylim(0, right_ax_ylim)
    #
    # plt.show()
    #
    # energyTable.to_csv('energyTable.csv')
