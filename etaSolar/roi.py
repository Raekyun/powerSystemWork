import pandas as pd
import numpy as np
from scipy.optimize import linprog
from matplotlib import pyplot as plt
from matplotlib import dates as mdates


def dailyEssOpt(pv_to_opt, smp_to_opt, t_to_opt, soc0):
    c = -smp_to_opt
    A_eq = np.ones([1, t_to_opt]) / etaD / battCap
    b_eq = (soc0 - socMin)/100
    A_ub = np.zeros([2*t_to_opt, t_to_opt])
    b_ub = np.zeros(2*t_to_opt)
    for t in range(t_to_opt):
        A_ub[t, :t + 1] = 1/etaD/battCap
        b_ub[t] = (-socMin + soc0)/100
        A_ub[t + t_to_opt, t] = 1
        b_ub[t + t_to_opt] = 0.7 * pvCap - pv_to_opt[t]

    x_bounds = []
    for k in range(t_to_opt):
        x_bounds.append((0, pcsCap))
    res = linprog(c=c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=x_bounds)
    dchgOpt = res.x
    socOpt = np.zeros(t_to_opt)
    for t in range(t_to_opt):
        socOpt[t] = soc0 - (1 / etaD * dchgOpt[:t + 1]).sum() / battCap * 100
    revOpt = smp_to_opt * dchgOpt
    return dchgOpt, socOpt, revOpt


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


def socStabilizer(socInit, socB, pv, smp_to_bal):
    """
    SOC 안정화를 수행하는 block
    :param soc0: t=16 시점의 soc
    :param socB: BMS가 soc balancing을 수행하는 최소 soc
    :param pv: 해당 일의 pv 발전량 패턴. unit = Wh. time scale = 1h
    :return:
        dchg: 방전패턴
        t: balancing 종료 후 시간
        soc0: balancing 종료 후 soc
    """
    t = t_start
    if socInit <= socB:
        return [], t, []
    dchg = []
    socBal = []
    eTot = (socInit - socB) * battCap * etaD / 100
    pMax = min(0.7*pvCap - pv[t], pcsCap)
    socTemp = socInit
    if smp_to_bal[0] == np.sort(smp_to_bal)[-1] or smp_to_bal[0] == np.sort(smp_to_bal)[-2]:
        dchgTemp = pMax
        dchg.append(dchgTemp)
        socTemp = socTemp - dchgTemp/etaD/battCap * 100
        socBal.append(socTemp)
        eTot = (socTemp - socB) * battCap * etaD / 100
        t += 1
        pMax = min(0.7*pvCap - pv[t], pcsCap)
        if eTot <= 0:
            return dchg, t, socBal

    while eTot >= pMax:
        dchgTemp = pMax
        dchg.append(dchgTemp)
        socTemp = socTemp - dchgTemp/etaD/battCap * 100
        socBal.append(socTemp)
        eTot = (socTemp - socB) * battCap * etaD / 100
        t += 1
        pMax = min(0.7*pvCap - pv[t], pcsCap)
        if t >= 24:
            break
    dchgLast = eTot
    dchg.append(dchgLast)
    soc0 = socTemp - dchgLast/etaD/battCap * 100
    socBal.append(soc0)
    t += 1
    return dchg, t, socBal


if __name__ == "__main__":
    # 데이터 입력 부분
    data = pd.read_excel('청명 20201001~31 데이터.xlsx', index_col=0, usecols="A,E:AB", engine='openpyxl')
    smp = pd.read_excel('시간별SMP 20201001~31.xlsx', index_col=0, usecols="A:Y", engine='openpyxl')
    dailySoc = [80.00, 70.90, 58.90, 75.70, 80.00, 80.00, 80.00, 79.90, 80.00, 80.00, 75.40,
                80.00, 80.00, 80.00, 80.00, 59.50, 80.00, 80.00, 80.00, 80.00, 26.10, 37.40,
                80.00, 80.00, 80.00, 79.90, 80.00, 80.00, 80.00, 80.00, 80.00]
    socB = 70
    socMin = 0
    socMax = 80
    t_start = 16
    t_total = 4
    pvCap = 2747.52
    pcsCap = 2000
    battCap = 7488.24
    etaC = 0.9644
    etaD = 0.916

    # 데이터 전처리 및 dataframe 생성
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

    # balancing
    days_to_opt = len(dailySoc)
    energyTable['dchg_opt'] = 0
    energyTable['soc_opt'] = 0
    energyTable['dchg_uniform'] = 0
    energyTable['soc_uniform'] = 0
    for d in range(days_to_opt):
        socInit = dailySoc[d]
        targetDate = energyTableDailySum.index[d].strftime('%Y-%m-%d')
        pv_to_bal = energyTable['pv'][targetDate].values/1000
        smp_to_bal = energyTable['smp'][targetDate].iloc[t_start:].values
        [dchgBal, t0, socBal] = socStabilizer(socInit, socB, pv_to_bal, smp_to_bal)

        t_to_opt = 24 - t0
        pv_to_opt = pv_to_bal[t0:]
        smp_to_opt = energyTable['smp'][targetDate].iloc[t0:].values
        if socBal:
            soc0 = socBal[-1]
        else:
            soc0 = socInit
        [dchgOpt, socOpt, rev] = dailyEssOpt(pv_to_opt, smp_to_opt, t_to_opt, soc0)

        dchg_daily_opt = np.append(dchgBal, dchgOpt)
        soc_daily_opt = np.append(socBal, socOpt)
        energyTable.loc[:, 'dchg_opt'].iloc[d*24+t_start:(d+1)*24] = dchg_daily_opt * \
                                                                     energyTable['essD'][f'2020-10-{d+1}'].sum() \
                                                                     / dchg_daily_opt.sum()
        energyTable.loc[:, 'soc_opt'].iloc[d*24+t_start:(d+1)*24] = soc_daily_opt
        energyTable.loc[:, 'soc_opt'].iloc[d*24 + 15] = energyTableDailySum['soc'].iloc[d]

        eTot = socInit * battCap * etaD * 10
        dchg_daily_uniform = np.zeros(t_total)
        dchg_daily_uniform = dchg_daily_uniform + eTot/t_total
        soc_daily_uniform = []
        for t in range(t_total):
            soc_daily_uniform.append(socInit - t*eTot/t_total/etaD/battCap/10)
        energyTable.loc[:, 'dchg_uniform'].iloc[d*24+t_start+2:(d+1)*24-2] = dchg_daily_uniform
        energyTable.loc[:, 'soc_uniform'].iloc[d*24+t_start+2:(d+1)*24-2] = soc_daily_uniform

    energyTable['rev_etaSolar'] = energyTable['smp'] * energyTable['essD'] / 1000
    energyTable['rev_opt'] = energyTable['smp'] * energyTable['dchg_opt'] / 1000
    energyTable['rev_uniform'] = energyTable['smp'] * energyTable['dchg_uniform'] / 1000

    energyTableDailySum.index = energyTableDailySum.index.to_series().dt.date
    energyTableDailySum[['essD', 'essC']].plot(kind='bar')
    plt.show()

    date_to_plot = '2020-10-16'

    myFmt = mdates.DateFormatter('%d:%H')
    f, axarr = plt.subplots(4)
    axarr[1].set_ylim(0,pcsCap)
    energyTable['smp'].loc[date_to_plot].plot(legend=True, style='r', marker='o', ax=axarr[0])
    (energyTable['essD']/1000).loc[date_to_plot].plot(ax=axarr[1], kind='bar', legend=True).xaxis.set_major_formatter(myFmt)
    (energyTable['dchg_opt']/1000).loc[date_to_plot].plot(ax=axarr[2], kind='bar', legend=True).xaxis.set_major_formatter(myFmt)
    energyTable['soc_opt'].loc[date_to_plot].plot(ax=axarr[3], style='y', legend=True).xaxis.set_major_formatter(myFmt)
    plt.show()

    (energyTable['pv']/1000).loc[date_to_plot].plot(legend=True, color='red', kind='bar').xaxis.set_major_formatter(myFmt)
    plt.show()

    print(f'Total revenue with etaSolar ESS: {energyTable["rev_etaSolar"].sum()}')
    print(f'Total revenue with optimized ESS: {energyTable["rev_opt"].sum()}')
    print(f'Total increased revenue [KRW]: {energyTable["rev_opt"].sum() - energyTable["rev_etaSolar"].sum()}')
    print(f'Total increased revenue [%]: {(energyTable["rev_opt"].sum() - energyTable["rev_etaSolar"].sum()) / energyTable["rev_etaSolar"].sum() * 100}')



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
