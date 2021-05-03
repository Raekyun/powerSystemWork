# LP 한번에 푸는 방식으로 변경

import numpy as np
import numpy_financial as npf
import pandas as pd
from matplotlib import pyplot as plt


def getPvFromPvWatts(path, load, capacity=None, expectation=None):
    """
    pvWatts에서 제공하는 1년치 1시간 pv 데이터를 load의 index와 일치도록 만듦
    # TODO: PV index를 계절, 월, 일이 맞도록 재조정해야함
    :param load: unit [kWh/hour], type [df], len [8760]
    :param pvPath: type [str]
    :param capacity: capacity [kW]
    :param expectation: kWh/year/kW - 1kW당 1년간 기대 발전량
    :return: pv: unit [kWh/hour], type [df], len [8760]
    """
    pv = pd.read_csv(path, skiprows=16)
    pv = pv['AC System Output (W)']
    pv.name = 'pv'
    pv = pv.drop(8760)
    pv.index = load.index
    pv = pv.divide(1000)
    pv = pd.DataFrame(pv)
    if capacity is not None and expectation is not None:
        pv = pv / pv.sum() * capacity * expectation
    return pv


def getEnergyTable(load, pv, pcsCap, battCap, essMode, peakTimeStart, peakTimeEnd):
    """
    load, pv를 입력받아서 load, pv, ess, net, netNoExport, gridExport 반환
    :param load: pd.DataFrame
    :param pv: pd.DataFrame
    :return: pd.DataFrame
    """
    energyTable = load.copy()
    energyTable['pv'] = pv.pv
    energyTable['netLoadPv'] = energyTable.load - energyTable.pv
    energyTable['peak'] = [True if peakTimeStart <= h < peakTimeEnd else False for h in energyTable.index.hour]
    energyTable.loc[(energyTable.index.dayofweek >= 5), 'peak'] = False
    energyTable['ess'] = calcEss(energyTable, pcsCap, battCap, essMode)
    energyTable['net'] = energyTable.netLoadPv - energyTable.ess
    energyTable['netNoExport'] = energyTable.net.clip(lower=0)
    energyTable['gridExport'] = -energyTable.net.clip(upper=0)
    return energyTable


def calcEss(energyTable, pcsCap=None, battCap=None, essMode=None):
    df = energyTable.copy()
    if not pcsCap or not battCap or not essMode:
        ess = 0
    elif essMode == 'flat':
        ess = np.array([])
        for k in range(365):
            dailyNetLoadPv = df.netLoadPv.iloc[24 * k: 24 * (k + 1)].values
            dailyEss = dailyEssOptFlat(dailyNetLoadPv, pcsCap, 24)
            ess = np.append(ess, dailyEss)
    elif essMode == 'peak':
        ess = np.array([])
        for k in range(365):
            dailyNetLoadPv = df.netLoadPv.iloc[24 * k: 24 * (k + 1)].values
            dailyEss = dailyEssOptPeak(dailyNetLoadPv, pcsCap, battCap, 24)
            ess = np.append(ess, dailyEss)
    return ess


def dailyEssOptFlat(dailyNetLoadPv, pcsCap, n_t):
    from scipy.optimize import linprog
    c = np.zeros(n_t + 1)
    c[0] = 1
    a = np.zeros([n_t, n_t + 1])
    a[:, 0] = 1
    for i in range(n_t):
        for j in range(n_t + 1):
            if i == j:
                a[i, j + 1] = 1
    a *= -1
    a_ub = a
    b_ub = -dailyNetLoadPv
    a_eq = np.ones([1, n_t + 1])
    a_eq[0, 0] = 0
    b_eq = 0
    x_bounds = [(None, None)]
    for k in range(n_t):
        x_bounds.append((-pcsCap, pcsCap))
    res = linprog(c, a_ub, b_ub, a_eq, b_eq, x_bounds)
    dailyEss = np.delete(res.x, 0)
    return dailyEss

def dailyEssOptPeak(dailyNetLoadPv, pcsCap, battCap, n_t):
    from scipy.optimize import linprog
    c = np.zeros(n_t + 3)
    c[0] = 1
    a = np.zeros([3*n_t + 2, n_t + 3])
    a[0:n_t, 0] = -1
    for i in range(n_t):
        for j in range(n_t):
            if i == j:
                a[i, j + 3] = -1
    a[n_t, 1] = 1
    a[n_t + 1:2*n_t + 1, 1] = -1
    for i in range(n_t):
        for j in range(n_t):
            if i >= j:
                a[n_t + 1 + i, j + 3] = 1
    a[2*n_t + 1, 2] = -1
    a[2*n_t + 2:, 2] = 1
    for i in range(n_t):
        for j in range(n_t):
            if i >= j:
                a[2*n_t + 2 + i, j + 3] = -1

    a_ub = a
    a_eq = np.ones([1, n_t + 3])
    a_eq[0,0:3] = 0
    b_eq = 0
    b_s0 = -dailyNetLoadPv
    b_s1 = np.append(battCap, np.zeros(n_t))
    b_s2 = np.append(battCap, np.zeros(n_t))
    b_ub = np.concatenate([b_s0, b_s1, b_s2], axis=0)
    x_bounds = [(None, None), (None, None), (None, None)]
    for k in range(n_t):
        x_bounds.append((-pcsCap, pcsCap))
    res = linprog(c, a_ub, b_ub, a_eq, b_eq, bounds = x_bounds)
    dailyEss = res.x[3:]
    return dailyEss


def calcFixedCharge(supplyCharge, marketFee, meteringCharge):
    return supplyCharge + marketFee + meteringCharge


def calcEnergyCharge(energyTable, price_retailPeak, price_retailOffpeak, lf_retail, price_environ, lf_environ,
                     price_networkPeak, price_networkOffpeak, price_market, lf_market):
    df = energyTable.copy()
    df['price_retail'] = [price_retailPeak if p else price_retailOffpeak for p in df.peak]
    df['lf_retail'] = lf_retail
    df['retailChargeBefore'] = df.load * df.price_retail * df.lf_retail
    df['retailChargeAfter'] = df.netNoExport * df.price_retail * df.lf_retail
    retailChargeBefore = df.retailChargeBefore.sum()
    retailChargeAfter = df.retailChargeAfter.sum()
    df['environChargeBefore'] = df.load * price_environ * lf_environ
    df['environChargeAfter'] = df.netNoExport * price_environ * lf_environ
    environChargeBefore = df.environChargeBefore.sum()
    environChargeAfter = df.environChargeAfter.sum()
    df['price_network'] = [price_networkPeak if p else price_networkOffpeak for p in df.peak]
    df['networkChargeBefore'] = df.load * df.price_network
    df['networkChargeAfter'] = df.netNoExport * df.price_network
    networkChargeBefore = df.networkChargeBefore.sum()
    networkChargeAfter = df.networkChargeAfter.sum()
    df['marketChargeBefore'] = df.load * price_market * lf_market
    df['marketChargeAfter'] = df.netNoExport * price_market * lf_market
    marketChargeBefore = df.marketChargeBefore.sum()
    marketChargeAfter = df.marketChargeAfter.sum()
    energyChargeBefore = retailChargeBefore + environChargeBefore + networkChargeBefore + marketChargeBefore
    energyChargeAfter = retailChargeAfter + environChargeAfter + networkChargeAfter + marketChargeAfter
    energyChargeSaving = energyChargeBefore - energyChargeAfter
    results = {
        'df': df,
        'before': energyChargeBefore,
        'after': energyChargeAfter,
        'saving': energyChargeSaving
    }
    return results


def calcFit(energyTable, price_fit):
    df = energyTable.copy()
    df['fitRevenue'] = df.gridExport * price_fit
    fitRevenue = df.fitRevenue.sum()
    results = {
        'df': df,
        'revenue': fitRevenue
    }
    return results


def calcDemandCharge(energyTable, price_monthlyDemand):
    df = energyTable.copy()
    monthlyDemand = df.resample('1M').max()
    monthlyDemand['price_monthlyDemand'] = price_monthlyDemand
    monthlyDemand['demandChargeBefore'] = monthlyDemand.load * price_monthlyDemand
    monthlyDemand['demandChargeAfter'] = monthlyDemand.netNoExport * price_monthlyDemand
    demandChargeBefore = monthlyDemand.demandChargeBefore.sum()
    demandChargeAfter = monthlyDemand.demandChargeAfter.sum()
    demandChargeSaving = demandChargeBefore - demandChargeAfter
    results = {
        'df': df,
        'before': demandChargeBefore,
        'after': demandChargeAfter,
        'saving': demandChargeSaving
    }
    return results


def getCustomerCostRevenue(load, pv, price_ppa, pcsCap, battCap, essMode):
    annualGen = pv.sum().values[0]
    energyTable = getEnergyTable(load=load,
                                 pv=pv,
                                 pcsCap=pcsCap,
                                 battCap=battCap,
                                 essMode=essMode,
                                 peakTimeStart=7,
                                 peakTimeEnd=23)

    energyCharge = calcEnergyCharge(energyTable=energyTable,  # price unit: $/kWh
                                    price_retailPeak=0.089454,
                                    price_retailOffpeak=0.056275,
                                    lf_retail=1.06837,
                                    price_environ=0.0202,
                                    # price_environ = price_environSrecs + price_environVeecs + price_environLrecs
                                    lf_environ=1.0996,
                                    price_networkPeak=0.0467,
                                    price_networkOffpeak=0.0247,
                                    price_market=0.001897,  # price_market = price_marketAnci + price_marketMarket
                                    lf_market=1.0996)

    demandCharge = calcDemandCharge(energyTable=energyTable,  # price unit: $/kVA/month
                                    price_monthlyDemand=9.99917)

    fit = calcFit(energyTable=energyTable, price_fit=0.1)
    results = {'cost': annualGen * price_ppa,
               'revenue': energyCharge['saving'] + demandCharge['saving'] + fit['revenue']}
    return results


if __name__ == "__main__":
    import time
    start = time.time()
    # calc load
    raw_load = pd.read_csv("load_15min.csv", index_col=1, parse_dates=True, infer_datetime_format=True, dayfirst=True)
    load = raw_load.drop(['Meter', 'ReadingDate', 'KW', 'KVA', 'PF', 'Quality', 'QualityText'], 1)
    load.columns = ['load']
    load.sort_index(inplace=True)
    indexToDrop = load.index[0:96].append(load.index[-192:])
    load.drop(indexToDrop, inplace=True)
    load = load.resample('1H').sum()
    # plt.plot(load)
    # plt.show()

    # calc annual charges and fit
    pvCap = 500
    pcsCap = 150
    battCap = 500
    essMode = 'peak'
    pv = getPvFromPvWatts(path=f'pv_{pvCap}.csv',
                          load=load,
                          capacity=pvCap,
                          expectation=1340)  # expectation means generation expectation [kWh] per 1kW pv for 1 year
    annualGen = pv.sum().values[0]
    energyTable = getEnergyTable(load=load,
                                 pv=pv,
                                 pcsCap=pcsCap,
                                 battCap=battCap,
                                 essMode=essMode,
                                 peakTimeStart=7,
                                 peakTimeEnd=23)

    fixedCharge = calcFixedCharge(supplyCharge=8800,  # price unit: $/year
                                  marketFee=1.326,
                                  meteringCharge=2500)

    energyCharge = calcEnergyCharge(energyTable=energyTable,  # price unit: $/kWh
                                    price_retailPeak=0.089454,
                                    price_retailOffpeak=0.056275,
                                    lf_retail=1.06837,
                                    price_environ=0.0202,
                                    # price_environ = price_environSrecs + price_environVeecs + price_environLrecs
                                    lf_environ=1.0996,
                                    price_networkPeak=0.0467,
                                    price_networkOffpeak=0.0247,
                                    price_market=0.001897,  # price_market = price_marketAnci + price_marketMarket
                                    lf_market=1.0996)

    demandCharge = calcDemandCharge(energyTable=energyTable,  # price unit: $/kVA/month
                                    price_monthlyDemand=9.99917)

    fit = calcFit(energyTable=energyTable, price_fit=0.1)

    # calc customer's cashflow
    price_ppa = 0.09
    customerCashflow = pd.DataFrame(0, index=np.arange(16), columns=['cost', 'revenue', 'cashflow'])

    customerCost = np.zeros(16)
    customerCost[1] = annualGen * price_ppa  # customer's annual cost = annual PPA buying cost
    for i in range(1, 15):
        customerCost[i + 1] = customerCost[i] * 0.995
    customerCashflow['cost'] = customerCost

    customerRevenue = np.zeros(16)  # customer's annual revenue = bill saving + fit
    for i in range(1, 16):
        customerRevenue[i] = getCustomerCostRevenue(load=load,
                                                    pv=pv * 0.995 ** (i - 1),
                                                    price_ppa=0.09,
                                                    pcsCap=pcsCap,
                                                    battCap=battCap,
                                                    essMode=essMode)['revenue']
    customerCashflow['revenue'] = customerRevenue
    customerCashflow['cashflow'] = customerCashflow.revenue - customerCashflow.cost

    customerNpv = round(npf.npv(0.035, customerCashflow.cashflow), 3)
    print(f'customerNpv: {customerNpv}')

    # calc GAIA's cashflow
    capex = pvCap * 1300 + pcsCap * 500 + battCap * 750
    opex = capex * 0.015

    gaiaCashflow = pd.DataFrame(data=0,
                                index=np.arange(16),
                                columns=['revenue', 'depreciation', 'grossProfit', 'cost',
                                         'operatingIncome', 'nonOperatingIncome', 'nonOperatingExpense',
                                         'tax', 'incomeBeforeTax', 'netIncome', 'cashflow'])
    gaiaRevenue = np.zeros(16)
    gaiaRevenue[1] = annualGen * price_ppa
    for i in range(1, 15):
        gaiaRevenue[i + 1] = gaiaRevenue[i] * 0.995
    gaiaCashflow['revenue'] = gaiaRevenue

    for i in range(1, 6):
        gaiaCashflow['depreciation'].iloc[i] = capex / 5

    gaiaCashflow['grossProfit'] = gaiaCashflow.revenue - gaiaCashflow.depreciation

    gaiaCashflow['cost'] = opex
    gaiaCashflow['cost'].iloc[0] = 0

    gaiaCashflow['operatingIncome'] = gaiaCashflow.grossProfit - gaiaCashflow.cost

    gaiaCashflow[
        'incomeBeforeTax'] = gaiaCashflow.operatingIncome + gaiaCashflow.nonOperatingIncome - gaiaCashflow.nonOperatingExpense

    gaiaCashflow['tax'] = gaiaCashflow.incomeBeforeTax * 0.3
    gaiaCashflow.loc[(gaiaCashflow.tax < 0), 'tax'] = 0

    gaiaCashflow['netIncome'] = gaiaCashflow.incomeBeforeTax - gaiaCashflow.tax

    gaiaCashflow['cashflow'] = gaiaCashflow.netIncome + gaiaCashflow.depreciation
    gaiaCashflow['cashflow'].iloc[0] = -capex

    gaiaNpv = round(npf.npv(0.035, gaiaCashflow.cashflow), 3)
    gaiaIrr = round(npf.irr(gaiaCashflow.cashflow), 3)
    print(f'gaiaNpv: {gaiaNpv}, gaiaIrr: {gaiaIrr}')
    round(gaiaCashflow.transpose(), 1).to_csv('../../../Google Drive/My Drive/제안 지원/2021/호주_GAIA/BTM_Wendouree/gaiaCashflow.csv')
    round(customerCashflow.transpose(), 1).to_csv('../../../Google Drive/My Drive/제안 지원/2021/호주_GAIA/BTM_Wendouree/customerCashflow.csv')
    round(energyTable, 1).to_csv('../../../Google Drive/My Drive/제안 지원/2021/호주_GAIA/energyTable.csv')

    # drawings
    energyTable['pv']['2019-09-01'].plot()
    plt.show()
    energyTable['load']['2019-09-01'].plot()
    plt.show()
    energyTable.loc[:, ['load', 'pv', 'netLoadPv']]['2019-09-01'].plot()
    plt.show()
    energyTable.loc[:, ['netLoadPv', 'ess', 'net']]['2019-09-01'].plot()
    plt.show()
    energyTable.loc[:, ['load', 'pv', 'ess']].plot()
    plt.show()

    demand_load = energyTable['load'].resample('1M').max()
    demand_netLoadPv = energyTable['netLoadPv'].resample('1M').max()
    demand_net = energyTable['net'].resample('1M').max()
    label = demand_load.index.values
    x = 1 + np.arange(len(label))
    plt.bar(x-0, demand_load, width=0.2, label='load')
    plt.bar(x+0.2, demand_netLoadPv, width=0.2, label='netLoadPv')
    plt.bar(x+0.4, demand_net, width=0.2, label='net')
    plt.xlabel('month')
    plt.ylabel('load [kWh]')
    plt.legend()
    plt.title('Comparison of peak demands')
    plt.show()
    # results

    totalRevenue = energyCharge['saving'] + demandCharge['saving'] + fit['revenue']
    print(f'totalGeneration: {annualGen} \ntotalRevenue: {totalRevenue}')
    print("time :", time.time() - start)