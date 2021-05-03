import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


def get_pv_from_PVWatts(pvPath, load, capacity=None, expectation=None):
    """
    pvWatts에서 제공하는 1년치 1시간 pv 데이터를 load의 index와 일치도록 만듦
    # TODO: PV index를 계절, 월, 일이 맞도록 재조정해야함
    :param load: unit [kWh/hour], type [df], len [8760]
    :param pvPath: type [str]
    :param capacity: capacity [kW]
    :param expectation: kWh/year/kW - 1kW당 1년간 기대 발전량
    :return: pv: unit [kWh/hour], type [df], len [8760]
    """
    pv = pd.read_csv(pvPath, skiprows=16)
    pv = pv['AC System Output (W)']
    pv.name = 'pv'
    pv = pv.drop(8760)
    pv.index = load.index
    pv = pv.divide(1000)
    pv = pd.DataFrame(pv)
    if capacity is not None and expectation is not None:
        pv = pv / pv.sum() * capacity * expectation
    return pv


def get_energy_table(load, pv, peakTime1, peakTime2):
    """
    load, pv를 입력받아서 load, pv, net load, reverse flow를 계산하여 반환
    :param load: pd.Series
    :param pv: pd.Series
    :return: pd.dataFrame (load, pv, net load, reverse flow)
    """
    merged = load.copy()
    merged['pv'] = pv.pv
    merged['net'] = merged.load - merged.pv
    merged['net_no_exports'] = merged.net.clip(lower=0)
    merged['grid_exports'] = -merged.net.clip(upper=0)
    merged['peak'] = [True if peakTime1 <= h < peakTime2 else False for h in merged.index.hour]
    net_posi = merged['net'] > 0
    print(net_posi.head(10))
    peak = merged['peak']
    merged['net_peak_posi'] = merged[(net_posi & peak)].net
    merged['net_off_peak_posi'] = merged[(net_posi) & (peak == False)].net
    return merged


def add_retail_charge(energy_table, price_peak, price_off_peak, lf_retail):
    energy_table['price_retail'] = [price_peak if p else price_off_peak for p in energy_table.peak]
    energy_table['lf_retail'] = lf_retail
    energy_table['retail_charge_before'] = energy_table.load * energy_table.price_retail * energy_table.lf_retail
    energy_table['retail_charge_after'] = energy_table.net_no_exports * energy_table.price_retail \
                                          * energy_table.lf_retail
    energy_table['retail_charge_saving'] = energy_table.retail_charge_before - energy_table.retail_charge_after
    return energy_table


def add_fit(energy_table, price_fit):
    energy_table['fit_revenue'] = energy_table.grid_exports * price_fit
    return energy_table


def add_environ(energy_table, price_environ, lf_environ):
    energy_table['environ_charge_before'] = energy_table.load * price_environ * lf_environ
    energy_table['environ_charge_after'] = energy_table.net_no_exports * price_environ * lf_environ
    energy_table['environ_charge_saving'] = energy_table.environ_charge_before - energy_table.environ_charge_after
    return energy_table


def add_network_energy_charge(energy_table, price_network_peak, price_network_off_peak):
    energy_table['price_network_energy'] = [price_network_peak if p else price_network_off_peak
                                            for p in energy_table.peak]
    energy_table['network_energy_charge_before'] = energy_table.load * energy_table.price_network_energy
    energy_table['network_energy_charge_after'] = energy_table.net_no_exports * energy_table.price_network_energy
    energy_table['network_energy_charge_saving'] = energy_table.network_energy_charge_before - \
                                                   energy_table.network_energy_charge_after
    return energy_table


def get_demand_charge_table(energy_table, price_demand):
    demand_monthly = energy_table.resample('1M').max()
    demand_monthly['price_demand'] = price_demand
    demand_monthly['demand_charge_before'] = demand_monthly.load * price_demand
    demand_monthly['demand_charge_after'] = demand_monthly.net_no_exports * price_demand
    demand_monthly['demand_charge_saving'] = demand_monthly.demand_charge_before - demand_monthly.demand_charge_after
    return demand_monthly


def add_market_charge(energy_table, price_anci, price_market, lf_market):  # ancillary + market fee
    price_market = price_anci + price_market
    energy_table['market_charge_before'] = energy_table.load * price_market * lf_market
    energy_table['market_charge_after'] = energy_table.net_no_exports * price_market * lf_market
    energy_table['market_charge_saving'] = energy_table.market_charge_before - energy_table.market_charge_after
    return energy_table


if __name__ == "__main__":
    # calc load

    raw_load = pd.read_csv("load_15min.csv", index_col=1, parse_dates=True, infer_datetime_format=True, dayfirst=True)
    load = raw_load.drop(['Meter', 'ReadingDate', 'KW', 'KVA', 'PF', 'Quality', 'QualityText'], 1)
    load.columns = ['load']
    load.sort_index(inplace=True)
    indexToDrop = load.index[0:96].append(load.index[-192:])
    load.drop(indexToDrop, inplace=True)
    load = load.resample('1H').sum()
    plt.plot(load)
    plt.show()

    # get billing costs and savings
    pv_cap = 500
    bat_cap = 500
    pv = get_pv_from_PVWatts(f'pv_{pv_cap}.csv', load, pv_cap, 1340)
    energy_table = get_energy_table(load, pv, 7, 23)
    billing_table = add_retail_charge(energy_table, 0.089454, 0.056275, 1.06837)
    billing_table = add_fit(billing_table, 0.08)
    billing_table = add_environ(billing_table, 0.0202, 1.0996)
    billing_table = add_network_energy_charge(billing_table, 0.0467, 0.0247)
    demand_monthly = get_demand_charge_table(energy_table, 119.99/12)
    supply_charge = 8800  # 연간 supply charge
    billing_table = add_market_charge(billing_table, 0.001529, 0.000368, 1.0996)
    market_fee = 0.003633 * 365
    metering_charge = 1250 * 2

    # get years saving
    annual_charge_before = billing_table.retail_charge_before.sum() + billing_table.environ_charge_before.sum()\
                           + billing_table.network_energy_charge_before.sum()\
                           + billing_table.market_charge_before.sum() + demand_monthly.demand_charge_before.sum()
    annual_charge_saving = billing_table.retail_charge_saving.sum() + billing_table.environ_charge_saving.sum() + billing_table.network_energy_charge_saving.sum() + billing_table.market_charge_saving.sum() + demand_monthly.demand_charge_saving.sum() + billing_table.fit_revenue.sum()
    bat_capa = energy_table.net_off_peak_posi.resample('1D').sum().mean()
    pcs_capa = bat_capa / energy_table.net_off_peak_posi.count() / 365 / 2
    bat_revenue = bat_capa * 365 * ((0.089454 - 0.056275) * 1.06837 + (0.04067 - 0.0247))

    # drawings
    load_pv = energy_table.loc[:, ['load', 'pv']]
    plt.plot(load_pv)
    plt.show()

    demand_monthly.loc[:, 'load'].plot(kind='bar')
    plt.show()



