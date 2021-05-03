import pandas as pd
from NHCH.utils import plot_Series, plot_DataFrame
from NHCH.read_data import read_load, read_pv_PVWatts
import matplotlib.pyplot as plt


# NHCH lat lon : 20.01, -155.66

# roi
def roi(pv_cap, total_cost, energy_charge, period, derating_ratio, load=read_load('../data/NHCH_1H.csv')):
    """
    economic feasibility study
    :param derating_ratio: unit [%], type [float]
    :param period: unit [year], type [int]
    :param pv_cap: unit [kW], type [int]
    :param total_cost: unit [$], type [int]
    :param energy_charge: unit [$], type [float]
    :param load: unit [kW], type [df], len [8760]
    :return: total profit
    """
    total_revenue = 0
    total_revenue_noReverse = 0
    pv = read_pv_PVWatts(pv_cap)
    for y in range(period):
        pv_thisyear = pv * (1 - derating_ratio) ** y
        net_load_temp = load.load - pv_thisyear.pv
        net_load_noreverse_temp = net_load_temp.clip(0)
        annual_revenue = pv_thisyear.sum() * energy_charge
        annual_revenue_noReverse = (load.load.sum() - net_load_noreverse_temp.sum()) * energy_charge
        total_revenue += annual_revenue
        total_revenue_noReverse += annual_revenue_noReverse
    total_profit = total_revenue - total_cost
    total_profit_noReverse = total_revenue_noReverse - total_cost
    annual_revenue = pv.sum() * energy_charge
    return annual_revenue, total_profit, total_profit_noReverse


def resiliency_sizing(pv_range, reserved_day, mode='max', load=read_load('../data/NHCH_1H.csv')):
    for cap in pv_range:
        ess_weekly_charge_frames = []
        fuel_weekly_frames = []
        for starting_day in range(0, 7):
            merged = load.copy()
            fuel = fuel_calc(merged, 750, 16.3, 19, 21.8, 24.55, 27.4, 30.2, 33.1, 36.1, 39.3, 42.6, 46, 49.6)
            fuel_reserved = int(fuel.median() * 24 * reserved_day)
            if starting_day == 0:
                merged.drop(merged['2017-06-30'].index, inplace=True)
                fuel.drop(fuel['2017-06-30'].index, inplace=True)
            elif starting_day == 1:
                merged.drop(merged['2016-07-01'].index, inplace=True)
                fuel.drop(fuel['2016-07-01'].index, inplace=True)
            else:
                merged = merged.iloc[starting_day * 24: 8760 - (8 - starting_day) * 24]
                fuel = fuel.iloc[starting_day * 24: 8760 - (8 - starting_day) * 24]

            merged['pv'] = read_pv_PVWatts(cap)
            merged['net_load'] = merged.load - merged.pv
            merged['net_load_noReverse'] = merged.net_load.clip(0)
            merged['net_load_overMinPower'] = merged.net_load.clip(200)
            merged['pv_curtailment'] = (200 - merged.net_load).clip(0)
            merged['ess'] = merged.net_load
            merged.loc[merged['pv_curtailment'] > 0]['ess'] = merged[merged['pv_curtailment'] > 0].pv_curtailment
            merged['ess_cumsum_by_week'] = merged.groupby(merged.index.week)['ess'].cumsum()
            ess_daily_charge = merged.pv_curtailment.resample('1D').sum()
            ess_weekly_charge = ess_daily_charge.resample('7D').max()


            merged_weekly = merged.resample('7D').sum()
            fuel_weekly = fuel.resample('7D').sum()
            fuel_weekly['fuel_required'] = fuel_weekly.fuel - fuel_reserved
            fuel_weekly['pv_curtailment_converted_to_fuel'] = merged_weekly.pv_curtailment / 13
            fuel_weekly['fuel_remained'] = fuel_weekly.fuel - fuel_reserved - merged_weekly.pv_curtailment / 13
            ess_weekly_charge_frames.append(ess_weekly_charge)
            fuel_weekly_frames.append(fuel_weekly)

        ess_weekly_charge_concat = pd.concat(ess_weekly_charge_frames).sort_index()
        fuel_weekly_concat = pd.concat(fuel_weekly_frames).sort_index()

        print(f'PV capacity: {cap}\n'
              f'Fuel_remained_max: {fuel_weekly_concat.fuel_remained.max()}\n'
              f'Fuel_remained_median: {fuel_weekly_concat.fuel_remained.median()}')

        if mode == 'max':
            if (fuel_weekly_concat.fuel_remained <= 0).all():

                print(f'Minimum PV capacity required: {cap}kW\n'
                      f'Minimum PCS capacity required: {merged.pv_curtailment.max()}\n'
                      f'Minimum Batt capacity required: {ess_weekly_charge_concat.max()}kWh\n')
                fuel_weekly_concat.plot()
                plt.show()
                return merged, ess_weekly_charge_concat, fuel_weekly_concat
                break
        elif mode == 'median':
            if fuel_weekly_concat.fuel_remained.median() <= 0:
                fuel_weekly_sorted = fuel_weekly_concat.sort_values(by='fuel_remained')
                median_index = fuel_weekly_sorted[fuel_weekly_sorted['fuel_remained'] >
                                                  fuel_weekly_sorted['fuel_remained'].median()].iloc[0].name
                pcs_cap = merged.pv[f'{median_index}':f'{median_index + pd.Timedelta("7 days")}'].max()
                batt_cap = ess_weekly_charge_concat[f'{median_index}']
                # diesel_generation =
                print(f'Minimum PV capacity required: {cap}kW\n'
                      f'Minimum PCS capacity required: {pcs_cap}\n'
                      f'Minimum Batt capacity required: {batt_cap}kWh\n'
                      f'median_index: {median_index}')
                fuel_weekly_concat.plot()
                plt.show()
                return merged, ess_weekly_charge_concat, fuel_weekly_concat
                break
    print('Can not find proper capacity')


def fuel_calc(load, diesel_rated, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12):
    fuel_temp = []
    for l in load.load.values.squeeze():
        if 4 / 16 * diesel_rated <= l < 5 / 16 * diesel_rated:
            fuel_temp.append(c1)
        elif 5 / 16 * diesel_rated <= l < 6 / 16 * diesel_rated:
            fuel_temp.append(c2)
        elif 6 / 16 * diesel_rated <= l < 7 / 16 * diesel_rated:
            fuel_temp.append(c3)
        elif 7 / 16 * diesel_rated <= l < 8 / 16 * diesel_rated:
            fuel_temp.append(c4)
        elif 8 / 16 * diesel_rated <= l < 9 / 16 * diesel_rated:
            fuel_temp.append(c5)
        elif 9 / 16 * diesel_rated <= l < 10 / 16 * diesel_rated:
            fuel_temp.append(c6)
        elif 10 / 16 * diesel_rated <= l < 11 / 16 * diesel_rated:
            fuel_temp.append(c7)
        elif 11 / 16 * diesel_rated <= l < 12 / 16 * diesel_rated:
            fuel_temp.append(c8)
        elif 12 / 16 * diesel_rated <= l < 13 / 16 * diesel_rated:
            fuel_temp.append(c9)
        elif 13 / 16 * diesel_rated <= l < 14 / 16 * diesel_rated:
            fuel_temp.append(c10)
        elif 14 / 16 * diesel_rated <= l < 15 / 16 * diesel_rated:
            fuel_temp.append(c11)
        elif 15 / 16 * diesel_rated <= l < diesel_rated:
            fuel_temp.append(c12)
        else:
            fuel_temp.append(0)
    fuel = pd.DataFrame(fuel_temp, load.index, ['fuel'])
    return fuel


if __name__ == '__main__':
    pv_range = range(1000, 3600, 50)
    merged, ess_weekly_charge_concat, fuel_weekly_concat =\
        resiliency_sizing(pv_range, 3, 'median')
    cdf_index = range(359)
    ess_sorted = ess_weekly_charge_concat.sort_values(ascending=False)
    ess_sorted.index = cdf_index
    ess_sorted.plot()
    plt.show()
    fuel_sorted = fuel_weekly_concat.fuel_remained.sort_values(ascending=False)
    fuel_sorted.index = cdf_index
    fuel_sorted.plot()
    plt.show()

    # merged, ess_weekly, fuel_weekly = resiliency_sizing(pv_range, 3, day, 'median')

    # batt_cap_cdf = abs(merged.ess_cumsum_by_week).resample('7D').max()
    # batt_cap_cdf.name = 'batt_cap_cdf'
    # batt_cap_cdf.sort_values(ascending=False, inplace=True)
    # index = range(0, 52)
    # batt_cap_cdf.index = index
    # batt_cap_cdf.plot()
    # plt.show()
