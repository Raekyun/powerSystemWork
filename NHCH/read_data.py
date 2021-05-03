import pandas as pd

def read_load(path):
    load = pd.read_csv(path, index_col=0, parse_dates=True)
    load = load.drop(load.iloc[8760].name)
    load.columns = ['load']
    load.index = load.index.tz_localize('US/Hawaii')
    return load


# PVWatts PV data csv to dataframe
def read_pv_PVWatts(cap, load=read_load('../data/NHCH_1H.csv')):
    """
    load의 index와 일치도록 만듦
    #TODO: PV index를 계절, 월, 일이 맞도록 재조정해야함
    :param load: unit [kW], type [df], len [8760]
    :param cap: unit [kW], type [int], pv raw data가 data 폴더에 포맷에 맞춰서 있어야 함
    :return: pv
    """
    pv = pd.read_csv(f'../data/pvwatts_hourly_{cap}.csv', skiprows=16)
    pv = pv['AC System Output (W)']
    pv.name = 'pv'
    pv = pv.drop(8760)
    pv.index = load.index
    pv = pv.divide(1000)
    pv = pd.DataFrame(pv)
    return pv