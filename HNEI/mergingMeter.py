import numpy as np
import pandas as pd
from GAIA.roi_v1 import getPvFromPvWatts
from GAIA.roi_v1 import getEnergyTable
import matplotlib.pyplot as plt

load_FC_raw = pd.read_csv('FC_15min_Load.csv', index_col=1, parse_dates=True, infer_datetime_format=True, dayfirst=True)
load_FC = load_FC_raw.iloc[:,1].resample('1H').sum()['2018']

for capacity in range(350, 400, 50):
    pv_cap = capacity
    print(pv_cap)
    load_RC_raw = pd.read_csv('RC_15min_Load.csv', index_col=1, parse_dates=True, infer_datetime_format=True, dayfirst=True)
    load_RC = load_RC_raw.iloc[:,1].resample('1H').sum()['2018'].to_frame('load')
    pv = getPvFromPvWatts(f'pv_{pv_cap}.csv', load_RC)

    energyTable = getEnergyTable(load_RC, pv, 0, 0, None, 10, 11)
    net = energyTable['net']
    netNoExport = energyTable['netNoExport']

    energyCharge_FC = load_FC['2018-01'].sum() * 0.2
    energyCharge_RC = netNoExport['2018-01'].sum() * 0.2
    demandCharge_FC = load_FC['2018-01'].max() * 10.25
    demandCharge_RC = netNoExport['2018-01'].max() * 10.25

    load_merged = load_FC + net
    load_merged_noexport = load_merged.clip(0)
    load_merged_export = -load_merged.clip(upper=0)
    energyCharge_merged = load_merged_noexport['2018-01'].sum() * 0.2
    demandCharge_merged = load_merged_noexport['2018-01'].max() * 10.25
    print(f'reverse flow gap = {energyTable.gridExport["2018-01"].sum() - load_merged_export["2018-01"].sum()}')
    print(f'energyCharge saving = {energyCharge_RC + energyCharge_FC - energyCharge_merged}')
    print(f'demandCharge saving = {demandCharge_RC + demandCharge_FC - demandCharge_merged}')
