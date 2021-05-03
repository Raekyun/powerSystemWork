# 오전 시간대도 최적화 범위에 포함

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

def preprocessEpsisData(df):
    columns = []
    for k in range(1, 28):
        columns.append(f'{k}')
    df.columns = columns
    df = df.drop(columns=df.columns[-3:])
    df = df.sort_index()
    index = pd.date_range(df.index[0], df.index[-1], freq='H')
    index = pd.date_range(index[0] + index[0].freq, index[-1] + index[-1].freq * 24, freq='H')
    data = df.values
    data = data.flatten()
    columns = ['smp']
    df = pd.DataFrame(data, index, columns)
    return df

df = pd.read_csv('smp_land_2020_2021.csv', index_col=0, parse_dates=True)
energyTable = preprocessEpsisData(df)

pcsCap = 500
battCap = 2000

energyTable['smp_for_calc'] = [energyTable['smp'][i] if (17 <= i.hour or i.hour <=10) else 0 for i in energyTable.index]
energyTable['ess_uniform'] = [pcsCap if 17 <= h < 17 + battCap/pcsCap else 0 for h in energyTable.index.hour]
# energyTable['ess_uniform'] = [battCap/18 if (16 <= h or h <= 9) else 0 for h in energyTable.index.hour]
energyTable['revenue_uniform'] = energyTable['smp'] * energyTable['ess_uniform']

analysisPeriod = int(energyTable.index.shape[0]/24)
ess = np.array([])
for k in range(analysisPeriod):
    dailySmp = energyTable['smp_for_calc'].iloc[k*24:(k+1)*24]
    dailyEss = dailyEssOpt(dailySmp, pcsCap, battCap)
    ess = np.append(ess, dailyEss)

energyTable['ess_opt'] = ess
energyTable['revenue_opt'] = energyTable['smp'] * energyTable['ess_opt']

print(f'Total revenue with uniform ESS: {energyTable["revenue_uniform"].sum()}')
print(f'Total revenue with optimized ESS: {energyTable["revenue_opt"].sum()}')

date = '2020-01-01'
right_ax_ylim = energyTable.loc[:, 'smp'][date].max() + 5
ax = energyTable[['ess_uniform', 'ess_opt', 'smp']].loc[date].plot(secondary_y='smp', ylim=[0, 1500])
ax.right_ax.set_ylim(0, right_ax_ylim)

plt.show()

energyTable.to_csv('energyTable.csv')