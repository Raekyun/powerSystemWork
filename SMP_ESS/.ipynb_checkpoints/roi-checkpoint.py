import pandas as pd
import numpy as np
from scipy.optimize import linprog

def dailyEssOpt(dailySmp, pcsCap, battCap):
    c = -dailySmp
    a_eq = np.ones([1, 24])
    b_eq = battCap
    x_bounds = []
    for k in range(24):
        x_bounds.append((0, pcsCap))
    res = linprog(c=c, A_eq=a_eq, b_eq=b_eq, bounds=x_bounds)
    return res.x

energyTable = pd.read_csv('smp_jeju_one_week.csv', index_col=0, parse_dates=True)

pcsCap = 500
battCap = 2000

energyTable['smp'] = [energyTable['smp'][i] if (17 <= i.hour or i.hour==0) else 0 for i in energyTable.index]
energyTable['ess_uniform'] = [pcsCap if 17 <= h < 17 + battCap/pcsCap else 0 for h in energyTable.index.hour]
energyTable['revenue_uniform'] = energyTable['smp'] * energyTable['ess_uniform']

ess = np.array([])
for k in range(7):
    dailySmp = energyTable['smp'].iloc[k*24:(k+1)*24]
    dailyEss = dailyEssOpt(dailySmp, pcsCap, battCap)
    ess = np.append(ess, dailyEss)

energyTable['ess_opt'] = ess
energyTable['revenue_opt'] = energyTable['smp'] * energyTable['ess_opt']

