from scipy.optimize import linprog
from mip import Model, xsum, maximize, BINARY

p = [20, 40, 20, 15, 30]
w1 = [5, 4, 3, 7, 8]
I1 = range(len(w1))
w2 = [1, 7, 9, 4, 6]
I2 = range(len(w2))
w3 = [8, 10, 2, 1, 10]
I3 = range(len(w3))
c1 = 25
c2 = 25
c3 = 25


m = Model()
x = [m.add_var(var_type=BINARY) for i in I1]
m.objective = maximize(xsum(p[i]*x[i] for i in I1))

m += xsum(w1[i] * x[i] for i in I1) <= c1
m += xsum(w2[i] * x[i] for i in I2) <= c2
m += xsum(w3[i] * x[i] for i in I3) <= c3

selected = [i for i in I1 if x[i].x >= 0.99]