import pandas as pd
import matplotlib.pyplot as plt


def plot_Series(series, name):
    series.name = name
    df = pd.DataFrame(series)
    df.plot()
    plt.show()


def plot_DataFrame(df):
    df.plot()
    plt.show()
