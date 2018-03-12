#! /usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

def runtime_plot(df):
    # select data of master process
    masters_filter = df["rank"] == 0
    masters = df[masters_filter]
    
    # select data for type 0
    type0_filter = masters["type"] == 0
    type0 = masters[type0_filter]

    # select data for type 1
    type1_filter = masters["type"] == 1
    type1 = masters[type1_filter]

    # select data for type 2
    type2_filter = masters["type"] == 2
    type2 = masters[type2_filter]

    # create plot
    plt.scatter(type0["size"], type0["runtime"], marker=".")
    plt.scatter(type1["size"], type1["runtime"], marker="+")
    plt.xlabel('Number of processes')
    plt.ylabel('Time in seconds')
    plt.legend(['type 0', 'type 1'])
    plt.title('Runtime plot')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # read the .csv file
    df = pd.read_csv('plot_data.csv')

    runtime_plot(df)
    
