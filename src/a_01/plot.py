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
    markersize = 15
    markersymbol = "x"
    plt.scatter(type0["size"], type0["runtime"], marker=markersymbol, s=markersize)
    plt.scatter(type1["size"], type1["runtime"], marker=markersymbol, s=markersize)
    plt.scatter(type2["size"], type2["runtime"], marker=markersymbol, s=markersize)
    plt.xlabel('Number of processes')
    plt.ylabel('Time in seconds')
    plt.legend(['type 0', 'type 1', 'type 2'])
    plt.title('Runtime plot')

def speedup_plot(df):
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

    # get the serial time for each type
    type0_serial = list(type0[type0["size"] == 1]["runtime"])[0]
    type1_serial = list(type1[type1["size"] == 1]["runtime"])[0]
    type2_serial = list(type2[type2["size"] == 1]["runtime"])[0]

    # create plot
    markersize = 15
    markersymbol = "x"
    plt.scatter(type0["size"], type0_serial / type0["runtime"], marker=markersymbol, s=markersize)
    plt.scatter(type1["size"], type1_serial / type1["runtime"], marker=markersymbol, s=markersize)
    plt.scatter(type2["size"], type2_serial / type2["runtime"], marker=markersymbol, s=markersize)
    plt.xlabel('Number of processes')
    plt.ylabel('Speedup')
    plt.legend(['type 0', 'type 1', 'type 2'])
    plt.title('Speedup plot')

if __name__ == '__main__':
    # read the .csv file
    df = pd.read_csv('plot_data.csv')
    
    plt.subplot(2, 1, 1)
    runtime_plot(df)

    plt.subplot(2, 1, 2)
    speedup_plot(df)

    plt.tight_layout()
    plt.show()
    
