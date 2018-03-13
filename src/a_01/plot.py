#! /usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

def runtime_plot(df, axis):
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
    axis.scatter(type0["size"], type0["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(type1["size"], type1["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(type2["size"], type2["runtime"], marker=markersymbol, s=markersize)
    axis.set_xlabel('Number of processes')
    axis.set_ylabel('Time in seconds')
    axis.legend(['type 0', 'type 1', 'type 2'])
    axis.set_title('Runtime plot')

def speedup_plot(df, axis):
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
    axis.scatter(type0["size"], type0_serial / type0["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(type1["size"], type1_serial / type1["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(type2["size"], type2_serial / type2["runtime"], marker=markersymbol, s=markersize)
    axis.set_xlabel('Number of processes')
    axis.set_ylabel('Speedup')
    axis.legend(['type 0', 'type 1', 'type 2'])
    axis.set_title('Speedup plot')

def distribution_plot(df, axis):
    # select data for run on 48 nodes
    data = df[df["size"] == 48]

    # select data for type 0
    type0 = data[data["type"] == 0]
    type0 = type0.sort_values("rank")

    # select data for type 1
    type1 = data[data["type"] == 1]
    type1 = type1.sort_values("rank")

    # select data for type 2
    type2 = data[data["type"] == 2]
    type2 = type2.sort_values("rank")

    # create plot
    axis.plot(type0["rank"], type0["calctime"] + type0["mpitime"], "-o")
    axis.plot(type1["rank"], type1["calctime"] + type1["mpitime"], "-o")
    axis.plot(type2["rank"], type2["calctime"] + type2["mpitime"], "-o")
    axis.set_xlabel('Process ID')
    axis.set_ylabel('Calctime + MPITime in seconds')
    axis.legend(['type 0', 'type 1', 'type 2'])
    axis.set_title('Time plot')


if __name__ == '__main__':
    # read the .csv file
    df = pd.read_csv('plot_data.csv')
    
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    runtime_plot(df, ax1)

    ax2 = plt.subplot2grid((2, 2), (0, 1))
    speedup_plot(df, ax2)

    ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
    distribution_plot(df, ax3)

    plt.tight_layout()
    plt.show()
    
