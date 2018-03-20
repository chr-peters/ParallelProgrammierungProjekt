#! /usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

def runtime_plot(df, df_mpi, axis):
    # create plot
    markersize = 15
    markersymbol = "x"
    axis.scatter(df[df["type"]=="static"]["size"], df[df["type"]=="static"]["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(df[df["type"]=="static-1"]["size"], df[df["type"]=="static-1"]["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(df[df["type"]=="dynamic"]["size"], df[df["type"]=="dynamic"]["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(df[df["type"]=="guided"]["size"], df[df["type"]=="guided"]["runtime"], marker=markersymbol, s=markersize)
    mpi_masters = df_mpi[df_mpi["rank"]==0]
    axis.scatter(mpi_masters[mpi_masters["type"]==2]["size"][:24], mpi_masters[mpi_masters["type"]==2]["runtime"][:24], marker=markersymbol, s=markersize)
    axis.set_xlabel('Number of processes')
    axis.set_ylabel('Time in seconds')
    axis.legend(['static', 'static-1', 'dynamic', 'guided', 'mpi master-slave'])
    axis.set_title('Runtime plot')

def speedup_plot(df, df_mpi, axis):
    static_serial = list(df[df["type"]=="static"]["runtime"])[0]
    static1_serial = list(df[df["type"]=="static-1"]["runtime"])[0]
    dynamic_serial = list(df[df["type"]=="dynamic"]["runtime"])[0]
    guided_serial = list(df[df["type"]=="guided"]["runtime"])[0]

    # create plot
    markersize = 15
    markersymbol = "x"
    axis.scatter(df[df["type"]=="static"]["size"], static_serial / df[df["type"]=="static"]["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(df[df["type"]=="static-1"]["size"], static1_serial / df[df["type"]=="static-1"]["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(df[df["type"]=="dynamic"]["size"], dynamic_serial / df[df["type"]=="dynamic"]["runtime"], marker=markersymbol, s=markersize)
    axis.scatter(df[df["type"]=="guided"]["size"], guided_serial / df[df["type"]=="guided"]["runtime"], marker=markersymbol, s=markersize)
    mpi_masters = df_mpi[df_mpi["rank"]==0]
    mpi_serial = list(mpi_masters[mpi_masters["type"]==2]["runtime"])[0]
    axis.scatter(mpi_masters[mpi_masters["type"]==2]["size"][:24], mpi_serial / mpi_masters[mpi_masters["type"]==2]["runtime"][:24], marker=markersymbol, s=markersize)
    axis.set_xlabel('Number of processes')
    axis.set_ylabel('Speedup')
    axis.legend(['static', 'static-1', 'dynamic', 'guided', 'mpi master-slave'])
    axis.set_title('Speedup plot')

if __name__ == '__main__':
    # read the .csv file
    df = pd.read_csv('plot_data.csv')
    df_mpi = pd.read_csv('../a_01/plot_data.csv')
    
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    runtime_plot(df, df_mpi, ax1)

    ax2 = plt.subplot2grid((1, 2), (0, 1))
    speedup_plot(df, df_mpi, ax2)

    plt.tight_layout()
    plt.show()
    
