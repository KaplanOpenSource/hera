#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:38:03 2024

@author: nirb

This function plots how variable x changes in time using probe data

"""
import matplotlib.pyplot as plt

def plotProbe(file, times, vectors, el,pos, maxpos):
    # plot the openFOAM probes data
    
    # file - for the title
    # times, vectors - the data to plot
    # el - which probes to plot
    # We need pos and max pos to know where to plot this subfigure in the figure
    
    dirr = file.replace("/"," ")
    simulation = dirr.split(" ")[len(dirr.split(" "))-5]

    ele=[]
    for i in range(len(el)):
        ele.append([float(elem[el[i]]) for elem in vectors])
    
    jplen=len(times)
    if jplen>100000:
       jp=jplen//10000
    else:
       jp=10    
    #plt.figure() #(figsize=(10,5))
    col=int(maxpos**.5)
    row=col
    if row*row<maxpos:
       col +=1
    if row*row<maxpos:
       row +=1
    plt.subplot(row, col, pos)
    
    for i in range(len(el)):
        plt.plot(times[::jp],ele[i][::jp], label=str(i))
    plt.legend()
    plt.title(simulation)
#    plt.show()

    return

def getProbe(file, jumps = 1):
    # get openfoam probe data and convert it to python data
    f = open(file+"U","r")
    lines = f.readlines()
    times = []
    vectors = []
    probes=0
    while (lines[probes][:5]=="# Pro"):
        probes+=1
    if len(lines)>100000:
        jumps = len(lines)//100000
    for i in range(len(lines)):
        if i % jumps ==0:
            if(i>probes):
                #first 20 chars of line are times
                times.append(float(lines[i][0:14]))
                #remaining chars are vector
                vector = lines[i][14:]
                #remove ( and )
                vector = vector.replace("\n", "")
                vector = vector.replace("(", "")
                vector = vector.replace(")", "")
                #split by space
                vectors.append(vector.split(" "))
    return times, vectors

if __name__ == "__main__":

    dirs=[]
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/0/")
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/500000/")
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2/postProcessing/probes/247000/")
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcesddsing/probes/500000/")
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306122200a2b/postProcessing/probes/0/")
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/lestest/postProcessing/probes/1000/")
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/postProcessing/probes/0/")
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simple/postProcessing/probes/0/")
    #dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simple1/postProcessing/probes/0/")
    
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2simple/postProcessing/probes/0/")

#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simple/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simplea1/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb1/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb2/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2lmnr/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2lmnr2/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2lmnr1/postProcessing/probes/0/")

#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/postProcessing/probes/0/")
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800d2/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a21/postProcessing/probes/0/")

    el = [0,18,33, 45, 69, 75]  #a2
    el = [0,18,33, 48 , 69, 81]  # b2
    el = [3,21,36, 48 , 72, 81]  # c2


    for i in range(len(dirs)):
        print('probing ',i,' out of ',len(dirs))
        times, vectors = getProbe(dirs[i])
        plotProbe(dirs[i], times, vectors, el,i+1,len(dirs))
        
    plt.show()
    
