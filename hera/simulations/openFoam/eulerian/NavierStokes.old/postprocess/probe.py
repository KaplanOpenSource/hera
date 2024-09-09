#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:38:03 2024

@author: nirb

This function plots how variable x changes in time using probe data

"""
import matplotlib.pyplot as plt
import datetime
import math

def plotProbe(file, times, vectors, el,pos, maxpos, offset=0):
    # plot the openFOAM probes data
    
    # file - for the title
    # times, vectors - the data to plot
    # el - which probes to plot
    # We need pos and max pos to know where to plot this subfigure in the figure
    
    dirr = file.replace("/"," ")
    simulation = dirr.split(" ")[len(dirr.split(" "))-5]

#    print('1))',el)
    elbkup=el.copy()
    for i in range(len(el)):
        while(math.fabs(float(vectors[-1][elbkup[i]]))> 10000):
            elbkup[i]+=3
        elbkup[i]+=3*offset   
#    print('el2))',elbkup)
    
    ele=[]
    for i in range(len(elbkup)):
        ele.append([float(elem[elbkup[i]]) for elem in vectors])
    
    jplen=len(times)
    if jplen>100000:
       jp=jplen//10000
    else:
       jp=10    
    #plt.figure() #(figsize=(10,5))
    col=int(maxpos**.5)
    row=col
    if row*col<maxpos:
       col +=1
    if row*col<maxpos:
       row +=1
    plt.subplot(row, col, pos)
    
    for i in range(len(elbkup)):
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
    
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2simple/postProcessing/probes/0/")

#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simple/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simplea1/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb1/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb2/postProcessing/probes/0/")
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2lmnr/postProcessing/probes/0/")
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2lmnr2/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2lmnr1/postProcessing/probes/0/")

#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/postProcessing/probes/0/")
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2/postProcessing/probes/0/")
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800d2/postProcessing/probes/0/")
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/postProcessing/probes/0/")
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a21/postProcessing/probes/0/")

##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a22/postProcessing/probes/0/")


#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a23/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a24/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a25/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b25/postProcessing/probes/0/")

    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2a/postProcessing/probes/0/")
    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2a/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2a/postProcessing/probes/0/")
##    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2a2/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800d2a/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800d2az2/postProcessing/probes/0/")
#####    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2az1/postProcessing/probes/0/")
    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2az1/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2az0/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2z0simple/postProcessing/probes/0/")
#     dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2az1/postProcessing/probes/0/")
    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2az2/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2az0/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2a1/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2a2/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2aLES/postProcessing/probes/0/")
    
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/syn1/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/syn2/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/syn3/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/syn4/postProcessing/probes/0/")

    dirs=[]
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2az0/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2az0test/postProcessing/probes/0/")
    
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2az0/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200c2az0/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200b2az0simple/postProcessing/probes/0/")
    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2az0/postProcessing/probes/0/")
    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb1/postProcessing/probes/0/")
    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb2/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2simple/postProcessing/probes/0/") # 150k???

    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2az1test/postProcessing/probes/0/") # 150k???
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2a2/postProcessing/probes/0/") # 150k???
    
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a1z0/postProcessing/probes/0/")
    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a1z0t/postProcessing/probes/0/")
    # dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a1z0u/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a2z0/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200a4z0/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306131200cube/postProcessing/probes/0/")


    el = [0,18,33, 45, 69, 75]  #a2
    el = [0,18,33, 45, 60, 75]  #a2
    el = [0,18,33, 48 , 69, 81]  # b2
#    el = [3,21,36, 48 , 72, 81]  # c2
#    el=[0, 18,33,45,72,84] # syn
#    el=[6, 24, 39,51,75,84] # syn
    el = [0,15,30, 45, 60, 75]  #dynamic
#    el = [6, 24, 39, 54, 66, 87]

    print('start @ ',datetime.datetime.now())
    for i in range(len(dirs)):
        times, vectors = getProbe(dirs[i])
        plotProbe(dirs[i], times[1:], vectors[1:], el,i+1,len(dirs), offset=1)
        # plotProbe(dirs[i], times[1:], vectors[1:], el,i+1,len(dirs), offset=1)

#        plotProbe(dirs[i], times, vectors, el,4*i+1,len(dirs)*4, offset=0)
#        plotProbe(dirs[i], times, vectors, el,1+4*i+1,len(dirs)*4, offset=1)
#        plotProbe(dirs[i], times, vectors, el,2+4*i+1,len(dirs)*4, offset=2)
#        plotProbe(dirs[i], times, vectors, el,3+4*i+1,len(dirs)*4, offset=3)
        print('probed ',i,' out of ',len(dirs), '(',len(times),')')        
    plt.show()
    
