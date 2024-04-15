#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:38:03 2024

@author: nirb

This function plot how variable x change in time using probe data
"""
import matplotlib.pyplot as plt
import math

def plotProbe(file, times, vectors, el,pos, maxpos):
    
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
    plt.subplot(math.floor(maxpos**.5), math.ceil(maxpos**.5), pos)
    
    for i in range(len(el)):
        plt.plot(times[::jp],ele[i][::jp], label=str(i))
    plt.legend()
    plt.title(simulation)
#    plt.show()

    return

def getProbe(file):
    
    f = open(file+"U","r")
    lines = f.readlines()
    
    times = []
    vectors = []
    
    probes=0
    while (lines[probes][:5]=="# Pro"):
        probes+=1
    
    for i in range(len(lines)):
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
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/0/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/500000/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2/postProcessing/probes/247000/"
    # dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/500000/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306122200a2b/postProcessing/probes/0/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/lestest/postProcessing/probes/1000/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2/postProcessing/probes/0/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simple/postProcessing/probes/0/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simple1/postProcessing/probes/0/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2simple/postProcessing/probes/0/"
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simplea/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simplea1/postProcessing/probes/0/")
#    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb1/postProcessing/probes/0/")
    dirs.append("/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800a2simpleb2/postProcessing/probes/0/")

    el = [0,18,33, 45, 69, 75]  #a2
#    el = [0,18,33, 48 , 69, 81]  # b2


    for i in range(len(dirs)):
        times, vectors = getProbe(dirs[i])        
        plotProbe(dirs[i], times, vectors, el,i+1,len(dirs))
        
    plt.show()
