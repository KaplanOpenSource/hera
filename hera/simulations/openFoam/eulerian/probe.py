#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:38:03 2024

@author: nirb

This function plots how variable x changes in time using probe data
"""
import matplotlib.pyplot as plt

def plotProbe(file, times, vectors, el):
    
    dirr = file.replace("/"," ")
    simulation = dirr.split(" ")[len(dirr.split(" "))-5]

    ele=[]
    for i in range(len(el)):
        ele.append([float(elem[el[i]]) for elem in vectors])
    
    jp=10
    plt.figure() #(figsize=(10,5))
    for i in range(len(el)):
        plt.plot(times[::jp],ele[i][::jp], label=str(i))
    plt.legend()
    plt.title(simulation)
    plt.show()
    
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
    #for i in range(100):
        #start from data
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

    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/0/"
    dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800c2/postProcessing/probes/100000/"
    # dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306121800b2/postProcessing/probes/500000/"
    # dir = "/data5/NOBACKUP/nirb/Simulations/Haifa/wrf202306122200a2b/postProcessing/probes/0/"

    times, vectors = getProbe(dir)        

    el = [0,18,33, 45, 69, 75]

    plotProbe(dir, times, vectors, el)

