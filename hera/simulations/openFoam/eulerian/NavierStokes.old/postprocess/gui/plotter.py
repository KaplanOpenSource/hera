# cd "../pymeteoGIS"
# %run ofplot2.py

usehera=False
usehera=True
import paraview.simple as pvsimple

import sys
import matplotlib.pyplot as plt

if usehera:
    from hera.simulations.openFoam.datalayer.pvOpenFOAMBase import paraviewOpenFOAM
    from pynumericalmodels.OpenFOAM.postprocess.gui.measurements import chooseobservation
    from pynumericalmodels.OpenFOAM.postprocess.gui.morphology import choosearea, meanheight
    from pynumericalmodels.OpenFOAM.postprocess.gui.ml import ml
else:
#    sys.path.append('pynumericalmodels/OpenFOAM/postprocess/gui')
#    sys.path.append('pynumericalmodels/OpenFOAM/postprocess')
    sys.path.append('..')
#    from pynumericalmodels.OpenFOAM.postprocess.pvOpenFOAMBase import pvOFBase
    from pvOpenFOAMBase import pvOFBase
    from measurements import chooseobservation
    from morphology import choosearea, meanheight
    from ml import ml


from paraview import servermanager


#from scipy.interpolate import griddata
import json
# import io
# sys.path.insert(0,'..')
# from pvOpenFOAMBase import pvOFBase

from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLineEdit, QLabel, QListWidget, QCheckBox, QMessageBox\
    , QComboBox, QFileDialog
# import PyQt5.QtGui
from PyQt5.QtCore import pyqtSlot
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# import mlab
import math


import datetime

import re
from scipy.interpolate import griddata
from scipy import stats
import sklearn.metrics
from sklearn.cross_decomposition import CCA
import random

from threading import Thread


# import psutil
if usehera is False:
    bse = pvOFBase()
axisscale = 1. # 0.001


def shrinkdir(filename):
    lastslash = filename.rfind('/')
    return filename[:lastslash+1]

def shrinktitle(filename):
    lastslash = filename.rfind('/')
    previousslash = filename[:lastslash].rfind('/')
    return filename[previousslash+1:lastslash]

def rmse(a,b, zeros=True):
    if zeros:
        return round(((a - b) ** 2).mean() ** .5,4)
    else:
        b1=b[a!=0]
        a1=a[a!=0]
        a2=a1[b1!=0]
        b2=b1[b1!=0]
        return round(((a2 - b2) ** 2).mean() ** .5,4)
        


def logreader(filename, time1):
    slash = filename.rfind("/")+1
    directory = filename[:slash]

    try:
        logfile1 = directory + "log.simpleFoam"
        fh = open(logfile1, 'r')
    except:
        logfile1 = directory + "log"
        fh = open(logfile1, 'r')

    logtxt = fh.read()
    idx = [m.start() for m in re.finditer('Time', logtxt)]
    time2 = int(float(time1)*3-2)
    try:
        idx2 = logtxt[idx[time2+1]:idx[time2+1]+20]
        result = idx2[7:idx2.find('s') - 1]
    except:
        result = "no record"
    return result


def logall(filename):
    resultitter = []
    resulttime = []
    slash = filename.rfind("/")+1
    directory = filename[:slash]
    try:
        logfile1 = directory + "log.simpleFoam"
        fh = open(logfile1, 'r')
        parallel = False
    except:
        logfile1 = directory + "log"
        fh = open(logfile1, 'r')
        parallel = True

    logtxt = fh.read()
    idx = [m.start() for m in re.finditer('Time', logtxt)]
    idxitter = idx[1::3]
    idxtime = idx[2::3]
    for i in range(len(idxtime)):   # I want to have array of itteration number and array of time
                                    # than I will have array of correlation of current itteration
                                    #  with last one on height as gui
        idx2 = logtxt[idxitter[i]+7:idxitter[i]+15]
        resultitter.append(int(re.findall("\d+", idx2)[0]))
        idx3 = logtxt[idxtime[i]:idxtime[i]+15]
        if len(re.findall("\d+\.\d+", idx3)) == 0:
            resulttime.append(float(re.findall("\d+", idx3)[0]))
        else:
            resulttime.append(float(re.findall("\d+\.\d+", idx3)[0]))
    print ('logall3', len(idxtime))

    return resultitter, resulttime

def stat(a,b, kind='rmse', nantozero=True, verbose=None):
    # return np.sqrt(np.mean((predictions-targets)**2))
    # np.corrcoef(db[i][sel], db[-1][sel])[1, 0]
    try:
        if nantozero:
            a[np.isnan(a)]=0.0
            b[np.isnan(b)]=0.0
        else:
            a=a[~np.isnan(a)]
            b=b[~np.isnan(b)]
            minlen=min(len(a),len(b))
            a=a[:minlen]
            b=b[:minlen]
        if (len(a)<1) or (len(b)<1):
            c=-999
        else:
            if kind=='rmse':
                c = np.sqrt(np.mean((a-b)**2))
            elif kind=='mae':
                c = np.mean(np.abs(a-b))
            elif kind=='bias':
                c = np.mean(a-b)
            elif kind=='r2':
                # coefficient of determination
                c = sklearn.metrics.r2_score(a.ravel(), b.ravel())
            elif kind=='r':
                # Pearson correlation coefficient 
                c = np.corrcoef(a, b)[1, 0]
        #        if len(c.shape)>1:
        #            c=c[1,0]
            elif kind=='nmse':
                # normalized mean square error, chang and hanna 2004
                c = np.sqrt(np.mean((a-b)**2))/(np.mean(a)*np.mean(b))
            elif kind=='fb':
                # Fractional bias, chang and hanna 2004
                c = (np.mean(a)-np.mean(b))/(0.5*(np.mean(a)+np.mean(b)))
            else:
                print('unknown kind ',kind)
                c=-998
    except:
        c=-997
        print('stat except', a, b, len(a), len(b), kind)

    return round(c,4)


def corrall(filename, sel, axisIndex , axisValue, clipxmin, clipxmax, clipymin, clipymax, clipzmin, clipzmax
            , heightcontour, reader=None,end=None, parallel=True, last=None):
    corrlast = []
    corrnear = []
    print(filename, sel)
    if usehera:
        times = get_times(filename)        
        if end is not None:
            times = times[:end]
        if last is not None:
            times=times[-last:]
    else:
        if reader is None:
    #        reader = bse.ReadCase(filename, filename, CaseType='Reconstructed Case')  # ' Reconstructed Decomposed Case')
            if parallel:
                reader = bse.ReadCase(filename, filename, CaseType='Decomposed Case')  # 'Decomposed  Reconstructed Case')
            else:
                reader = bse.ReadCase(filename, filename, CaseType='Reconstructed Case')  # 'Decomposed  Reconstructed Case')            
        if end is None:
            times = reader.TimestepValues[:]
        else:
            times = reader.TimestepValues[:end]
        if last is not None:
            times=times[-last:]
        print('TimestepValues: ', reader.TimestepValues)

    if heightcontour:
        db = get_slice_height(filename, times,  axisIndex, axisValue,
                              clipxmin=clipxmin, clipxmax=clipxmax,
                              clipymin=clipymin, clipymax=clipymax,
                              clipzmin=clipzmin, clipzmax=clipzmax)
    else:
        db = get_slice(filename, times, axisIndex, axisValue,
                       clipxmin=clipxmin, clipxmax=clipxmax,
                       clipymin=clipymin, clipymax=clipymax,
#                       clipzmin=clipzmin, clipzmax=clipzmax, reader=reader)
                       clipzmin=clipzmin, clipzmax=clipzmax, parallel=parallel)

    for i in range(len(times)-1):
        corrlast.append(stat(db[i][sel].values, db[-1][sel].values, kind='r2'))
        corrnear.append(stat(db[i][sel].values, db[i+1][sel].values, kind='r2'))
        print('i is', i, corrlast[-1], corrnear[-1], times[i] )

    return  corrlast, corrnear, times[:-1]
    
    
def makegrid(mat, field, axis=2, ticks=400, method='linear'):
    if axis == 0:  # 'x'
        x1 = mat['y']
        y1 = mat['z']
        z1 = mat[field]
    elif axis == 1:  # 'y'
        x1 = mat['x']
        y1 = mat['z']
        z1 = mat[field]
    elif axis == 2:  # 'z'
        x1 = mat['x']
        y1 = mat['y']
        z1 = mat[field]
    else:
        print('Something wrong')
        raise ValueError('Axis is not 0-2 (X/Y/Z).')

    xmin = min(x1)
    xmax = max(x1)
    ymin = min(y1)
    ymax = max(y1)
    # https://stackoverflow.com/questions/45889329/compromise-between-scipy-interpolate-rbf-and-scipy-interpolate-griddata
    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(xmin, xmax, ticks), np.linspace(ymin, ymax, ticks)
    xi, yi = np.meshgrid(xi, yi)
    # Normalize data and grid.
    x1_new, xi_new = (x1 - xmin) / (xmax - xmin), (xi - xmin) / (xmax - xmin)
    y1_new, yi_new = (y1 - ymin) / (ymax - ymin), (yi - ymin) / (ymax - ymin)

    # Interpolate new data with griddata.
    if usehera:
        # Python 3.6
        vals = np.array([x1_new, y1_new]).T        
    else:
        # Python 2.7
        vals = zip(*[x1_new, y1_new])
#    print('test111')
    print ('method=',method, field)
#    print('vals', vals)
    zi1 = griddata(vals, z1.values, (xi_new, yi_new), method=method)  # 'nearest' , linear
#    print('test2222')
    return xmin, xmax, ymin, ymax, zi1

def makegrid3d(mat, field, ticks=40):
    x1 = mat['x']
    y1 = mat['y']
    z1 = mat['z']
    v1 = mat[field]

    xmin = min(x1)
    xmax = max(x1)
    ymin = min(y1)
    ymax = max(y1)
    zmin = min(z1)
    zmax = max(z1)
    # https://stackoverflow.com/questions/45889329/compromise-between-scipy-interpolate-rbf-and-scipy-interpolate-griddata
    # Set up a regular grid of interpolation points
    xi, yi, zi = np.linspace(xmin, xmax, ticks), np.linspace(ymin, ymax, ticks), np.linspace(zmin, zmax, ticks)
    xi, yi ,zi= np.meshgrid(xi, yi,zi)
    # Normalize data and grid.
    x1_new, xi_new = (x1 - xmin) / (xmax - xmin), (xi - xmin) / (xmax - xmin)
    y1_new, yi_new = (y1 - ymin) / (ymax - ymin), (yi - ymin) / (ymax - ymin)
    z1_new, zi_new = (z1 - zmin) / (zmax - zmin), (zi - zmin) / (zmax - zmin)

    # Interpolate new data with griddata.
    if usehera:
        # Python 3.6
        vals = np.array([x1_new, y1_new, z1_new]).T        
    else:
        # Python 2.7
    #    vals = zip(*[x1_new, y1_new])
        vals = zip(*[x1_new, y1_new, z1_new])
#    print('test111')
    vi1 = griddata(vals, v1, (xi_new, yi_new, zi_new), method='nearest')  # 'nearest' , linear
#    print('test2222')
    return xmin, xmax, ymin, ymax, zmin, zmax, vi1

def get_times(filename):
    bse = paraviewOpenFOAM(shrinkdir(filename))
    times = bse.reader.TimestepValues
    return times

def get_slice(filename, itteration, axis_index, axis_value, offset=0.0
              , clipxmin=0.0, clipxmax=0.0, clipymin=0.0, clipymax=0.0, clipzmin=0.0, clipzmax=0.0, reader=None, parallel=False):

    if usehera:
        # dirname=r'/data4bk/nirb/Simulations/Mala/mala2b/'
        bse = paraviewOpenFOAM(shrinkdir(filename))
        # bse = paraviewOpenFOAM(dirname)
        reader=bse.reader
        slice_1 = pvsimple.Slice(Input = bse.reader,guiName="mySlice")
    else:
        if reader is None:
            if parallel:
                if usehera is False:
                    reader = bse.ReadCase(filename, filename, CaseType='Decomposed Case')  # ' Reconstructed Decomposed Case')
            else:
                if usehera is False:
                    reader = bse.ReadCase(filename, filename, CaseType='Reconstructed Case')  # ' Reconstructed Decomposed Case')
    #        print('Initial timelist print:', len(reader.TimestepValues), reader.TimestepValues[-1])
            print ('first time takes time')
        else:
            print ('save time and memory')
        slice_1 = pvsimple.Slice(reader)
    # textboxaxisValue = float(self.axispos.text())
    slice_1.SliceType = 'Plane'
    slice_1.SliceOffsetValues = offset #[0]
    normal = [0.0, 0.0, 0.0]
    normal[axis_index] = 1.0  # int(self.listaxis.currentRow())
    slice_1.SliceType.Origin = [x * axis_value for x in normal]  # normal*textboxaxisValue # [0.0, 0.0, 2.0]
    slice_1.SliceType.Normal = normal  # [0.0, 0.0, 1.0]
    slice_1.UpdatePipeline()

    if clipxmin < clipxmax:
        # if axis_index == 0:
        #     raise ValueError('Tried to set boundaries in the main axis')
        clip_xmin = pvsimple.Clip(slice_1)
        clip_xmin.ClipType = "Plane"
        clip_xmin.ClipType.Normal = [1.0, 0.0, 0.0]
        clip_xmin.ClipType.Origin = [clipxmin, 0.0, 0.0]
        if usehera:
            clip_xmin.Invert = 0
        else:
            clip_xmin.InsideOut = 0
        clip_xmin.UpdatePipeline()
        clip_xmax = pvsimple.Clip(clip_xmin)
        clip_xmax.ClipType = "Plane"
        clip_xmax.ClipType.Normal = [1.0, 0.0, 0.0]
        clip_xmax.ClipType.Origin = [clipxmax, 0.0, 0.0]
        if usehera:
            clip_xmin.Invert = 1
        else:
            clip_xmin.InsideOut = 1
        clip_xmax.UpdatePipeline()
        slice_1 = clip_xmax
    if clipymin < clipymax:
        # if axis_index == 1:
            # raise ValueError('Tried to set boundaries in the main axis')
        clip_ymin = pvsimple.Clip(slice_1)
        clip_ymin.ClipType = "Plane"
        clip_ymin.ClipType.Normal = [0.0, 1.0, 0.0]
        clip_ymin.ClipType.Origin = [0.0, clipymin, 0.0]
        if usehera:
            clip_ymin.Invert = 0
        else:
            clip_ymin.InsideOut = 0
        clip_ymin.UpdatePipeline()
        clip_ymax = pvsimple.Clip(clip_ymin)
        clip_ymax.ClipType = "Plane"
        clip_ymax.ClipType.Normal = [0.0, 1.0, 0.0]
        clip_ymax.ClipType.Origin = [0.0, clipymax, 0.0]
        if usehera:
            clip_ymin.Invert = 1
        else:
            clip_ymin.InsideOut = 1
        clip_ymax.UpdatePipeline()
        slice_1 = clip_ymax
    if clipzmin < clipzmax and axis_index != 2:
        # if axis_index == 2:
            # raise ValueError('Tried to set boundaries in the main axis')
        clip_zmin = pvsimple.Clip(slice_1)
        clip_zmin.ClipType = "Plane"
        clip_zmin.ClipType.Normal = [0.0, 0.0, 1.0]
        clip_zmin.ClipType.Origin = [0.0, 0.0, clipzmin]
        if usehera:
            clip_zmin.Invert = 0
        else:
            clip_zmin.InsideOut = 0
        clip_zmin.UpdatePipeline()
        clip_zmax = pvsimple.Clip(clip_zmin)
        clip_zmax.ClipType = "Plane"
        clip_zmax.ClipType.Normal = [0.0, 0.0, 1.0]
        clip_zmax.ClipType.Origin = [0.0, 0.0, clipzmax]
        if usehera:
            clip_zmin.Invert = 1
        else:
            clip_zmin.InsideOut = 1
        clip_zmax.UpdatePipeline()
        slice_1 = clip_zmax

    if usehera:
        x=bse.to_pandas("mySlice",timelist=itteration)
        db=[]
        for z in x:
            # y=next(z)
            # print('z',z,x, itteration, parallel,filename)
            db.append(z['mySlice'])
        if len(db)==1:
            db=db[0]
        return db        
        # for x in  bse.to_pandas("mySlice",timelist=itteration):
            # print(x)        
    else:
        itr = bse.to_pandas(reader, slice_1, timelist=itteration)
        del reader
        print('type list111', itteration)
        if isinstance(itteration, list):
            print('type list', len(itteration))
            db = []
            for i in range(len(itteration)):
                dbb = itr.next()
                db.append(dbb[dbb.keys()[0]])
                if len(dbb[dbb.keys()[0]]) == 0:
                    print('There are no records in get slice 1, please check coordinates limits or slice location')
            return db
        else:
            dbb = itr.next()
            db = dbb[dbb.keys()[0]]
            if len(db) == 0:
                print('There are no records in get slice 2, please check coordinates limits or slice location')
            return db


def get_slice_height(filename, itteration, axis_index, axis_value, offset=0.0,
                     clipxmin=0.0, clipxmax=0.0, clipymin=0.0, clipymax=0.0, clipzmin=0.0, clipzmax=0.0, reader = None):
    if reader is None:
        reader = bse.ReadCase('case name', filename, CaseType='Decomposed Case')  # 'Reconstructed Case')
        print('Initial timelist print:', reader.TimestepValues)

    contour1 = pvsimple.Contour(reader)
    contour1.ContourBy = ['POINTS', 'HeightFromTopo']
    try:
        realoffset = offset + axis_value
        print ('realoffset easy')
    except:
        print ('realoffset except', len(offset))
        realoffset = []
        for i in range(len(offset)):
            realoffset.append(offset[i] + axis_value)

    contour1.Isosurfaces = realoffset
    contour1.PointMergeMethod = 'Uniform Binning'

    if clipxmin < clipxmax:
        if axis_index == 0:
            raise ValueError('Tried to set boundaries in the main axis')
        clip_xmin = pvsimple.Clip(contour1)
        clip_xmin.ClipType = "Plane"
        clip_xmin.ClipType.Normal = [1.0, 0.0, 0.0]
        clip_xmin.ClipType.Origin = [clipxmin, 0.0, 0.0]
        if usehera:
            clip_xmin.Invert = 0
        else:
            clip_xmin.InsideOut = 0
        clip_xmin.UpdatePipeline()
        clip_xmax = pvsimple.Clip(clip_xmin)
        clip_xmax.ClipType = "Plane"
        clip_xmax.ClipType.Normal = [1.0, 0.0, 0.0]
        clip_xmax.ClipType.Origin = [clipxmax, 0.0, 0.0]
        if usehera:
            clip_xmin.Invert = 1
        else:
            clip_xmin.InsideOut = 1
        clip_xmax.UpdatePipeline()
        contour1 = clip_xmax
    if clipymin < clipymax:
        if axis_index == 1:
            raise ValueError('Tried to set boundaries in the main axis')
        clip_ymin = pvsimple.Clip(contour1)
        clip_ymin.ClipType = "Plane"
        clip_ymin.ClipType.Normal = [0.0, 1.0, 0.0]
        clip_ymin.ClipType.Origin = [0.0, clipymin, 0.0]
        if usehera:
            clip_ymin.Invert = 0
        else:
            clip_ymin.InsideOut = 0
        clip_ymin.UpdatePipeline()
        clip_ymax = pvsimple.Clip(clip_ymin)
        clip_ymax.ClipType = "Plane"
        clip_ymax.ClipType.Normal = [0.0, 1.0, 0.0]
        clip_ymax.ClipType.Origin = [0.0, clipymax, 0.0]
        if usehera:
            clip_ymin.Invert = 1
        else:
            clip_ymin.InsideOut = 1
        clip_ymax.UpdatePipeline()
        contour1 = clip_ymax
    if clipzmin < clipzmax:
        if axis_index == 2:
            raise ValueError('Tried to set boundaries in the main axis')
        print('testclipz')
        clip_zmin = pvsimple.Clip(contour1)
        clip_zmin.ClipType = "Plane"
        clip_zmin.ClipType.Normal = [0.0, 0.0, 1.0]
        clip_zmin.ClipType.Origin = [0.0, 0.0, clipzmin]
        if usehera:
            clip_zmin.Invert = 0
        else:
            clip_zmin.InsideOut = 0
        clip_zmin.UpdatePipeline()
        clip_zmax = pvsimple.Clip(clip_zmin)
        clip_zmax.ClipType = "Plane"
        clip_zmax.ClipType.Normal = [0.0, 0.0, 1.0]
        clip_zmax.ClipType.Origin = [0.0, 0.0, clipzmax]
        if usehera:
            clip_zmin.Invert = 1
        else:
            clip_zmin.InsideOut = 1
        clip_zmax.UpdatePipeline()
        contour1 = clip_zmax

    itr = bse.to_pandas('case name', contour1, timelist=itteration)
    try:
        db = []
        for i in range(len(itteration)):
            db.append(itr.next())
    except:
        db = itr.next()
    return db


def get_clip(filename, itteration,
             clipxmin=0.0, clipxmax=0.0, clipymin=0.0, clipymax=0.0, clipzmin=0.0, clipzmax=0.0, reader=None, parallel=False):
    if usehera:
        bse = paraviewOpenFOAM(shrinkdir(filename))
        # bse = paraviewOpenFOAM(dirname)
        # reader=bse.reader
        clip_1 = pvsimple.Clip(Input = bse.reader,guiName="myClip1")
        # clip_1 = reader
    else:
        if reader is None:
            if parallel:
                reader = bse.ReadCase(filename, filename, CaseType='Decomposed Case')  # ' Reconstructed Decomposed Case')
            else:
                reader = bse.ReadCase(filename, filename, CaseType='Reconstructed Case')  # ' Reconstructed Decomposed Case')
            print('Initial timelist print:', len(reader.TimestepValues), reader.TimestepValues[-1])
            print ('first time takes time')
            clip_1 = reader
        else:
            print ('save time and memory')
            clip_1 = reader
    
    if clipxmin < clipxmax:
        # if axis_index == 0:
        #     raise ValueError('Tried to set boundaries in the main axis')
        clip_xmin = pvsimple.Clip(clip_1)
        clip_xmin.ClipType = "Plane"
        clip_xmin.ClipType.Normal = [1.0, 0.0, 0.0]
        clip_xmin.ClipType.Origin = [clipxmin, 0.0, 0.0]
        if usehera:
            clip_xmin.Invert = 0
        else:
            clip_xmin.InsideOut = 0
        clip_xmin.UpdatePipeline()
        clip_xmax = pvsimple.Clip(clip_xmin)
        clip_xmax.ClipType = "Plane"
        clip_xmax.ClipType.Normal = [1.0, 0.0, 0.0]
        clip_xmax.ClipType.Origin = [clipxmax, 0.0, 0.0]
        if usehera:
            clip_xmin.Invert = 1
        else:
            clip_xmin.InsideOut = 1
        clip_xmax.UpdatePipeline()
        clip_1 = clip_xmax
    if clipymin < clipymax:
        # if axis_index == 1:
            # raise ValueError('Tried to set boundaries in the main axis')
        clip_ymin = pvsimple.Clip(clip_1)
        clip_ymin.ClipType = "Plane"
        clip_ymin.ClipType.Normal = [0.0, 1.0, 0.0]
        clip_ymin.ClipType.Origin = [0.0, clipymin, 0.0]
        if usehera:
            clip_ymin.Invert = 0
        else:
            clip_ymin.InsideOut = 0
        clip_ymin.UpdatePipeline()
        clip_ymax = pvsimple.Clip(clip_ymin)
        clip_ymax.ClipType = "Plane"
        clip_ymax.ClipType.Normal = [0.0, 1.0, 0.0]
        clip_ymax.ClipType.Origin = [0.0, clipymax, 0.0]
        if usehera:
            clip_ymin.Invert = 1
        else:
            clip_ymin.InsideOut = 1
        clip_ymax.UpdatePipeline()
        clip_1 = clip_ymax
    if clipzmin < clipzmax:
        # if axis_index == 2:
            # raise ValueError('Tried to set boundaries in the main axis')
        clip_zmin = pvsimple.Clip(clip_1)
        clip_zmin.ClipType = "Plane"
        clip_zmin.ClipType.Normal = [0.0, 0.0, 1.0]
        clip_zmin.ClipType.Origin = [0.0, 0.0, clipzmin]
        if usehera:
            clip_zmin.Invert = 0
        else:
            clip_zmin.InsideOut = 0
        clip_zmin.UpdatePipeline()
        clip_zmax = pvsimple.Clip(clip_zmin)
        clip_zmax.ClipType = "Plane"
        clip_zmax.ClipType.Normal = [0.0, 0.0, 1.0]
        clip_zmax.ClipType.Origin = [0.0, 0.0, clipzmax]
        if usehera:
            clip_zmin.Invert = 1
        else:
            clip_zmin.InsideOut = 1
        clip_zmax.UpdatePipeline()
        clip_1 = clip_zmax

    if 5==6:
        clip_1 = pvsimple.CellCenters(clip_1)
    clip_1.UpdatePipeline()

    
    # db = []
    # for i in itr:
    #     dbb = itr.next()
    #     db.append(dbb[dbb.keys()[0]])
    #     if len(dbb[dbb.keys()[0]]) == 0:
    #         print('There is no records, please check coordinates limits or slice location')


    if usehera:
        x=bse.to_pandas("myClip1",timelist=itteration)
        db=[]
        for z in x:
            # y=next(z)
            db.append(z['myClip1'])
        if len(db)==1:
            db=db[0]
    else:
        # itr = bse.to_pandas(reader, clip_1, timelist=itteration)
        itr = bse.to_pandas(reader, clip_1, timelist=999999)
        if isinstance(itteration, list):
            # print('type list')
            db = []
            for i in range(len(itteration)):
                dbb = itr.next()
                db.append(dbb[dbb.keys()[0]])
                if len(dbb[dbb.keys()[0]]) == 0:
                    print('There are no records in get clip 1, please check coordinates limits or slice location')
        else:
            dbb = itr.next()
            db = dbb[dbb.keys()[0]]
            if len(db) == 0:
                print('There are no records in get clip 2, please check coordinates limits or slice location')

    return db


def density_estimation(m1, m2):
    X, Y = np.mgrid[min(m1):max(m1):100j, min(m2):max(m2):100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def flood(floor, i=0, j=0, max_floor=0):
    # print ('ijk',i,j,floor[i,j])
    if floor[i, j] == -999:
        floor[i, j] = max_floor
        if i>0:
            floor = flood(floor, i-1, j, max_floor)
            if j>0:
                floor = flood(floor, i - 1, j-1, max_floor)
            if j<floor.shape[1]-1:
                floor = flood(floor, i - 1, j+1, max_floor)
        if j>0:
            floor = flood(floor, i, j - 1, max_floor)
        if j<floor.shape[1]-1:
            floor = flood(floor, i, j+1, max_floor)
        if i<floor.shape[0]-1:
            floor = flood(floor, i + 1, j, max_floor)
            if j>0:
                floor = flood(floor, i + 1, j-1, max_floor)
            if j < floor.shape[1] - 1:
                floor = flood(floor, i + 1, j+1, max_floor)
    return floor


def flood_main(floor):
    sys.setrecursionlimit(floor.shape[0]*floor.shape[1])  # 10000 is an example, try with different values
    print('flood size:',floor.shape[0], floor.shape[1])
    max_floor = 1
    floor[floor == 0] = -999
    floor[floor != -999] = 0
    for i in range(floor.shape[0]):
        for j in range(floor.shape[1]):
            if floor[i, j] == -999:
                floor = flood(floor, i, j, max_floor)
                print('flood', i, j, max_floor)
                max_floor += 1
    print('flood fin')
    return floor


# if floor[i - 1, j - 1] == -999: floor = flood(floor, i - 1, j - 1, max_floor)
# if floor[i, j - 1] == -999: floor = flood(floor, i, j - 1, max_floor)
# if floor[i + 1, j - 1] == -999: floor = flood(floor, i + 1, j - 1, max_floor)
# if floor[i - 1, j] == -999: floor = flood(floor, i - 1, j, max_floor)
# if floor[i + 1, j] == -999: floor = flood(floor, i + 1, j, max_floor)
# if floor[i - 1, j + 1] == -999: floor = flood(floor, i - 1, j + 1, max_floor)
# if floor[i, j + 1] == -999: floor = flood(floor, i, j + 1, max_floor)
# if floor[i + 1, j + 1] == -999: floor = flood(floor, i + 1, j + 1, max_floor)


def remove_left(s):
    return s[1:]


def remove_right(s):
    return s[:-1]


def read_measurement(textboxfile1value, textboxtime1value, mainsimulation, showfieldindex):

    if float(textboxtime1value).is_integer():
        print('probe int')
    else:
        textboxtime1value = 1
        print('probe0:',textboxtime1value)
    probefile = textboxfile1value[:textboxfile1value.rfind('/')] + '/postProcessing/probes/0/U'
    with open(probefile) as f:
        a = [line.split() for line in f]
    items = a.index(['#', 'Time']) - 1  # if items is 115 that means that we have 0-114 items, first time is 117
    b = a[int(textboxtime1value + items + 1)]
    b[1::3] = [remove_left(s) for s in b[1::3]]
    b[3::3] = [remove_right(s) for s in b[3::3]]
    c = np.array(b[1:]).reshape(int((len(b) - 1) / 3), 3)

    d = np.zeros(len(c))
    for i in range(len(d)):
        # print ('i', i)
        # print ('c', c[i, 0], c[i, 1], c[i, 2])
        if abs(float(c[i, 0])) < 1e10:
            d[i] = (float(c[i, 0]) ** 2 + float(c[i, 1]) ** 2 + float(c[i, 2]) ** 2) ** 0.5
        else:
            d[i] = 0

    model = d.tolist()

    if mainsimulation == 'observationE':
        modelE = np.zeros(len(model) - 1)
        for i in range(len(modelE)):
            modelE[i] = d[i + 1] #/ d[0]
        model = modelE
        print('mainsimulation-observationE:', mainsimulation)
    elif mainsimulation == 'observationB':
        model = [float(i) for i in c[:, showfieldindex]]  # c[:,0]
        print('mainsimulation-observationB:', mainsimulation)
    elif mainsimulation == 'Michelstadt':
        model = [float(i) for i in c[:, showfieldindex]]  # c[:,0]
        print('mainsimulation-Michel:', mainsimulation)
    model = np.asarray(model)

    if mainsimulation == 'observationE':
        coordinates = np.zeros([len(d)-1, 3])
        for i in range(len(d)-1):
            coordinates[i][0] = float(a[i+1][3][1:])
            coordinates[i][1] = float(a[i+1][4][:])
            coordinates[i][2] = float(a[i+1][5][:-1])
    else:
        coordinates = np.zeros([len(d), 3])
        for i in range(len(d)):
            coordinates[i][0] = float(a[i][3][1:])
            coordinates[i][1] = float(a[i][4][:])
            coordinates[i][2] = float(a[i][5][:-1])

    return model, coordinates


def log_profile(x, a, b, c, d):
#        u[i] = uh * beta / k * math.log(((z[i] - h) + d) / z0)
    return a * np.log(b * (x + d)) + c


def exp_profile(x, a, b, c):
#        u[i] = uh * math.exp(beta * (z[i] - h) / lvar)
#        return Uh * np.exp()
    return a * np.exp(-b * x) + c


def logexp_profile(x, a, b, c, d, e, g):
#        u[i] = uh * math.exp(beta * (z[i] - h) / lvar)
#        return Uh * np.exp()
    return ((np.arctan(g * x )/(math.pi/2) + 1) / 2) * (a * np.log(b * (x + c))) + ((np.arctan(-g*x)/(math.pi/2) + 1) / 2) * (d * np.exp(-e * x))

def test_profile(x, a, b, c):
#        u[i] = uh * beta / k * math.log(((z[i] - h) + d) / z0)
    return a * np.log(b * x + c)

def find_index(xcoord, ycoord, zcoord, x, y, z):
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    dz = z[1]-z[0]
    xindex = (np.abs(x - xcoord)).argmin()
#    if xcoord<x[xindex]:
    x0=xindex
#    else:
#        x0=xindex-1
    xindex=x0+(xcoord-x[x0])/dx
    yindex = (np.abs(y - ycoord)).argmin()
#    if ycoord<y[yindex]:
    y0=yindex
#    else:
#        y0=yindex-1
    yindex=y0+(ycoord-y[y0])/dy
    zindex = (np.abs(z - zcoord)).argmin()
#    if zcoord<z[zindex]:
    z0=zindex
#    else:
#        z0=zindex-1
    zindex=z0+(zcoord-z[z0])/dz
    return xindex,yindex, zindex

def traj(xcoord, ycoord, zcoord, x, y, z, gridux,griduy, griduz, duration = 200, dt = 0.1, verbose = False):
    if 5==7:
        a=random.randint(-5,5)
        xcoord0=xcoord
        ycoord0=ycoord
        zcoord0=zcoord

        xcoord=xcoord0
        ycoord=ycoord0
        zcoord=zcoord0
        
        verbose = True
        duration = 2000000
        dt = 0.1

    xcoordbkup = xcoord
    ycoordbkup = ycoord
    zcoordbkup = zcoord
    
    trajx = []
    trajy = []
    trajz = []
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    dz = z[1]-z[0]
    d=min(dx,dy,dz)
#    print (dx,dy,dz)

    for t in range(duration):
        xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
        trajx.append(xindex)
        trajy.append(yindex)
        trajz.append(zindex)
        if zindex<-1:
            print('zindex<-1',zcoordbkup, zcoord, dt, dz)
            break
        xindex = max(0, xindex)
        yindex = max(0, yindex)
        zindex = max(0, zindex)
#        xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
#        oldu,oldv,oldw = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
#        oldspeed=max(math.fabs(oldu),math.fabs(oldv),math.fabs(oldw))
        u,v,w = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
        speed=max(math.fabs(u),math.fabs(v),math.fabs(w))
        if speed==0.0:
            xcoord1=xcoord
            ycoord1=ycoord
            zcoord1=zcoord

            tries=0
            while (speed==0) and (tries<10):
                tries+=1
                u = random.randint(-5,5)
                xcoord = xcoordbkup + dt * u / dx
                ycoord = ycoordbkup + dt * v / dy
                zcoord = zcoordbkup + dt * w / dz               
                xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
                u,v,w = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
                speed=max(math.fabs(u),math.fabs(v),math.fabs(w))
                print(tries, u, v, w, 'bounce0',t,xcoord,xcoordbkup,ycoord,ycoordbkup,zcoord,zcoordbkup)
                

        if speed==0.0:
            zcoord = zcoordbkup
            ycoord = ycoord1
            xcoord = xcoord1
            print('bounce1',t,xcoord,xcoordbkup,ycoord,ycoordbkup,zcoord,zcoordbkup)
            xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
            u,v,w = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
            speed=max(math.fabs(u),math.fabs(v),math.fabs(w))
        if speed==0.0:
            print('bounce2',t,xcoord,xcoordbkup,ycoord,ycoordbkup,zcoord,zcoordbkup)
            zcoord = zcoord1
            ycoord = ycoordbkup
            xcoord = xcoord1
            xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
            u,v,w = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
            speed=max(math.fabs(u),math.fabs(v),math.fabs(w))
        if speed==0.0:
            zcoord = zcoord1
            ycoord = ycoord1
            xcoord = xcoordbkup
            print('bounce3',t,xcoord,xcoordbkup,ycoord,ycoordbkup,zcoord,zcoordbkup)
            xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
            u,v,w = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
            speed=max(math.fabs(u),math.fabs(v),math.fabs(w))
        if speed==0.0:
            print('bounce4',t,xcoord,xcoordbkup,ycoord,ycoordbkup,zcoord,zcoordbkup)
            zcoord = zcoord1
            ycoord = ycoordbkup
            xcoord = xcoordbkup
            xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
            u,v,w = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
            speed=max(math.fabs(u),math.fabs(v),math.fabs(w))
        if speed==0.0:
            print('bounce5',t,xcoord,xcoordbkup,ycoord,ycoordbkup,zcoord,zcoordbkup)
            xcoord = xcoord1
            ycoord = ycoordbkup
            zcoord = zcoordbkup
            xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
            u,v,w = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
            speed=max(math.fabs(u),math.fabs(v),math.fabs(w))
        if speed==0.0:
            print('bounce6',t,xcoord,xcoordbkup,ycoord,ycoordbkup,zcoord,zcoordbkup)
            xcoord = xcoordbkup
            ycoord = ycoord1
            zcoord = zcoordbkup
            xindex, yindex, zindex = find_index(xcoord, ycoord, zcoord, x, y, z)
            u,v,w = interpvel3d(gridux, griduy, griduz, xindex, yindex, zindex)
            speed=max(math.fabs(u),math.fabs(v),math.fabs(w))
            
        dt=d/speed
#        if t == 0:
#            print('speed:', u, v, w)
        xcoordbkup = xcoord
        ycoordbkup = ycoord
        zcoordbkup = zcoord
        xcoord += dt * u / dx
        ycoord += dt * v / dy
        zcoord += dt * w / dz
        if zcoord < z.min():
#            zcoord -= dt * w / dz
            zcoord = z.min()
            
#        xindexnew = (np.abs(x - xcoord)).argmin()
#        if xcoord<x[xindexnew]:
#            x0=xindexnew
#        else:
#            x0=xindexnew-1
#        xindexnew=x0+(xindexnew-x[x0])/dx
#        
#        yindexnew = (np.abs(y - ycoord)).argmin()
#        if ycoord<y[yindexnew]:
#            y0=yindexnew
#        else:
#            y0=yindexnew-1
#        yindexnew=y0+(ycoord-y[y0])/dy
#        
#        zindexnew = (np.abs(z - zcoord)).argmin()
#        if zcoord<z[zindexnew]:
#            z0=zindexnew
#        else:
#            z0=zindexnew-1
#        zindexnew=z0+(zcoord-z[z0])/dz
#        if gridux[xindexnew, yindexnew, zindexnew] == 0:
#            if gridux[xindexnew, yindex, zindexnew] == 0:
#                xcoord -= dt * u / dx
#            if gridux[xindex, yindexnew, zindexnew] == 0:
#                ycoord -= dt * v / dy
        if verbose:
            print (t,dt, xcoord, ycoord, zcoord, xindex, yindex, zindex, u, v, w)
        if u==0.0:
            print('fin traj u==0')
            break
#        if zcoord<0:
#            break
        if len(trajx)>2:
            if ((xindex==trajx[-2]) and (yindex==trajy[-2]) and (zindex==trajz[-2])):
                print('fin traj stay in place')
                break

#    print('end traj', datetime.datetime.now())
    if 5==7:
        plt.figure()
        plt.imshow(gridux[:,:,int(trajz[0])].T, origin='lower')  # jet, Paired
        plt.colorbar()
        plt.plot(trajx, trajy, marker="+")
        plt.xlabel('x')
        plt.ylabel('y')
#        plt.plot(trajx2, trajy2, marker="*")
        plt.show()
        plt.figure()
        plt.plot(trajx, trajz, marker="+")
        plt.xlabel('x')
        plt.ylabel('z')        
        
        
    return (trajx, trajy, trajz, t)
    # plt.figure()
    # plt.imshow(gridux[:, :, zindex], origin='lower', cmap=plt.get_cmap('tab20c'))  # jet, Paired
    # plt.plot(trajx, trajy, marker="o")
    # plt.colorbar()
    # plt.show()
    
def morphology(textboxfile1value,textboxtime1value , x,y,radius=100):
    
# This funcrion will return the morphology parameter of an area base on the mesh

# Parameters:
# x coordinate
# y coordinate    
# radius, number of meter from each side, rectangle shape

# return
# he - median height
# lambdaP
# lambdaF
# LC
       
#    textboxfile1value=[r'/ibdata2/nirb/openFOAM/angle/klembenchmark/cylinder.foam']
#    x=219000
#    y=632000
#    textboxfile1value =r'/ibdata2/nirb/openFOAM/porous/hadasfine/windAroundCube.foam' ## self.filename1.text()
#    textboxtime1value = 200 # float(self.itteration1.text())
#    x=179000
#    y=664000
#    radius=1000

    db = get_clip(textboxfile1value, [textboxtime1value],
               clipxmin=x-radius, clipxmax=x+radius,
               clipymin=y-radius, clipymax=y+radius,
               clipzmin=0.0, clipzmax=0.0)


    # ap3 = db[['x', 'y', 'z', sel[0].data()]]  # at the ground (1 meter)
    ap3 = db[['x', 'y', 'z', 'U_x']]  # at the ground (1 meter)
    points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
    # values = ap3[sel[0].data()]
    valuesux = ap3['U_x']

#    if xmin==xmax:
    xmin=ap3['x'].min()
    xmax=ap3['x'].max()
#    if ymin==ymax:
    ymin=ap3['y'].min()
    ymax=ap3['y'].max()
#    if zmin==zmax:
    zmin=ap3['z'].min()
    zmax=ap3['z'].max()
    # ticks = 600j
    grid_points = 42000000  # add 00
    meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)

    resolution = (meters / grid_points)**(1./3)
    if resolution==0:
        print('no clip boundaries in one of the axis')
    # resolution = 0.1
    ticksx = int((xmax-xmin) / resolution) * 1j
    ticksy = int((ymax-ymin) / resolution) * 1j
    ticksz = int((zmax-zmin) / resolution) * 1j
    print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)

    xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
    gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                       method='nearest')  # Nearest for keeping zeros that are buildings
    
 
    xi1=xi[:,0,0]
    yi1=yi[0,:,0]
    zi1=zi[0,0,:]
#    ximin=0
#    ximax=len(xi1)
#    for i in range(len(xi1)):
#        if xi1[i]<x-radius:
#            ximin = i
#    for i in range(len(xi1)):
#        if xi1[len(xi1)-1-i]>x+radius:
#            ximax = i
    dx= xi1[4]-xi1[3]
    dy= yi1[4]-yi1[3]
    he=0
    for i in range(len(zi1)-1):
        he += zi1[i]*(np.count_nonzero(gridux[:,:,i]==0)-np.count_nonzero(gridux[:,:,i+1]==0))
    he /= np.count_nonzero(gridux[:,:,0]==0)
    
    lambdap=(1-float(np.count_nonzero(gridux[:,:,0]==0))/(gridux[:,:,0].shape[0]*gridux[:,:,0].shape[1]))*dx*dy

    lambdaf=lambdap
    
    lc = he*(1-lambdap)/lambdaf

    return he, lambdaf, lambdap, lc

def windprofile(z, uref=3, href=24, he=24, lambdap=0.3, lambdaf=0.3,beta=0.3, verbose=False):
    
	# This function will return the wind velocity at height z 

	# parameters:
	# z is the predicted height [m] above the ground
	# uref is the velocity reference [m/s] at height href
	# href is the reference height [m] of uref
	# he is the mean building height [m] above the ground
	# lambdap is the building percent when looking from top
	# lambdaf is the building percent when looking from the side
	# beta is the ratio between ustar and Uh

	# return
	# u is the predicted velocity [m/s] that is calculated in this function

	# variables:
	# l is the mixing length [m]

	###############################
	##    z=24
	#    uref=3
	#    href=24
	#    he=24
	#    lambdap=0.3
	#    lambdaf=0.3
	###############################

         if lambdaf==0:
#             print('lambdaf==0')
             z0 = 0.03 # open flat terrain
             if z>0:
                 u = uref * math.log((z+z0)/z0) / math.log((href+z0)/z0)  # I add z0 instead of d to avoid negetive velocities and z0 is neglectable to the vertical resolution, so it's OK
             else:
                 u = np.nan
         else:
             k = 0.41 # von Karman constant
    	    # beta = beta# For closed uniform natural canopies, Raupach 1996
    	    # beta = ustar / uh
    	    # LA commercial area has lambdap=0.28, lambdaf=0.27 he=24.5m, LC=66, Coceal 2004
             lc = he*(1-lambdap)/lambdaf
    	    # calculating the mixing length at the canopy only to know d
             l = 2. * beta**3 * lc # mixing length in the canopy
             d = l/k  # displacement height
             z0 = l/k*math.exp(-k/beta) # surface roughness
    
             if href>he: 
                uh  = uref * k /beta / math.log(((href-he)+d)/z0)
             else: # href<he
                uh  = uref / math.exp(beta*(href-he)/l)
#    	       print ('there is no sense of giving uref below building mean height because we have to give mean velocity, above he, one point represent the mean value, below, it is not represent')
#             print('uh=',uh)
             if z <= he: # below the canopy
                u = uh * math.exp(beta*((z-he))/l)
             else: # above the canopy
                u = uh*beta/k*math.log(((z-he)+d)/z0)
             if z<0:
                 u=np.nan
         if verbose is True: 
             print('windprofile', round(z,4), round(he,4), round(href,4), round(z0,4), round(d,4), round(u,4), round(uh,4))
         return u

         
def histogram(array, title):
    sortarray = np.sort(array, axis=None)
    l = len(array)
    plt.figure()
    plt.hist(sortarray[l/10:-l/10+1],10)
    plt.title(title)
    plt.show()
         
def corrfunc(array, x, y, title):
    cf=np.zeros(int(math.ceil(((x.max()-x.min())**2.+(y.max()-y.min())**2.)**0.5)+1.))
    sf=np.zeros(math.ceil(((x.max()-x.min())**2.+(y.max()-y.min())**2.)**0.5)+1.)
    sfm=np.zeros(math.ceil(((x.max()-x.min())**2.+(y.max()-y.min())**2.)**0.5)+1.)
    sfu=np.zeros(math.ceil(((x.max()-x.min())**2.+(y.max()-y.min())**2.)**0.5)+1.)
    if len(array)<20000:
        for i in range(len(array)):
            if i%200==0:
                print('corrfunc ', i, len(array))
            for j in range(len(array)):
    #            print(x[i],y[j],math.ceil((x[i]**2.+y[j]**2.)**0.5))
                r=int(math.ceil(((x[i]-x.min())**2.+(y[j]-y.min())**2.)**0.5))
                cf[r]+=array[i]*array[j]
                sf[r]+=math.fabs(array[i]-array[j])
                sfm[r]+=array[i]-array[j]
                if y[i]>=y[j]:
                    sfu[r]+=array[i]-array[j]
        plt.figure()
        plt.plot(cf/cf.max())
        plt.title('cf * '+title)
        plt.show()
#        plt.figure()
#        plt.plot(sf/sf.max())
#        plt.title('sf - '+title)
#        plt.show()        
#        plt.figure()
#        plt.plot(sfm/sfm.max())
#        plt.title('sfm - '+title)
#        plt.show()
#        plt.figure()
#        plt.plot(sfu/sfu.max())
#        plt.title('sfu - '+title)
#        plt.show()
         
def chi(uw0, uw1, main=True):
    precentile=.9
    uw0a=uw0.copy().ravel()
    uw1a=uw1.copy().ravel()
    count = len(uw0a)
    if main:
        uw0a=np.sort(uw0a)
        uw1a=np.sort(uw1a)
#        uw0a=uw0a[int(((1-precentile)/2.)*count):int(((1.+precentile)/2.)*count)]
#        uw1a=uw1a[int(((1-precentile)/2.)*count):int(((1.+precentile)/2.)*count)]
        uw0a=uw0a[:int(((1.+precentile)/2.)*count)]
        uw1a=uw1a[:int(((1.+precentile)/2.)*count)]
    lower=min(uw0a.min(),uw1a.min())
    upper=max(uw0a.max(),uw1a.max())
    binss=np.linspace(lower,upper,1000)
    hist0,bins=np.histogram(uw0a,bins=binss)
    hist1,bins=np.histogram(uw1a,bins=binss)
#    print(bins,hist0,hist1)
#    shist0=hist0.sum()
#    shist1=hist1.sum()
#    print(shist0)
#    hist0/=float(shist0)
#    hist1/=float(shist1)
    chi2=0.
    for i in range(len(hist0)):
        h0=hist0[i]
        h1=hist1[i]
        if ((h0+h1)!=0.0):
            chi2+=(h0-h1)**2./float(h0+h1)
#            print('chichi',i, chi2,h0,h1)
    chi2/=2.0
    chi2/=count
    return round(chi2,3)
         
def interpvel(matu,matw,i,j):
#    u = matu[int(i),int(j)]
#    w = matw[int(i),int(j)]
    fraci=i-int(i)
    fracj=j-int(j)

    if i+1<matu.shape[0] and j<matu.shape[1]: 
        if np.isnan(matu[int(i+1),int(j)]):
            u = matu[int(i),int(j)]
        else:
            u = matu[int(i+1),int(j)]*fraci+matu[int(i),int(j)]*(1-fraci)
    else:
        u = np.NaN #matu[min(int(i),matu.shape[0]-1),min(int(j),matu.shape[1]-1)]

    if j+1<matu.shape[1] and i<matu.shape[0]: 
        if np.isnan(matw[int(i),int(j+1)]):
            w = matw[int(i),int(j)]
        else:      
            w = matw[int(i),int(j+1)]*fracj+matw[int(i),int(j)]*(1-fracj)
    else:
        w = np.NaN #matw[min(int(i),matu.shape[0]-1),min(int(j),matu.shape[1]-1)]
        
    return [u, w]


def interpvel3d(matu, matv, matw, i, j, k):
#    u = matu[int(i),int(j)]
#    w = matw[int(i),int(j)]
    i0=int(i)
    j0=int(j)
    k0=int(k)
    i1=int(i)+1
    j1=int(j)+1
    k1=int(k)+1
    fraci=i-i0
    fracj=j-j0
    frack=k-k0
    # we assume that corner is a wall, need to fix it someday
    if i1<matu.shape[0]-1 and k1<matu.shape[2]-1 and j1<matu.shape[1]-1 and i0>=0 and j0>=0 and k0>=0: 
        u000=matu[i0,j0,k0]
        u001=matu[i1,j0,k0]
        u010=matu[i0,j1,k0]
        u011=matu[i1,j1,k0]
        u100=matu[i0,j0,k1]
        u101=matu[i1,j0,k1]
        u110=matu[i0,j1,k1]
        u111=matu[i1,j1,k1]
        v000=matv[i0,j0,k0]
        v001=matv[i1,j0,k0]
        v010=matv[i0,j1,k0]
        v011=matv[i1,j1,k0]
        v100=matv[i0,j0,k1]
        v101=matv[i1,j0,k1]
        v110=matv[i0,j1,k1]
        v111=matv[i1,j1,k1]
        w000=matw[i0,j0,k0]
        w001=matw[i1,j0,k0]
        w010=matw[i0,j1,k0]
        w011=matw[i1,j1,k0]
        w100=matw[i0,j0,k1]
        w101=matw[i1,j0,k1]
        w110=matw[i0,j1,k1]
        w111=matw[i1,j1,k1]
        if np.isnan(u000) or u000==0:
            u00=u001
            v00=v001
            w00=w001
        elif np.isnan(u001) or u001==0:
            u00=u000
            v00=v000
            w00=w000
        else:
            u00 = u000*fraci+u001*(1.-fraci)
            v00 = v000*fraci+v001*(1.-fraci)
            w00 = w000*fraci+w001*(1.-fraci)
            
        if np.isnan(u010) or u010==0:
            u01=u011
            v01=v011
            w01=w011
        elif np.isnan(u011) or u011==0:
            u01=u010
            v01=v010
            w01=w010
        else:          
            u01 = u010*fraci+u011*(1.-fraci)
            v01 = v010*fraci+v011*(1.-fraci)
            w01 = w010*fraci+w011*(1.-fraci)
            
        if np.isnan(u100) or u100==0:
            u10=u101
            v10=v101
            w10=w101
        elif np.isnan(u101) or u101==0:
            u10=u100
            v10=v100
            w10=w100
        else:
            u10 = u100*fraci+u101*(1.-fraci)
            v10 = v100*fraci+v101*(1.-fraci)
            w10 = w100*fraci+w101*(1.-fraci)
        
        if np.isnan(u110) or u110==0:
            u11=u111
            v11=v111
            w11=w111
        elif np.isnan(u111) or u111==0:
            u11=u110
            v11=u110
            w11=u110
        else:                        
            u11 = u110*fraci+u111*(1.-fraci)
            v11 = v110*fraci+v111*(1.-fraci)
            w11 = w110*fraci+w111*(1.-fraci)
            
        if np.isnan(u00) or u00==0:
            u0=u01
            v0=v01
            w0=w01
        elif np.isnan(u01) or u01==0:
            u0=u00
            v0=v00
            w0=w00
        else:                                   
            u0 = u00*fracj+u01*(1.-fracj)
            v0 = v00*fracj+v01*(1.-fracj)
            w0 = w00*fracj+w01*(1.-fracj)
        
        if np.isnan(u10) or u10==0:
            u1=u11
            v1=v11
            w1=w11
        elif np.isnan(u11) or u11==0:
            u1=u10
            v1=v10
            w1=w10
        else:                                   
            u1 = u10*fracj+u11*(1.-fracj) 
            v1 = v10*fracj+v11*(1.-fracj) 
            w1 = w10*fracj+w11*(1.-fracj) 
            
        if np.isnan(u0) or u0==0:
            u=u1
            v=v1
            w=w1
        elif np.isnan(u1) or u1==0:
            u=u0
            v=u0
            w=u0
        else:                                   
            u  = u0*frack+u1*(1.-frack)
            v  = v0*frack+v1*(1.-frack)
            w  = w0*frack+w1*(1.-frack)
    else:
        u = np.NaN 
        v = np.NaN 
        w = np.NaN 

        
    return [u, v, w]

        
def rpaintinside(mat,cur,point, border):
#    print('recursive flood fill')

    if (not((mat[cur[0],cur[1]]==border) | (mat[cur[0],cur[1]]==point))):
        mat[cur[0],cur[1]]=point    
        if (cur[0]<mat.shape[0]-2):
            print(mat[cur[0]+1,cur[1]],point,border)
            paintinside(mat,[cur[0]+1,cur[1]],point,border)
        if cur[0]>0:
            paintinside(mat,[cur[0]-1,cur[1]],point,border)
        if (cur[1]<mat.shape[1]-2):
            paintinside(mat,[cur[0],cur[1]+1],point,border)
        if cur[1]>0:
            paintinside(mat,[cur[0],cur[1]-1],point,border)

def paintinside(mat,cur,point, border):
#    print('non rec flood fill')
    matbkup=mat.copy()
    matbkup[cur[0],cur[1]]=point
    found = 9999
    while (found>0):
        found=0
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if (matbkup[i,j]==point):
                    if i>0:
                       if (matbkup[i-1,j]==0):
                           found+=1
                           matbkup[i-1,j]=point
                    if j>0:
                       if (matbkup[i,j-1]==0):
                           found+=1
                           matbkup[i,j-1]=point                        
                    if (i<matbkup.shape[0]-2):
                       if (matbkup[i+1,j]==0):
                           found+=1
                           matbkup[i+1,j]=point                        
                    if (j<matbkup.shape[1]-2):
                       if (matbkup[i,j+1]==0):
                           found+=1
                           matbkup[i,j+1]=point
#        print('found ',found)
    mat[matbkup==point]=point    


def eddie_mask_gradient(psi1):
    
    if 5==6:
        itter=1000000000
        files=[r'/ibdata2/nirb/openFOAM/ml/windAroundurbanMichelstadt2zomegablock1/windAroundCube.foam']
        axisindex = 1
        axispos = -250000#100000
        xminc=0#300000
        xmaxc=0#300000
        yminc=0
        ymaxc=0
        zminc=0
        zmaxc=0
        currentRow=1

        textboxfile1value = files[0]
        db = get_slice(textboxfile1value, itter, axisindex, axispos,
                       clipxmin=xminc, clipxmax=xmaxc,
                       clipymin=yminc, clipymax=ymaxc,
                       clipzmin=zminc, clipzmax=zmaxc)

        ticks = 400
        xmin, xmax, ymin, ymax, u = makegrid(db, 'U_x', currentRow, ticks=ticks)
        xmin, xmax, ymin, ymax, w = makegrid(db, 'U_z', currentRow, ticks=ticks)
        u[u==0]= np.NaN
        w[w==0]= np.NaN
    
        xi = np.linspace(0,ticks,ticks)
        zi = np.linspace(0,ticks,ticks)
        X,Z=np.meshgrid(xi,zi)
        w[np.isnan(w)]=0
        u[np.isnan(u)]=0
    
        intx=integrate.cumtrapz(w,X,axis=1,initial=0)[0]
        inty=integrate.cumtrapz(u,Z,axis=0,initial=0)
        psi1=intx-inty     



#        w[np.isnan(w)]=0
#        u[np.isnan(u)]=0
#        u[np.abs(u)<0.000001]=0
#        w[np.abs(w)<0.000001]=0
#        intx=integrate.cumtrapz(w,xi,axis=1,initial=0)[0]
#        inty=integrate.cumtrapz(u,zi,axis=0,initial=0)
#        psi1=intx-inty
#        w[w==0]=np.nan
#        u[u==0]=np.nan
              
    zi1=psi1.copy()
#    zi1[np.abs(u)<0.0001] = np.nan

# multiple neibour find negative (sign change)
    psi1u=np.roll(psi1,1,axis=0) #comes from above   
    psi1d=np.roll(psi1,-1,axis=0) #comes from below    
    psi1l=np.roll(psi1,1,axis=1) #comes from left
    psi1lu=np.roll(psi1l,1,axis=0) #comes from left and above
    psi1ld=np.roll(psi1l,-1,axis=0) #comes from left and below
    psi1r=np.roll(psi1,-1,axis=1) #comes from right
    psi1ru=np.roll(psi1r,1,axis=0) #comes from right and above
    psi1rd=np.roll(psi1r,-1,axis=0) #comes from right and below
    psi1m=np.zeros_like(psi1)-1
    psi1m[(psi1>psi1u) & (psi1>psi1d) & (psi1>psi1l) & (psi1>psi1r) & (psi1>psi1lu) & (psi1>psi1ld) & (psi1>psi1ru) & (psi1>psi1rd)]=-2 # all max
    psi1m[(psi1<psi1u) & (psi1<psi1d) & (psi1<psi1l) & (psi1<psi1r) & (psi1<psi1lu) & (psi1<psi1ld) & (psi1<psi1ru) & (psi1<psi1rd)]=-3 # all min
    
#        print('1-',np.sum(psi1m==-2))
#        print(np.where(psi1m==-2))
#        print('2-',np.sum(psi1m==-3))
#        print(np.where(psi1m==-3))
    suspect1 = np.where((psi1m==-2) | (psi1m==-3))
    print('suspected points: ', len(suspect1[0]), suspect1)
    suspect=[]
    for i in range(len(suspect1[0])):
        suspect.append([suspect1[0][i],suspect1[1][i],psi1[suspect1[0][i],suspect1[1][i]]])
    suspect=np.asarray(suspect)
    suspect = suspect[(suspect[:, 2]).argsort()]
    countindex = suspect              
    
    zi2=np.zeros_like(zi1)
    zi22=np.zeros([suspect.shape[0], zi1.shape[0], zi1.shape[1]])
    counter=0
#        zi1[np.isnan(u)]=np.nan
        
########### option 2 start
    for k in range(len(suspect)):
        zi2=np.zeros_like(zi1) # test overlap
        counter=k+1
        zi2[int(suspect[k][0]),int(suspect[k][1])]=counter
        print (k, 'CCOOUUNNTTEERR',counter, suspect[k][0], suspect[k][1])
        loop = 0 # debug print variable     
        while True:
            loop +=1
            found = False       
            countindex = np.where(zi2==counter)
#            print('gloop ',counter, loop, len(np.where(zi2==-counter)[0]), len(np.where(zi2==counter)[0]))#,countindex)
#            zi2[int(suspect[k][0]),int(suspect[k][1])]=-counter # test overlap
            for i in range(len(countindex[0])):
                zi2[countindex[0][:],countindex[1][:]]=-counter # test overlap
                if countindex[0][i]>0:
                    if zi2[countindex[0][i]-1,countindex[1][i]]!=-counter: # test overlap
                        if zi1[countindex[0][i]-1,countindex[1][i]]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i]-1,countindex[1][i]]=counter
                                found=True               
                    if countindex[1][i]>0:
                        if zi2[countindex[0][i]-1,countindex[1][i]-1]!=-counter:# test overlap
                            if zi1[countindex[0][i]-1,countindex[1][i]-1]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i]-1,countindex[1][i]-1]=counter
                                found=True
                    if countindex[1][i]<zi1.shape[1]-1:
                        if zi2[countindex[0][i]-1,countindex[1][i]+1]!=-counter:# test overlap
                            if zi1[countindex[0][i]-1,countindex[1][i]+1]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i]-1,countindex[1][i]+1]=counter
                                found=True
                if countindex[0][i]<zi1.shape[0]-1:
                    if zi2[countindex[0][i]+1,countindex[1][i]]!=-counter:# test overlap
                        if zi1[countindex[0][i]+1,countindex[1][i]]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i]+1,countindex[1][i]]=counter
                                found=True
                    if countindex[1][i]>0:
                        if zi2[countindex[0][i]+1,countindex[1][i]-1]!=-counter:# test overlap
                            if zi1[countindex[0][i]+1,countindex[1][i]-1]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i]+1,countindex[1][i]-1]=counter
                                found=True
                    if countindex[1][i]<zi1.shape[1]-1:
                        if zi2[countindex[0][i]+1,countindex[1][i]+1]!=-counter:# test overlap
                            if zi1[countindex[0][i]+1,countindex[1][i]+1]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i]+1,countindex[1][i]+1]=counter
                                found=True
                if countindex[1][i]>0:
                    if zi2[countindex[0][i],countindex[1][i]-1]!=-counter:# test overlap
                        if zi1[countindex[0][i],countindex[1][i]-1]<=zi1[countindex[0][i],countindex[1][i]]:
                            zi2[countindex[0][i],countindex[1][i]-1]=counter
                            found=True
                if countindex[1][i]<zi1.shape[1]-1:
                    if zi2[countindex[0][i],countindex[1][i]+1]!=-counter:# test overlap
                        if zi1[countindex[0][i],countindex[1][i]+1]<=zi1[countindex[0][i],countindex[1][i]]:
                            zi2[countindex[0][i],countindex[1][i]+1]=counter
                            found=True
                                
                                
            if found is False:
                break
        zi22[k,:,:]=zi2[:,:]
    zi2 = np.zeros_like(zi1)
    for i in range(zi2.shape[0]):
        for j in range(zi2.shape[1]):
            if np.sum(zi22[:,i,j]!=0)!=1:
                zi2[i,j]=0
            else:
                zi2[i,j]=np.sum(zi22[:,i,j])
########### option 2 end            
        
#        plt.figure()
#        ax = plt.subplot(1,2,1)
#        plt.imshow(u,origin='lower')
#        plt.colorbar()
#        plt.title('u')
#
##        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
#
#        plt.subplot(1,2,2, sharex=ax, sharey=ax)
#        plt.imshow(w,origin='lower')
#        plt.colorbar()
#        plt.title('w')
#
##        plt.subplot(1,3,3, sharex=ax, sharey=ax)
#        
#        plt.figure()
#        zi1[np.isnan(w)]=np.nan       
#        plt.imshow(zi1,origin='lower')
#        plt.colorbar()
#        plt.title('psi1')
#
#    plt.figure()       
##        zi2[np.isnan(w)]=np.nan
##        zi2[w==0]=np.nan
#    plt.imshow(np.abs(zi2),origin='lower')
#    plt.colorbar()
#    plt.title('flood gradient')

#    plt.figure()       
##        zi2[np.isnan(w)]=np.nan
#    plt.imshow(np.abs(zi22[36,:,:]),origin='lower')
#    plt.colorbar()
#    plt.title('flood gradient 36')
        
    return zi2



def eddies(files, itter, axisindex, axispos, xminc, xmaxc, yminc, ymaxc, zminc, zmaxc, currentRow, parallel=True):
        from skimage import measure
        from shapely.geometry import Point
        from shapely.geometry.polygon import Polygon
        
        separate = False
            
        textboxtime1value = itter
        lu=[]
        lw=[]
        lpsi1m=[]
#        lfuse=[]
        for filesindex in range(len(files)):
            textboxfile1value = files[filesindex] #self.filename1.text()
            print('eddie:',filesindex, files, itter, axisindex, axispos, xminc, xmaxc, yminc, ymaxc, zminc, zmaxc, currentRow)
            db = get_slice(textboxfile1value, textboxtime1value, axisindex, axispos,
                           clipxmin=xminc, clipxmax=xmaxc,
                           clipymin=yminc, clipymax=ymaxc,
                           clipzmin=zminc, clipzmax=zmaxc, parallel=parallel)

            cma = plt.cm.OrRd
            cma.set_bad(alpha=0.0)

            ticks = 600
            xmin, xmax, zmin, ymax, u = makegrid(db, 'U_x', currentRow, ticks=ticks, method='linear')
            xmin, xmax, zmin, ymax, w = makegrid(db, 'U_z', currentRow, ticks=ticks, method='linear')
            mask = np.zeros_like(u)
            mask[u==0]=np.nan
            mask[np.isnan(u)]=np.nan
#            u = np.round(u,5)
#            w = np.round(w,5)
            
            lu.append(u)
            lw.append(w)
            u[u==0]= np.NaN
            w[w==0]= np.NaN
        
            xi = np.linspace(0,ticks,ticks)
            zi = np.linspace(0,ticks,ticks)
            X,Z=np.meshgrid(xi,zi)
            w[np.isnan(w)]=0
            u[np.isnan(u)]=0
                      
            u1 = np.sign(u)
            w1 = np.sign(w)
        
            test0 = np.zeros_like(u)
            test0[(u1==1)&(w1==1)]=1
            test0[(u1==1)&(w1==-1)]=2
            test0[(u1==-1)&(w1==1)]=3
            test0[(u1==-1)&(w1==-1)]=4

#            plt.figure()
#            test0[np.isnan(mask)]=np.nan
#            plt.imshow(test0,origin='lower')
#            plt.title('directions '+ shrinktitle(files[filesindex]))
#            plt.colorbar()
#            
#            plt.figure()                      
#            plt.streamplot(X,Z,u, w, density=15, color=u, linewidth=2, cmap=plt.cm.autumn)
#            plt.colorbar()    
#            plt.title('stramlines')
                                   
# https://stackoverflow.com/questions/49557329/compute-stream-function-from-x-and-y-velocities-by-integration-in-python            
            # integrate to make an intial guess
            intx=integrate.cumtrapz(w,X,axis=1,initial=0)[0]
            inty=integrate.cumtrapz(u,Z,axis=0,initial=0)
            psi1=intx-inty     
            
#            intx=integrate.cumtrapz(w,X,axis=1,initial=0)
#            inty=integrate.cumtrapz(u,Y,axis=0,initial=0)[:,0][:,None]
#            psi2=intx-inty
            
#            psi1=psi2.copy()

#            psi1[np.isnan(unan)]==np.NaN  
#            psi1[u==0]=np.NaN
                 
#            plt.figure()
#            psi1[np.isnan(mask)]=np.nan
#            plt.imshow(psi1, origin='lower')
#            plt.title('psi1 ' + shrinktitle(files[filesindex]))
#            plt.colorbar()
#            plt.show()
#            
#            eddie_mask = eddie_mask_gradient(psi1)
#
#            plt.figure()
#            eddie_mask[np.isnan(mask)]=np.nan
#            plt.imshow(eddie_mask, origin='lower')
#            plt.title('eddie gradient mask' + str(filesindex))
#            plt.colorbar()
#            plt.show()
#            
#            uw = (u*u+w*w)**0.5
#            uwu=np.roll(uw,1,axis=0) #comes from above   
#            uwd=np.roll(uw,-1,axis=0) #comes from below    
#            uwl=np.roll(uw,1,axis=1) #comes from left
#            uwlu=np.roll(uw,1,axis=0) #comes from left and above
#            uwld=np.roll(uw,-1,axis=0) #comes from left and below
#            uwr=np.roll(uw,-1,axis=1) #comes from right
#            uwru=np.roll(uw,1,axis=0) #comes from right and above
#            uwrd=np.roll(uw,-1,axis=0) #comes from right and below
#            uwm=np.zeros_like(u)-1
#            uwm[(uw>uwu) & (uw>uwd) & (uw>uwl) & (uw>uwr) & (uw>uwlu) & (uw>uwld) & (uw>uwru) & (uw>uwrd)]=-2 # all max
#            uwm[(uw<uwu) & (uw<uwd) & (uw<uwl) & (uw<uwr) & (uw<uwlu) & (uw<uwld) & (uw<uwru) & (uw<uwrd)]=-3 # all min
#            plt.figure()
#            uwm[u==0]=np.NaN
#            uwm[np.isnan(mask)]=np.nan
#            plt.imshow(uwm, origin='lower')
#            plt.title('uw minmax')
#            plt.colorbar()
##            plt.streamplot(xx,zz,u, w, density=15, color=u, linewidth=2, cmap=plt.cm.autumn)
#            suspectuw = np.where((uwm==-2) | (uwm==-3))
#            plt.plot(suspectuw[1],suspectuw[0],'r*')
#            plt.show()
#            
            
#            print('before traj')
#            tmp =  self.trajectory(u,w, dx,dz, [20,208], u) #175
#            print('aftertraj')
#            plt.streamplot(xx,zz,u, w, density=15, color=u, linewidth=2, cmap=plt.cm.autumn)
            
#            plt.figure()
#            plt.imshow(u, origin='lower')
#            plt.title('u')
#            plt.show()
#            plt.streamplot(xx,zz,u, w, density=15, color=u, linewidth=2, cmap=plt.cm.autumn)
#   
            psi1u=np.roll(psi1,1,axis=0) #comes from above   
            psi1d=np.roll(psi1,-1,axis=0) #comes from below    
            psi1l=np.roll(psi1,1,axis=1) #comes from left
            psi1lu=np.roll(psi1l,1,axis=0) #comes from left and above
            psi1ld=np.roll(psi1l,-1,axis=0) #comes from left and below
            psi1r=np.roll(psi1,-1,axis=1) #comes from right
            psi1ru=np.roll(psi1r,1,axis=0) #comes from right and above
            psi1rd=np.roll(psi1r,-1,axis=0) #comes from right and below
            psi1m=np.zeros_like(u)-1
            psi1m[(psi1>psi1u) & (psi1>psi1d) & (psi1>psi1l) & (psi1>psi1r) & (psi1>psi1lu) & (psi1>psi1ld) & (psi1>psi1ru) & (psi1>psi1rd)]=-2 # all max
            psi1m[(psi1<psi1u) & (psi1<psi1d) & (psi1<psi1l) & (psi1<psi1r) & (psi1<psi1lu) & (psi1<psi1ld) & (psi1<psi1ru) & (psi1<psi1rd)]=-3 # all min
        
#            psi1m[np.isnan(unan)]==np.NaN                              
#            psi1m[u==0]=np.NaN
            plt.figure(100*(filesindex+1))
            
            print('1-',np.sum(psi1m==-2))
            print(np.where(psi1m==-2))
            print('2-',np.sum(psi1m==-3))
            print(np.where(psi1m==-3))
            suspect = np.where((psi1m==-2) | (psi1m==-3))
            print('suspected points: ', len(suspect[0]), suspect)

            
            
            
            
            
#            uwgrad=np.gradient(uw)[0]          

# loop al suspected point
            vortexlist=[]
#            suspect=[[],[]]
            for s in range(len(suspect[0])):
#                em = -eddie_mask[int(suspect[0][s]), int(suspect[1][s])]
                
                psi1em = psi1.copy()
#                psi1em[eddie_mask!=em] = np.nan
#                psi1em[int(suspect[0][s]), int(suspect[1][s])] = psi1[int(suspect[0][s]), int(suspect[1][s])]
#                print('mask;', em, len(np.isnan(psi1)), len(np.isnan(psi1em)))
# go for two grid above suspect and get it's value
                velocity = w[int(suspect[0][s]), int(suspect[1][s])]**2.+u[int(suspect[0][s]), int(suspect[1][s])]**2.
                print("checking suspected point ",s, suspect[0][s], suspect[1][s], psi1[int(suspect[0][s]), int(suspect[1][s])], velocity)
                if ((int(suspect[0][s])+1>psi1m.shape[0]-1) or (int(suspect[1][s])+1>psi1m.shape[1]-1) or (suspect[1][s]==0) or (suspect[0][s]==0)):  # not enoush space in the domain for contour
                    continue
                radius=2
                found = True
                confull= []
                point = Point(suspect[0][s], suspect[1][s])
                while found:
                    found = False
                    above = interpvel(psi1em,psi1em,int(suspect[0][s])+radius, int(suspect[1][s]))[0]
    #                print('above ', above)
    # get all contours with this value

                    contours = measure.find_contours(psi1em, above)
#                    print('contours', radius, above, len(contours))
                    for c in range(len(contours)):
                        # filter only closed contours (first and last points are the same)
                        if (contours[c][-1]==contours[c][0])[0] and (contours[c][-1]==contours[c][0])[1]: 
                            # check if suspected point inside the contour
                            con = []
                            conx= []
                            cony= []
                            for p in range(len(contours[c])):
                                con.append((contours[c][p][0],contours[c][p][1]))
                                conx.append(contours[c][p][1])
                                cony.append(contours[c][p][0])
                            contouri = Polygon(con)
#                            print(s, c, radius, ' inside? ', contouri.contains(point))                            
                            if contouri.contains(point):
                                # we want to check if there are few points in the contour 
                                count_points = 0
                                for op in range(len(suspect[0])):
                                    otherpoint = Point(suspect[0][op], suspect[1][op])
                                    if contouri.contains(otherpoint):
                                        count_points +=1
                                if separate:
                                    if count_points==1:
                                        found = True
                                        confull = con
                                        radius += .1
                                else:
                                    if count_points>0:
                                        found = True
                                        confull = con
                                        radius += .1
                                    
#                                plt.plot(conx,cony,'r')
# if yes, raise a flag and move another grid above
# plot the maximum contour
# fill the maximum contour
                if len(confull)>0:
                    print('conxcony',s,radius, above, suspect[0][s], suspect[1][s])
                    plt.figure(100*(filesindex+1))
#                    plt.plot(conx,cony,'r')
                    vlist=[]
                    contourifull = Polygon(confull)
                    for i in range(psi1.shape[0]):
                        for j in range(psi1.shape[1]):
                            point = Point(i, j)
                            if contourifull.contains(point):
                                if psi1m[i,j]>=0:
                                    vlist.append(psi1m[i,j])
                                psi1m[i,j]=s

                    # start - if vortex has several centers, we have to choose the one with the lowest velocity
                    vortexlist.append([s,velocity])
                    minvelocity =  velocity
                    minvortex = s
                    vlist=list(set(vlist))
                    print('vlist',vlist)
                    print('vortexlist',vortexlist)
                    for j in range(len(vortexlist)):
                        if vortexlist[j][0] in vlist:
                            if minvelocity>vortexlist[j][1]:
                                minvelocity = vortexlist[j][1]
                                minvortex=vortexlist[j][0]
                    for i in range(len(vlist)):
                        psi1m[psi1m==vlist[i]]=minvortex
                        print('change old vortex ',vlist[i], 'to ',minvortex)
                    if minvortex!=s:
                        psi1m[psi1m==s]=minvortex
                        print('change vortex ',s, 'to ', minvortex)
                    # end - combine several vortecies                            
#                    histogram(psi1[psi1m==s],str(filesindex)+', '+str(s)+' psi')
#                    histogram(uw[psi1m==s],str(filesindex)+', '+str(s)+' uw')
#                    histogram(uwgrad[psi1m==s], str(filesindex)+', '+str(s)+' uwgrad')
#                    corrfunc(uw[psi1m==minvortex],X[psi1m==minvortex],Y[psi1m==minvortex],shrinktitle(textboxfile1value)+str(s)+' uw')           
                        
#            plt.figure()
#            plt.imshow(u, origin='lower')            
#            plt.title('u '+ shrinktitle(files[filesindex]))
#            plt.show()
#            
                    
            lpsi1m.append(psi1m)
            plt.figure(100*(filesindex+1))
            psi1m[u==0]=np.NaN
            psi1m[np.isnan(mask)]=np.NaN
            plt.imshow(psi1m, origin='lower')
            plt.title('psi1 minmax '+ shrinktitle(files[filesindex]))
            plt.colorbar()
#            plt.streamplot(xx,zz,u, w, density=15, color=u, linewidth=2, cmap=plt.cm.autumn)
            plt.plot(suspect[1],suspect[0],'r*')
            plt.show()
#            print('max Y', Y.max())
#            corrfunc(uw[(Y>350) & (X>350)],X[(Y>350) & (X>350)],Y[(Y>350) & (X>350)],shrinktitle(textboxfile1value)+'Y>150 uw')
#            
#            fuse=eddie_mask.copy()
#            fuse[psi1m==-1]=0
#            fuse[np.abs(u)<0.000001]=np.nan
#            fuse[np.isnan(u)]=np.nan
#            lfuse.append(np.abs(fuse))
#            plt.figure()
#            fuse[np.isnan(mask)]=np.nan
#            plt.imshow(fuse, origin='lower')
#            plt.title('fuse '+shrinktitle(files[filesindex]))
#            plt.colorbar()
#            plt.show()
#            
            
            
            print ('fin1')
#
            if 5==6:
            # https://www.nco.ncep.noaa.gov/pmb/products/gfs/ surface analysis gdas.t12z.atmf006.nc   
                filename = '/data3/nirb/gdas.t12z.atmf003.nc'
                filename = '/data3/nirb/gdas.t12z.sfcf003.nc'
                import netCDF4 as nc
                ds = nc.Dataset(filename)
                print(ds['pressfc'])
                pressfc = ds['pressfc'][0]
                pressfc[pressfc<95001]=np.NaN
                print(ds['pressfc'])
                w = ds['vgrd10m'][0]
                u = ds['ugrd10m'][0]
    #            u1=u[200:600,1200:1600]
    #            w1=w[200:600,1200:1600]
    #            xi1=np.linspace(1,400,400)
    #            yi1=np.linspace(1,400,400)
                xi=np.linspace(1,u.shape[1],u.shape[1])
                yi=np.linspace(1,u.shape[0],u.shape[0])
                X,Y=np.meshgrid(xi,yi)
                intx=integrate.cumtrapz(w,X,axis=1,initial=0)[0]
                inty=integrate.cumtrapz(u,Y,axis=0,initial=0)
                psi1=intx-inty
                
 
                intx=integrate.cumtrapz(w,X,axis=1,initial=0)
                inty=integrate.cumtrapz(u,Y,axis=0,initial=0)[:,0][:,None]
                psi2=intx-inty               
                
                psi1 = 0.5*(psi1+psi2)
                
                psi1h=psi1.copy()
                window = 1
                for i in range(0,psi1.shape[0]):
                    print(i,'/',psi1.shape[0])
                    for j in range(0,psi1.shape[1]):
                        tmp = 0
                        for i1 in range(max(-window,-i),min(window+1, psi1.shape[0]-i-1)):
                            for j1 in range(max(-window,-j),min(window+1, psi1.shape[1]-j-1)):
                                tmp += psi1[i+i1,j+j1]
                        psi1h[i,j] = tmp / (window * window)
                psi1=psi1h.copy()
                psi1u=np.roll(psi1,1,axis=0) #comes from above   
                psi1d=np.roll(psi1,-1,axis=0) #comes from below    
                psi1l=np.roll(psi1,1,axis=1) #comes from left
                psi1lu=np.roll(psi1l,1,axis=0) #comes from left and above
                psi1ld=np.roll(psi1l,-1,axis=0) #comes from left and below
                psi1r=np.roll(psi1,-1,axis=1) #comes from right
                psi1ru=np.roll(psi1r,1,axis=0) #comes from right and above
                psi1rd=np.roll(psi1r,-1,axis=0) #comes from right and below
                psi1m=np.zeros_like(u)
                psi1m[(psi1>psi1u) & (psi1>psi1d) & (psi1>psi1l) & (psi1>psi1r) & (psi1>psi1lu) & (psi1>psi1ld) & (psi1>psi1ru) & (psi1>psi1rd)]=1 # all max
                psi1m[(psi1<psi1u) & (psi1<psi1d) & (psi1<psi1l) & (psi1<psi1r) & (psi1<psi1lu) & (psi1<psi1ld) & (psi1<psi1ru) & (psi1<psi1rd)]=2 # all min
                suspect = np.where((psi1m==1))
                print('suspected points: ', len(suspect[0]), suspect)
    
                plt.figure()
                plt.imshow(psi1)
                plt.colorbar()
                plt.plot(suspect[1],suspect[0],'r*')
                plt.title('psi '+ shrinktitle(files[filesindex]))
                plt.show()
                plt.streamplot(X,Y,u, -w, density=15, color=u, linewidth=2, cmap=plt.cm.autumn)
        print('*****files*****')
        matchvortex=0
        matchvortextotal=0
        matchsize=0
        matchsizetotal=0 # all vortices, 
        for i in range(len(files)):
            print(i,files[i])
        if len(files)>1:
            for i in range(len(lu)):
                lpsi1m[i][lpsi1m[i]==-3]=-1 # there can be few centers in one vortex
                lpsi1m[i][lpsi1m[i]==-2]=-1
                uniquevortex=np.unique(lpsi1m[i][~np.isnan(lpsi1m[i])])
#                uniquevortex=np.unique(lfuse[i][~np.isnan(lfuse[i])])
                print('simulation', i, ' unique: ',uniquevortex)                
                for j in range(1,len(uniquevortex)):
                    vsize=np.sum(lpsi1m[i]==uniquevortex[j])
#                    vsize=np.sum(lfuse[i]==uniquevortex[j])
                    print('vortex id', uniquevortex[j] ,'size', vsize,'pos:',np.where(lpsi1m[i]==uniquevortex[j])[0].min(),np.where(lpsi1m[i]==uniquevortex[j])[1].min())
                    if i==0:  # sum the total size of all votices of basic simulation
                        matchsizetotal+=vsize
                    
            #let us assume there are only two simulations
            #few vortices in each simulation
            uniquev0=np.unique(lpsi1m[0][~np.isnan(lpsi1m[0])])
#            uniquev0=np.unique(lfuse[0][~np.isnan(lfuse[0])])
            uniquev0=uniquev0[uniquev0!=-1]
            print('start match')
            luw0=(lu[0]**2.+lw[0]**2.)**0.5
            luw1=(lu[1]**2.+lw[1]**2.)**0.5
            uwgrad0=np.gradient(luw0)[0]
#            plt.figure()
#            plt.imshow(uwgrad0,origin='lower')
#            plt.colorbar()
#            plt.title('gradient')
#            plt.show()
            uwgrad1=np.gradient(luw1)[0]    
            for i in range(0,len(uniquev0)):
                print('*', 'Matched ', matchvortex, 'out of ', matchvortextotal,' *')
                print('size', matchsize,'of', matchsizetotal)
                
                ratiosize = float(np.sum(lpsi1m[0]==uniquev0[i]))/matchsizetotal
#                ratiosize = float(np.sum(lfuse[0]==uniquev0[i]))/matchsizetotal
#                print(np.sum(lpsi1m[0]==uniquev0[i]),matchsizetotal, ratiosize)
#                print(np.sum(lfuse[0]==uniquev0[i]),matchsizetotal, ratiosize)
                if ratiosize>0.01:
                   matchvortextotal+=1 
                else:
                   print('sizematch',uniquev0[i], np.sum(lpsi1m[0]==uniquev0[i]))
                   continue 
                second=lpsi1m[1][lpsi1m[0]==uniquev0[i]]
#                second=lfuse[1][lfuse[0]==uniquev0[i]]
                second=np.asarray(list(set(second)))
                second=second[second!=-1]
                second=second[~np.isnan(second)]        
                if len(second[second!=-1])>0:
                    print('*****try match ',uniquev0[i],'with ',second[second!=-1])   
                good=0
                for j in range(0,len(second)):
                    print('sizematch:',uniquev0[i],second[j], np.sum(lpsi1m[0]==uniquev0[i]),np.sum(lpsi1m[1]==second[j]))
#                    print('sizematch',uniquev0[i],second[j], np.sum(lfuse[0]==uniquev0[i]),np.sum(lfuse[1]==second[j]))
#                    center0=np.where(luw0==luw0[lpsi1m[0]==uniquev0[i]].min())
#                    center1=np.where(luw1==luw1[lpsi1m[1]==second[j]].min())
#                    center0=np.where(luw0==luw0[lfuse[0]==uniquev0[i]].min())
#                    center1=np.where(luw1==luw1[lfuse[1]==second[j]].min())
#                    print('centers:',center0,center1)
                    chiv =  chi(luw0[lpsi1m[0]==uniquev0[i]],luw1[lpsi1m[1]==second[j]])
                    chivg = chi(uwgrad0[lpsi1m[0]==uniquev0[i]],uwgrad1[lpsi1m[1]==second[j]])
#                    chiv =  chi(luw0[lfuse[0]==uniquev0[i]],luw1[lfuse[1]==second[j]])
#                    chivg = chi(uwgrad0[lfuse[0]==uniquev0[i]],uwgrad1[lfuse[1]==second[j]])
                    print('chi square is',chiv)
                    print('chi square of grad is', chivg)
                    dualsize=np.sum((lpsi1m[0]==uniquev0[i]) & (lpsi1m[1]==second[j]))
#                    dualsize=np.sum((lfuse[0]==uniquev0[i]) & (lfuse[1]==second[j]))
                    print('dual size is', dualsize)
                    dualratio=float(dualsize)/(np.sum(lpsi1m[0]==uniquev0[i]))
#                    dualratio=float(dualsize)/max(np.sum(lfuse[0]==uniquev0[i]),np.sum(lfuse[1]==second[j]))
                    if (chiv<1.0) & (chivg<1.):
                        if ratiosize>0.01:
                            if dualratio>0.1:
#                                matchvortex += 1
##                                matchsize+=min(float(np.sum(lpsi1m[0]==uniquev0[i])),np.sum(lpsi1m[1]==second[j]))#dualsize
#                                matchsize+=float(np.sum(lpsi1m[0]==uniquev0[i]))
                                good=float(np.sum(lpsi1m[0]==uniquev0[i]))
                                print(' +++ good match ')
                            else:
                                print('--- bad dual ratio')
                        else:
                            print(' --- bad match ratio size')
                    else:
                        print(' --- bad match chi<0.4')
                    if np.isnan(chi(uwgrad0[lpsi1m[0]==uniquev0[i]],uwgrad1[lpsi1m[1]==second[j]])):
#                    if np.isnan(chi(uwgrad0[lfuse[0]==uniquev0[i]],uwgrad1[lfuse[1]==second[j]])):
                        print('nnaann')
                        print('min1',uwgrad0.min(),uwgrad1.min())
                        plt.figure()
                        plt.imshow(uwgrad0,origin='lower')
                        plt.title('uwgrad0'+ str(filesindex))
                        plt.colorbar()
                        plt.show()
                        plt.figure()
                        plt.imshow(uwgrad1,origin='lower')
                        plt.title('uwgrad1'+ str(filesindex))
                        plt.colorbar()
                        plt.show()
                if good>0:
                    matchvortex += 1
                    matchsize+=dualsize #good
                
            print('profile')
            print('r:', stat(np.mean(lu[0],axis=1),np.mean(lu[1],axis=1),kind='r'), 
                  'r^2:', stat(np.mean(lu[0],axis=1),np.mean(lu[1],axis=1),kind='r2'),
                  'rmse:', stat(np.mean(lu[0],axis=1),np.mean(lu[1],axis=1),kind='rmse'),
                  'rmse div:', stat(np.mean(uwgrad0,axis=1),np.mean(uwgrad1,axis=1),kind='rmse'))
            print('fine:',
                  'mean:',np.mean(lu[0]),
                  'std:', np.std(lu[0]),
                  'mean std:',np.mean(uwgrad0),
                  'std div:', np.std(uwgrad0))
            print('coarse:',
                  'mean:',np.mean(lu[1]),
                  'std:', np.std(lu[1]),
                  'mean div:', np.mean(uwgrad1),
                  'std div:', np.std(uwgrad1))
            maxspeed = np.abs((lu)).max()
            uwr2=stat(np.mean(luw0,axis=1),np.mean(luw1,axis=1),kind='r2')
            print('max diff',np.abs(np.mean(lu[0],axis=1)-np.mean(lu[1],axis=1)).max(), 'maxspeed ', maxspeed)
            print('uwr:', stat(np.mean(luw0,axis=1),np.mean(luw1,axis=1),kind='r'), 
                  'uwr^2:', uwr2,
                  'uwrmse:', stat(np.mean(luw0,axis=1),np.mean(luw1,axis=1),kind='rmse'))
            uwgradr2=stat(np.mean(uwgrad0,axis=1),np.mean(uwgrad1,axis=1),kind='r2')
            print('gradr:', stat(np.mean(uwgrad0,axis=1),np.mean(uwgrad1,axis=1),kind='r'), 
                  'gradr^2:', uwgradr2,
                  'gradrmse:', stat(np.mean(uwgrad0,axis=1),np.mean(uwgrad1,axis=1),kind='rmse'))
            chiuw = chi(luw0.ravel(),luw1.ravel())
            chigraduw = chi(uwgrad0.ravel(),uwgrad1.ravel())
            print('chi', chiuw, 'chigrad', chigraduw)
            nmse = stat(luw0.ravel(),luw1.ravel(),kind='nmse')
            fb = stat(luw0.ravel(),luw1.ravel(),kind='fb')
            print('Mean: NMSE:',nmse, 'FB:',fb)
            nmse = stat(uwgrad0.ravel(),uwgrad1.ravel(),kind='nmse')
            fb = stat(uwgrad0.ravel(),uwgrad1.ravel(),kind='fb')
            print('STD: NMSE:',nmse, 'FB:',fb)
            print('*', 'Matched ', matchvortex, 'out of ', matchvortextotal,' *')
            print('size', matchsize,'of', matchsizetotal)
            print('assume ', files[0], 'is the fine simulation')
            plt.figure()
            for i in range(len(lu)):
#                print(i,np.mean(lu[i],axis=1))
                plt.plot(np.mean(lu[i],axis=1), zi, label=shrinktitle(files[i]))
            plt.plot((np.mean(lu[0],axis=1)-np.mean(lu[1],axis=1)),zi,label='diff')
            plt.title('wind profile')
            plt.legend()
            plt.show()
            if (uwr2>0.0) & (uwgradr2>0.0) & (float(matchvortex)/float(matchvortextotal)>=0.5) & (matchsize/float(matchsizetotal)>=0.5) & (chiuw<1.) & (chigraduw<1.):
                print('Almost the same !!!')
            else:
                print('not the same.....')
#
#            plt.figure()
#            plt.scatter(luw0.ravel(), luw1.ravel(), s=1)
#            line11 = np.linspace(min(luw0.min(),luw1.min()), max(luw0.max(),luw1.max()))
#            plt.scatter(line11, line11, s=1, c='r')
#            plt.xlabel('0')
#            plt.ylabel('1')
#            plt.title('uw')
#            plt.show()
#            
#            plt.figure()
#            plt.scatter(uwgrad0.ravel(), uwgrad1.ravel(), s=1)
#            line11 = np.linspace(min(uwgrad0.min(),uwgrad1.min()), max(uwgrad0.max(),uwgrad1.max()))
#            plt.scatter(line11, line11, s=1, c='r')
#            plt.xlabel('0')
#            plt.ylabel('1')
#            plt.title('div')
#            plt.show()
            
        print ('fin3')
        # fin eddies


def diff_files(files, itter, field, xminc, xmaxc, yminc, ymaxc, zminc, zmaxc, axisindex, axispos, currentRow, heightcontour=False, parallel=True):
        print('diffdebug',parallel)
        textboxtime1value = itter[0]
        textboxtime2value = itter[-1]
        ticks=600
        textboxfile1value = files[0]          
        textboxfile2value = files[1]          
        if heightcontour:
            db = get_slice_height(textboxfile1value, textboxtime1value, axisindex, axispos,
                                  clipxmin=xminc, clipxmax=xmaxc,
                                  clipymin=yminc, clipymax=ymaxc,
                                  clipzmin=zminc, clipzmax=zmaxc)
            dbb = get_slice_height(textboxfile2value, textboxtime2value, axisindex, axispos,
                                  clipxmin=xminc, clipxmax=xmaxc,
                                  clipymin=yminc, clipymax=ymaxc,
                                  clipzmin=zminc, clipzmax=zmaxc)
        else:
            db = get_slice(textboxfile1value, textboxtime1value, axisindex, axispos,
                                  clipxmin=xminc, clipxmax=xmaxc,
                                  clipymin=yminc, clipymax=ymaxc,
                                  clipzmin=zminc, clipzmax=zmaxc, parallel=parallel)
            print('so far2')

            dbb = get_slice(textboxfile2value, textboxtime2value, axisindex, axispos,
                                  clipxmin=xminc, clipxmax=xmaxc,
                                  clipymin=yminc, clipymax=ymaxc,
                                  clipzmin=zminc, clipzmax=zmaxc, parallel=parallel)
        xmin1, xmax1, ymin1, ymax1, zi1 = makegrid(db , field, currentRow, ticks=ticks, method='linear')
        xmin2, xmax2, ymin2, ymax2, zi2 = makegrid(dbb, field, currentRow, ticks=ticks, method='linear')
        
        x1=np.linspace(xmin1,xmax1,ticks)
        x2=np.linspace(xmin2,xmax2,ticks)
        y1=np.linspace(ymin1,ymax1,ticks)
        y2=np.linspace(ymin2,ymax2,ticks)
        
        xmin = max(xmin1, xmin2)
        xmax = min(xmax1, xmax2)
        ymin = max(ymin1, ymin2)
        ymax = min(ymax1, ymax2)
        
        xi1min=np.argmax(x1>=xmin)
        xi1max=np.argmin(x1<xmax)
        xi2min=np.argmax(x2>=xmin)
        xi2max=np.argmin(x2<xmax)
        yi1min=np.argmax(y1>=ymin)
        yi1max=np.argmin(y1<ymax)
        yi2min=np.argmax(y2>=ymin)
        yi2max=np.argmin(y2<ymax)
        
        print(xmin1,xmin2,xmin)
        print(xmax1,xmax2,xmax)
        print(ymin1,ymin2,ymin)
        print(ymax1,ymax2,ymax)
#        print(x1)
#        print(x2)
#        print(y1)
#        print(y2)
        print(xi1min,xi1max,xi2min,xi2max,yi1min,yi1max,yi2min,yi2max)
        
        zmin = min(np.nanmin(zi1), np.nanmin(zi2))
        zmax = max(np.nanmax(zi1), np.nanmax(zi2))
        z1 = np.linspace(zmin,zmax,ticks)
        
        # Plot density map.
        plt.figure()
        ax1 = plt.subplot(221)
        plt.imshow(
            zi1, extent=(xmin1, xmax1, ymin1, ymax1), origin='lower',
            cmap=plt.get_cmap('GnBu_r'))
#        plt.tricontourf(db['x'], db['y'], db[self.listfields.currentItem().text()], 100)  # , cmap=cma

        if ~np.isnan(zmin):
            plt.clim(zmin, zmax)
        plt.title('db) '+shrinktitle(textboxfile1value)+',Cells:'+str(len(db))+', Time:'+logreader(textboxfile1value, textboxtime1value))
        plt.colorbar()

        # Plot density map.
        plt.subplot(222, sharex=ax1, sharey=ax1)
        plt.imshow(
            zi2, extent=(xmin2, xmax2, ymin2, ymax2), origin='lower',
            cmap=plt.get_cmap('GnBu_r'))
#        plt.tricontourf(dbb['x'], dbb['y'], dbb[self.listfields.currentItem().text()], 100)  # , cmap=cma
        
        if ~np.isnan(zmin):
            plt.clim(zmin, zmax)
        plt.title('dbb) '+shrinktitle(textboxfile2value)+',Cells:'+str(len(dbb))+', Time:'+logreader(textboxfile2value, textboxtime2value))
        plt.colorbar()
        
        
        
#        db = get_slice(textboxfile1value, textboxtime1value, axisindex, axispos,
#                           clipxmin=xmin, clipxmax=xmax,
#                           clipymin=ymin, clipymax=ymax,
#                           clipzmin=zmin, clipzmax=zmax, reader=self.reader1)
#        dbb = get_slice(textboxfile2value, textboxtime2value, axisindex, axispos,
#                            clipxmin=xmin, clipxmax=xmax,
#                            clipymin=ymin, clipymax=ymax,
#                            clipzmin=zmin, clipzmax=zmax)
#        xmin1, xmax1, ymin1, ymax1, zi1 = makegrid(db, self.listfields.currentItem().text(), self.listaxis.currentRow(), ticks=ticks)
#        xmin2, xmax2, ymin2, ymax2, zi2 = makegrid(dbb, self.listfields.currentItem().text(), self.listaxis.currentRow(), ticks=ticks)
#        
# 
        
        
        
        plt.subplot(223, sharex=ax1, sharey=ax1)
        plt.imshow(
            zi2-zi1, extent=(xmin, xmax, ymin, ymax), origin='lower',
            cmap=plt.get_cmap('GnBu_r'))
        plt.title('dbb-db')
        plt.colorbar()
        plt.subplot(224)
        plt.scatter(
            zi2, zi1, s=1)
        line11 = np.linspace(zmin, zmax)
        plt.scatter(line11, line11, s=1, c='r')
        plt.xlabel(shrinktitle(textboxfile2value))
        plt.ylabel(shrinktitle(textboxfile1value))
        zi1nan = np.isnan(zi1.ravel())
        zi2nan = np.isnan(zi2.ravel())

        corr = round(np.corrcoef(zi1.ravel()[~(zi2nan | zi1nan)], zi2.ravel()[~(zi2nan | zi1nan)])[1, 0],4)
        # zi10 = zi1.ravel()[~(zi2nan | zi1nan)]
        # zi20 = zi2.ravel()[~(zi2nan | zi1nan)]
        # zi20 = zi20[zi10!=0]
        # zi10 = zi10[zi10!=0]
        plt.title(field+': corr'+str(corr))

        plt.show()
        field1 = zi1.ravel()[~(zi2nan | zi1nan)]
        field2 = zi2.ravel()[~(zi2nan | zi1nan)]
        rmseval = rmse(field1, field2, zeros=False)
        print('corr rmse:',corr, rmseval)
        print('min:', field1.min(),field2.min())
        print('max:', field1.max(),field2.max())
        print('mean:', field1.mean(),field2.mean())
        print('mean diff:', (field1-field2).mean())
        print('nmse:', stat(field1,field2,'nmse'))
        print('rmse:', stat(field1,field2,'rmse'))
        print('r:', stat(field1,field2,'r'))
        print('r2:', stat(field1,field2,'r2'))


def vprofile(files, bounds, itter):
    plt.figure()
    for i in range(len(files)):
        db = get_clip(files[i], itter,
          clipxmin=bounds[0], clipxmax=bounds[1],
          clipymin=bounds[2], clipymax=bounds[3],
          clipzmin=0., clipzmax=0.)

        print('len(db) is', len(db))
        ap3 = db[['x', 'y', 'z', 'U_x']]  
        points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
        valuesux = ap3['U_x']
        xmin=ap3['x'].min() #cordxmin
        xmax=ap3['x'].max() #cordxmax
        ymin=ap3['y'].min() #cordymin
        ymax=ap3['y'].max() #cordymax
        zmin=ap3['z'].min()
        zmax=ap3['z'].max()
        # ticks = 600j
        grid_points = 42000000.  # add 00
        meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)
        resolution = (meters / grid_points)**(1./3)
        if resolution==0:
            print('no clip boundaries in one of the axis')
        # resolution = 0.1
        ticksx = int((xmax-xmin) / resolution / 2) * 1j 
        ticksy = int((ymax-ymin) / resolution / 2) * 1j 
        ticksz = int((zmax-zmin) / resolution * 4) * 1j 
        xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
        gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticksz.imag))
       
        windfine = []
        for k in range(int(ticksz.imag)):
            windfine.append(gridux[:,:,k][gridux[:,:,k]!=0.].mean())
#                        windmax.append(gridux[:,:,i][gridux[:,:,i]!=0.].max())
        plt.plot(windfine,z,label='windfine ' + shrinktitle(files[i]))
    plt.legend()
    plt.ylabel("z [m]")
    plt.xlabel('U [m/s]')
    plt.title('Ux vertical profile')
#        plt.ylim(0,800)
    plt.show()


def plot_file(files, fields, itter, xminc, xmaxc, yminc, ymaxc, zminc, zmaxc, axisindex, axispos, addtotitle="", parallel=True):

    if 5==7:
        files=[r'/ibdata2/nirb/openFOAM/vortex/tlv1/tlv1.foam']
        files=[u'/ibdata2/nirb/openFOAM/porous/windAroundBuildings-davidsona/windAroundBuildings-davidsona.OpenFOAM',
               u'/ibdata2/nirb/openFOAM/porous/windAroundBuildings-davidson0a/windAroundCube.foam']
        axisindex = 1
        axispos = 664200
        xminc=177750
        xmaxc=179700
        yminc=0
        ymaxc=0
        zminc=0
        zmaxc=40
        currentRow=1  
        itter=444444
        fields=['U_x']
        filesindex=0
        
        
        
        itter=40000
        files = [ r'/data4bk/nirb/Simulations/Mala/mala2b/mala2b.foam',
        ]
        areaname='mala'
        textboxtime1value = 300000.0
        axisindex = 1
        axis_index = 1
        axispos = 563735.0
        axis_value = 563735.0
        learnfrom = ''
        currentitem =  u'U_x'
        fields = [u'U_x']
        xminc = 0.0
        xmaxc = 0.0
        yminc = 0.0
        ymaxc = 0.0
        zminc = 0.0
        zmaxc = 0.0
        currentRow=1
        xmin=155200
        xmax=155400
        ymin=0
        ymax=0
        zmin=0
        zmax=0
        parallel=True    
        addtotitle=''


    for filesindex in range(len(files)):
        textboxfile1value = files[filesindex] #self.filename1.text()
        db = get_slice(textboxfile1value, itter, axisindex, axispos,
                       clipxmin=float(xminc), clipxmax=float(xmaxc),
                       clipymin=float(yminc), clipymax=float(ymaxc),
                       clipzmin=float(zminc), clipzmax=float(zmaxc), parallel=parallel)

        cma = plt.cm.OrRd
        cma.set_bad(alpha=0.0)
#         plt.figure()
#         for f in range(len(fields)):
#             if f==0:
#                 ax=plt.subplot(1,len(fields),f+1)
#             else:
#                 plt.subplot(1,len(fields),f+1, sharex=ax, sharey=ax)
#             if axisindex == 0:
# #                plt.tricontourf(db['y'], db['z'], db[self.listfields.currentItem().text()], 100)
#                 plt.tricontourf(db['y'], db['z'], db[fields[f]], 100)
#             if axisindex == 1:
# #                plt.tricontourf(db['x'], db['z'], db[self.listfields.currentItem().text()], 100)
#                 plt.tricontourf(db['x'], db['z'], db[fields[f]], 100)
#             if axisindex == 2:
# #                plt.tricontourf(db['x'], db['y'], db[self.listfields.currentItem().text()], 100)  # , cmap=cma
#                 plt.tricontourf(db['x'], db['y'], db[fields[f]], 100)  # , cmap=cma
#             plt.colorbar()
#             plt.title(addtotitle+shrinktitle(textboxfile1value)+'-'+fields[f])
#             plt.show()


        plt.figure()
        for f in range(len(fields)):
            if f==0:
                ax=plt.subplot(1,len(fields),f+1)
            else:
                plt.subplot(1,len(fields),f+1, sharex=ax, sharey=ax)

    #        plt.imshow(zi1, origin='lower', cmap=plt.get_cmap('tab20c'))  # jet, Paired
            xmin, xmax, ymin, ymax, zi1 = makegrid(db, fields[f], axisindex, ticks=5000, method='linear')
            # Plot density map.
            zi1[zi1 == 0] = np.NaN
    #        print('lendb',len(db))
    #        print(logreader(textboxfile1value, textboxtime1value))
            plt.imshow(zi1, origin='lower')  # jet, Paired
            plt.title(addtotitle+shrinktitle(textboxfile1value)+'-'+fields[f])
            plt.colorbar()
        plt.show()

        plt.figure()
        for f in range(len(fields)):
            if f==0:
                ax=plt.subplot(1,len(fields),f+1)
            else:
                plt.subplot(1,len(fields),f+1, sharex=ax, sharey=ax)
            if axisindex == 0:
                plt.scatter(db['y'], db['z'], c=db[fields[f]], s=4)
            if axisindex == 1:
                plt.scatter(db['x'], db['z'], c=db[fields[f]], s=4)
            if axisindex == 2:
                plt.scatter(db['x'], db['y'], c=db[fields[f]], s=4)
            plt.title(addtotitle+shrinktitle(textboxfile1value)+'-'+fields[f])
            # plt.colorbar()
            plt.show()
        
#        ffvalue = np.abs(np.fft.fftshift(np.fft.fft2(wi)))
#        print('zi2',ffvalue)
#        print('222')
#        
#        plt.figure()
##        plt.imshow(zi1, origin='lower', cmap=plt.get_cmap('tab20c'))  # jet, Paired
#        plt.imshow(ffvalue, origin='lower', vmin=0, vmax=200)  # jet, Paired
#        plt.title('fft')
#        plt.colorbar()
#        plt.show()

    
    # plt.figure()
    # plt.imshow(db['U_x'], origin='lower', cmap=plt.get_cmap('tab20c'))  # jet, Paired
    # # plt.title('db) Cells:'+str(len(db))+', Time:'+logreader(textboxfile1value, textboxtime1value))
    # plt.colorbar()
    # plt.show()


def profile3(files,textboxtime1value,axisindex,axispos,learnfrom,areaname,startarea,endarea,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui, verbose=False, parallel=True):
    
    print(files,textboxtime1value,axisindex,axispos,learnfrom,startarea,endarea,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui)

#        files=[r'/ibdata2/nirb/openFOAM/porous/hadasflatcoarse8/windAroundCube.foam']
#        textboxtime1value=900
#        files=[r'/ibdata2/nirb/openFOAM/porous/sinthetic0/sinthetic0.OpenFOAM']
#        textboxtime1value=1700

    hadassin=[[2000,	3000,	2000,	3000,	0.25,	0.45,	50.]]
    hadassin=[[2000,	3000,	2000,	3000,	0.25,	0.25,	50.]]
    hadassin=[[2000,	3000,	2000,	3000,	0.25,	0.25,	100.]]
    hadassin=[[2000,	3000,	2000,	3000,	0.08,	0.08,	100.]]

    hadassin=[[100,	150,	0,	150,	0.08,	0.08,	2.5]]

    hadassin=[[1,	2,	1,	2,	0.25,	0.25,	0.1]]
    hadassin=[[200,	400,	200,	400,	0.25,	0.25,	10]]
    hadassin=[[-100000,	1000000,	-1000000,	1000000,	0.25,	0.25,	0.10]]
    hadassin44=[[-100000,	1000000,	-1000000,	1000000,	0.44,	0.44,	0.10]]
    hadassin44big=[[-100000,	1000000,	-1000000,	1000000,	0.44,	0.44,	10.]]
    hadassin16=[[-100000,	1000000,	-1000000,	1000000,	0.16,	0.16,	0.10]]
    hadassin6=[[-100000,	1000000,	-1000000,	1000000,	0.06,	0.06,	0.10]]
    hadassinlien=[[-100000,	1000000,	-1000000,	1000000,	0.25,	0.25,	0.15]]

    
    locx=179000
    locy=664000
    radius = 900
    
    cordxmin=178987
    cordxmax=179387
    cordymin=663894
    cordymax=664294


#        hadas=hadas200
    hadas=np.asarray(choosearea('jer200'))
    hadas=np.asarray(choosearea('hadas200'))
    hadas=np.asarray(choosearea('yehuda200'))
    hadas=np.asarray(choosearea(areaname))

#        hadas=np.asarray(hadassin)

    mac006=[
            [0.25,3.5],
            [0.45,4.9],
            [0.60,5.4],
            [0.75,5.8],
            [0.95,5.9],
            [1.15,5.95],
            [1.35,6.2],
            [1.50,6.4],
            [1.65,7],
            [1.85,7.9],
            [1.95,8.6],
            [2.25,9.25],
            [2.55,10.],
            [2.90,10.20],
            [3.30,10.75],
            [3.55,11.35],
            [4.25,12.15],
            [4.95,12.65],
            [5.70,13.3],
            [6.40,13.65],
            [7.10,14.25],
            [8.85,15.40],
            [10.55,15.75],
            [14,16.6],
            ]
    mac006=np.asarray(mac006)
    mac006[:,0] = mac006[:,0]*4./7
    mac006[:,1] = mac006[:,1]/9.#*17.56
    
    
    
    mac016=[
            [0.10,0.6],
            [0.30,2.3],
            [0.45,2.75],
            [0.50,2.95],
            [0.70,3.05],
            [1.30,3.6],
            [1.50,4.35],
            [1.70,5.55],
            [1.85,7.05],
            [1.95,7.45],
            [2.30,8.1],
            [2.80,9.6],
            [3.25,10.35],
            [3.95,10.75],
            [4.50,11.35],
            [5.20,12.05],
            [5.85,12.3],
            [6.50,12.65],
            [8.25,13.6],
            [9.85,14.15],
            ]
    mac016=np.asarray(mac016)
    mac016[:,0] = mac016[:,0]*4./6.5
    mac016[:,1] = mac016[:,1]*1./9#*17.56

    mac044=[
            [1.25,1.10],
            [1.35,1.50],
            [1.55,1.75],
            [1.70,2.30],
            [1.90,6.40],
            [2.05,8.30],
            [2.45,9.30],
            [2.80,10.0],
            [3.15,10.3],
            [3.50,11.2],
            [4.15,11.9],
            [4.80,12.4],
            [5.45,13.2],
            [6.15,13.5],
            [6.90,13.9],
            [8.65,15.1],
            [10.4,15.7],
            [13.9,16.3]
            ]
    mac044=np.asarray(mac044)
    mac044[:,0] = mac044[:,0]*4./7.
    mac044[:,1] = mac044[:,1]/9.5#*14.1

    mac044big=[
            [1.25,1.10],
            [1.35,1.50],
            [1.55,1.75],
            [1.70,2.30],
            [1.90,6.40],
            [2.05,8.30],
            [2.45,9.30],
            [2.80,10.0],
            [3.15,10.3],
            [3.50,11.2],
            [4.15,11.9],
            [4.80,12.4],
            [5.45,13.2],
            [6.15,13.5],
            [6.90,13.9],
            [8.65,15.1],
            [10.4,15.7],
            [13.9,16.3]
            ]
    mac044big=np.asarray(mac044big)
    mac044big[:,0] = mac044big[:,0]*4./7.
    mac044big[:,1] = mac044big[:,1]/9.5#*14.1

    lien=[ #drag model 1
            [0.0,0.4],
            [0.08,0.6],
            [0.11,1.0],
            [0.15,1.3],
            [0.20,2.0],
            [0.30,2.9],
            [0.40,3.65],
            [0.50,3.9],
            [0.60,4.2],
            [0.70,4.2],
            [0.80,4.25],
            ]

    lien=[ # RANS
            [0.0,-0.2],
            [0.08,0.4],
            [0.11,0.7],
            [0.15,1.0],
            [0.20,2.5],
            [0.30,3.4],
            [0.40,3.8],
            [0.50,4.0],
            [0.60,4.2],
            [0.70,4.2],
            [0.80,4.3],
            ]
                 
    lien=np.asarray(lien)
    lien[:,0] = lien[:,0]*10./1.5
    lien[:,1] = lien[:,1]/1.#*14.1
    
    if (learnfrom==u'6'):
        mac044=mac006
        hadas=np.asarray(hadassin6)
    if (learnfrom==u'16'):
        mac044=mac016
        hadas=np.asarray(hadassin16)
    if (learnfrom==u'44'):
        mac044=mac044
        hadas=np.asarray(hadassin44)
    if (learnfrom==u'big'):
        mac044=mac044big
        hadas=np.asarray(hadassin44big)
    if (learnfrom==u'lien'):
        mac044=lien
        hadas=np.asarray(hadassinlien)
    
    if learnfrom!=u'':
        heindex = np.argwhere(mac044[:,0]>1.)[0][0]
        uh = mac044[heindex-1, 1] + (1-mac044[heindex-1., 0])*(mac044[heindex, 1]-mac044[heindex-1, 1])/(mac044[heindex, 0]-mac044[heindex-1, 0])
        
    if startarea==endarea:
        startarea=0 #41
        endarea=len(hadas) # 42
    plots=float(endarea-startarea)
    plotsx=math.ceil(plots**0.5)
    plotsy=math.ceil(plots/plotsx)
    readers=[]
    for fileindex in range(len(files)):
        reader = bse.ReadCase(files[fileindex], files[fileindex], CaseType='Decomposed Case')  # ' Reconstructed Decomposed Case')
        readers.append(reader)
    if verbose:
        plt.figure()
    print('proflie31',startarea,endarea)
    file1 = open(areaname+"scores.txt", "a")  # append mode
    file2 = open(areaname+"scores2.txt", "a")  # append mode
    for i in range(len(files)):
        file1.writelines(files[i]+'\n')
        file2.writelines(files[i]+'\n')
    file1.writelines('len areas (without limits) is '+str(len(hadas))+'\n')
    file2.writelines('len areas (without limits) is '+str(len(hadas))+'\n')
    file1.close()
    file2.close()
    
    scores=np.zeros(len(files)+1)
    scoresvalue=np.zeros(len(files)+1)  
    
    for area in range(startarea,min(endarea, len(hadas))): #range(len(hadas)): #41-42
#        for area in range(len(hadas)): #range(len(hadas)):
        print('new area that is checked:', area, ': ',hadas[area][0], ' ',hadas[area][1], ' ',hadas[area][2], ' ',hadas[area][3],' ',hadas[area][4],' ',hadas[area][5],' ',hadas[area][6])
        if verbose:
            plt.subplot(plotsx, plotsy, area-startarea+1)
        comparez = []
        compareu = []
        for fileindex in range(len(files)):
            print('now working on file ', files[fileindex] )
            db = get_clip(files[fileindex], [textboxtime1value],
              clipxmin=xmingui, clipxmax=xmaxgui,
              clipymin=ymingui, clipymax=ymaxgui,
              clipzmin=zmingui, clipzmax=zmaxgui,
              reader=readers[fileindex])

            cordxmin=hadas[area][0]
            cordxmax=hadas[area][1]
            cordymin=hadas[area][2]
            cordymax=hadas[area][3]
            locx=(cordxmax+cordxmin)/2
            locy=(cordymax+cordymin)/2
            radius = (cordxmax-cordxmin)/2
              
#            he, lambdaf, lambdap, lc = morphology(files[0],textboxtime1value, locx,locy,radius=radius)

            lambdaf = hadas[area][4] #0.43
            lambdap = hadas[area][5] #0.33
            he = hadas[area][6] #14.8


            dbsmall=db[0].copy()
            
            print('0:', len(dbsmall))
            dbsmall=dbsmall[dbsmall['x']>cordxmin]
            print('1:', len(dbsmall))
            dbsmall=dbsmall[dbsmall['x']<cordxmax]
            print('2:', len(dbsmall))
            dbsmall=dbsmall[dbsmall['y']>cordymin]
            print('3:', len(dbsmall))
            dbsmall=dbsmall[dbsmall['y']<cordymax]
            print('4', len(dbsmall))
#            print(dbsmall['y'].min(),dbsmall['y'].max())
            ap3 = dbsmall[['x', 'y', 'z', currentitem]]  
#                print('profile3debug0')
#                print('ap3',ap3)
            if len(ap3)==0:
                print ('len(ap3)==0')
                comparez.append([])
                compareu.append([np.nan])                
                continue

            points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
            valuesux = ap3[currentitem]
            xmin=ap3['x'].min() #cordxmin
            xmax=ap3['x'].max() #cordxmax
            ymin=ap3['y'].min() #cordymin
            ymax=ap3['y'].max() #cordymax
            zmin=ap3['z'].min()
            zmax=ap3['z'].max()
            # ticks = 600j
            grid_points = 42000000.  # add 00
            meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)
            resolution = (meters / grid_points)**(1./3)
            if resolution==0:
                print('no clip boundaries in one of the axis',meters,xmax,xmin,ymax,ymin,zmax,zmin)
                xmin = cordxmin
                xmax = cordxmax
                ymin = cordymin
                ymax = cordymax
                zmin = 0
                meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)
                resolution = 50.
                print('no clip2', xmin, xmax, ymin, ymax)
                
            # resolution = 0.1
            ticksx = max(int((xmax-xmin) / resolution / 2) ,1) * 1j
            ticksy = max(int((ymax-ymin) / resolution / 2) ,1) * 1j
            ticksz = max(int((zmax-zmin) / resolution * 4) ,1) * 1j
#            print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)
            xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
#                print ('debug00')
            gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            x = np.linspace(ap3['x'].min(), ap3['x'].max(), int(ticksx.imag))
            y = np.linspace(ap3['y'].min(), ap3['y'].max(), int(ticksy.imag))
            z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticksz.imag))

            zi = zi

            if learnfrom!=u'':
                morphindex = len(mac044[:, 1])-3
                morphindex = np.argwhere(mac044[:, 0]>3.)[0][0]

                print('morphindex ',morphindex)
                print('morphindex1 ',mac044[morphindex, 0])
#                print('morphindex ',zi[0,0,:])

                hemean = np.argwhere(zi[0,0,:]>mac044[morphindex, 0])[0][0]
                print('hemean ',hemean)

                windfinehemean = gridux[:,:,hemean][gridux[:,:,hemean]!=0.].mean()
                print('windfinehemean ',windfinehemean)
            windfine = []
#                windmax  = []

            for i in range(int(ticksz.imag)):
                if learnfrom!=u'':
                    windfine.append(gridux[:,:,i][gridux[:,:,i]!=0.].mean()*mac044[morphindex, 1]/windfinehemean)
                else:
                    windfine.append(gridux[:,:,i][gridux[:,:,i]!=0.].mean().mean())
#                print('zi',zi[0,0,:])

            z11=0
            for i in range(len(zi[0,0,:])):
                if he>zi[0,0,i]:
                    z11=i   
                    
#                z11 = int(min(len(z)*0.8,z11*3))
#                z11=100
            z11 = int(len(z)*0.91)
        #######################################################                
#                plt.plot(windfine,z,label=shrinktitle(files[fileindex]))

            zpoint, indices = np.unique(points3[:,2], return_inverse=True)
            upoint = np.zeros_like(zpoint)
            npoint = np.zeros_like(zpoint)
            for i in range(len(zpoint)):
                upoint[i]=np.mean(valuesux[indices==i])
                npoint[i]=np.sum(indices==i)
            zpoint2 = np.linspace(z.min(),z.max())

            npoint2 = np.zeros_like(zpoint2)
            upoint2 = np.zeros_like(zpoint2)
            dzpoint = (zpoint2[2]-zpoint2[1])/2.0
            for i in range(len(zpoint2)):
                for j in range(len(zpoint)):
                    if zpoint2[i]-dzpoint<zpoint[j] and zpoint2[i]+dzpoint>zpoint[j]:
                        upoint2[i] += npoint[j]*upoint[j]
                        npoint2[i] += npoint[j]
            zpoint3=[]
            upoint3=[]
            for i in range(len(npoint2)):
                if npoint2[i]>0:
                    zpoint3.append(zpoint2[i])
                    upoint3.append(upoint2[i]/npoint2[i])

            zpoint3 = np.asarray(zpoint3)
            upoint3 = np.asarray(upoint3)
#                print(he, 'zpoint3', zpoint3)
            comparez.append(z)
#            print('debugz1',z)
#            print('debugcomparez0',comparez[0])
            compareu.append(windfine)
#            print('debugcompareu',compareu)
            
            if learnfrom!=u'':
                hemean = np.argwhere(zpoint3>mac044[morphindex, 0]/1.)[0][0]  ### big big delete after use (/10.)
                upoint3=upoint3*mac044[morphindex, 1]/upoint3[hemean]
                      

                errorlower = 0
                errorupper = 0
                lowercount = 0
                uppercount = 0

                for i in range(len(mac044[:, 0])):
                    hemean = np.argwhere(zi[0,0,:]>mac044[i, 0])[0][0]
                    
                    if mac044[i, 0]<1:
                        errorlower += ((windfine[hemean] - mac044[i, 1])**2.)**0.5
                        lowercount += 1
                    else:
                        if mac044[i, 0]<3.:
                            errorupper += ((windfine[hemean] - mac044[i, 1])**2.)**0.5
                            uppercount += 1
                
                errorlower /= lowercount
                errorupper /= uppercount

                if verbose:
                    plt.plot(windfine,z,label='windfine ' + shrinktitle(files[fileindex])+ ' ' + str(round(errorlower,4))+ ', '+str(round(errorupper,4)))
            else:           
                print('windfine ' + shrinktitle(files[fileindex]) + ' p='+str(round(lambdap,2))+ ' f='+str(round(lambdaf,2))+ ' he='+str(round(he,2)))
#                plt.plot(upoint3,zpoint3,label='means ' + shrinktitle(files[fileindex])+ ' ' + str(errorlower)+ ', '+str(errorupper))
                if verbose:
                    plt.plot(windfine,z-z[0],label='windfine ' + shrinktitle(files[fileindex]) + ' p='+str(round(lambdap,2))+ ' f='+str(round(lambdaf,2))+ ' he='+str(round(he,2)))
#                    plt.plot(windmax,z-z[0],label='windmax ' + shrinktitle(files[fileindex]))

#                plt.plot(upoint,zpoint,label='points '+ shrinktitle(files[fileindex]))
#                gra = np.gradient(upoint)
#                gra = np.gradient(gra)
#                gra = np.sign(gra) + 3.
#                plt.plot(gra,zpoint,label='gradient sign '+ shrinktitle(files[fileindex]))

        print('profile stat')
        foundall=True
        for i in range(len(comparez)):
            if len(comparez[i])<2:
                foundall=False
        if foundall is True:
            ground=[]
            for i in range(len(comparez)):
                j=0
                while (np.isnan(compareu[i][j])) and (j<len(compareu[i])-1):
                    j+=1
                ground.append(j)
    
            if len(comparez)>0:
                if len(comparez[0])>0:
                    heindex = np.abs(comparez[0]-he).argmin()
                    he3index = np.abs(comparez[0]-(comparez[0][ground[0]]+he*3.)).argmin()
                    if he3index==ground[0]:
                        he3index=min(int(ground[0]+len(comparez[0])/20.),len(comparez[0])/2.)
            dduz=[]
            for i in range(len(comparez)):
                dduz.append([])
    
            for i in range(len(comparez)):
                for j in range(len(comparez[0])):
                    place = np.abs((comparez[i])-comparez[0][j]-(comparez[i][ground[i]]-comparez[0][ground[0]])).argmin()
                    value = compareu[i][place]
                    dduz[i].append(value)
    
            dduz=np.asarray(dduz)
            for i in range(len(comparez)):
#                print('dbgr2 grnd, he3index',i, ground[i], he3index, comparez[0][:he3index],dduz[i][:he3index])
                hhr2=stat(dduz[0][:he3index],dduz[i][:he3index],'r2', nantozero=False)
                print('hh r2 r nmse:', i,hhr2)
                scoresvalue[i]+=hhr2
                if (hhr2>0.3):
                       scores[i]+=1
            logexpu=np.zeros(he3index)
    #        print('logdebug0',ground[0],he3index)
            for j in range(0,he3index):
                logexpu[j]=windprofile(comparez[0][j]-comparez[0][ground[0]], uref=windfine[len(zi[0,0,:])-3], href=zi[0,0,len(zi[0,0,:])-3]-ground[0], he=he, lambdap=lambdap, lambdaf=lambdaf, beta=0.2, verbose=False)
            logexpr2=stat(dduz[0][:he3index],logexpu[:he3index],'r2', nantozero=False)
#            print('dbgr2','logexp', len(dduz[0][:he3index]), dduz[0][:he3index],logexpu[:he3index])
#            print('dbgr2a',comparez[0], comparez[0][j]-comparez[0][ground[0]])
            print('logexp r2=',logexpr2)
            scoresvalue[-1]+=logexpr2
            if (logexpr2>0.3):
                   scores[-1]+=1
                
        if len(ap3)>0:

            z = np.asarray(z)

            if learnfrom==u'':
                windmorphologytest0=[]
                ztest = np.linspace(0,z[-1],z[-1]+1)
                hearea = len(zi[0,0,:])-3

                for i in range(len(ztest)):
                    windmorphologytest0.append(windprofile(ztest[i], uref=windfine[hearea], href=zi[0,0,hearea], he=he, lambdap=lambdap, lambdaf=lambdaf, beta=0.2, verbose=False))
                if verbose:
                    plt.plot(windmorphologytest0, ztest, label='log-exp')
                    plt.plot(np.linspace(0,windfine[hearea],10), np.ones(10)*he)               
            
            if learnfrom!=u'':
                windmorphologytest0=[]
                ztest = np.linspace(0,1,11)
                print('exptest',ztest,mac044[morphindex, 1],mac044[morphindex, 0])
                for i in range(11):
                    windmorphologytest0.append(windprofile(ztest[i], uref=5, href=100, he=he, lambdap=lambdap, lambdaf=lambdaf, beta=0.3))
                if verbose:
                    plt.plot(windmorphologytest0, ztest, label='log-exp')
                
                windmorphology=[]
                for i in range(len(mac044[:, 0])):
                    windmorphology.append(windprofile(mac044[i, 0], uref=mac044[morphindex, 1], href=mac044[morphindex, 0], he=1., lambdap=lambdap, lambdaf=lambdaf, beta=0.3))
                
    #            ax2 = ax1.twinx()
                print(windmorphology)
                errorlower = 0
                errorupper = 0
                lowercount = 0
                uppercount = 0
                for i in range(len(mac044[:, 0])):
                    if mac044[i, 0]<1:
                        errorlower += ((windmorphology[i] - mac044[i, 1])**2.)**0.5
                        lowercount += 1
                    else:
                        if mac044[i, 0]<3.:
                            errorupper += ((windmorphology[i] - mac044[i, 1])**2.)**0.5
                            uppercount += 1
                    
                errorlower /= lowercount
                errorupper /= uppercount
                
                if verbose:
                    plt.plot(windmorphology, mac044[:, 0], label='morphology ' + str(round(errorlower,4))+ ', '+str(round(errorupper,4)))
            if verbose:
                plt.title(currentitem+', area '+str(area)+': '+str(hadas[area][0])+ ' '+str(hadas[area][1])+ ' '+str(hadas[area][2])+ ' '+str(hadas[area][3]))
#           z=np.linspace(5,495,(490/5)+1), he=50, plt.plot(3*0.3/0.41*np.log(0.3*z/0.3),z/he,label='z0=0.3')
            heindex=0
            for i in range(len(z)):
                if z[i]<he:
                    heindex=i
#            plt.plot(np.linspace(0,5,10), (z[heindex]+z[0])*np.ones(10)/he)
            if learnfrom!=u'':
                if verbose:                
                    plt.plot(np.linspace(0,2,10), np.ones(10))               
                    plt.plot(mac044[:,1],mac044[:,0],'*')
            print('**scores=', area, scores, scoresvalue)
            
            file1 = open(areaname+"scores.txt", "a")  # append mode
            file1.writelines(str(area)+" "+ ",".join(str(scores))+"\n")
            file1.close()
            file2 = open(areaname+"scores2.txt", "a")  # append mode
            file2.writelines(str(area)+" "+ ",".join(str(scores))+"\n")
            file2.close()
            
            if verbose:
                plt.legend()
                plt.xlabel('U [m/s]')
                plt.ylabel('z [m]')
#                plt.ylim(0,800)
                plt.show()
    
    print('****scores final=', scores, scoresvalue)
    
    def exponential2(x, a, b, c):
          return 2*np.exp(x/(2*c*c*50*(1-c)/b))
    def exponential(x, a, b):
          return a*np.exp(b*x)
    def log(x, a, b, c):
          return a*np.log((x-c)/b)
#                  return a*np.log(b*x*x-c)**2
    def power_law(x, a, b):
          return a*np.power(x, b)
#            parsexp2, covexp2 = curve_fit(f=exponential2, xdata=z[:10], ydata=windfine[:10], p0=[0, 0], bounds=(-np.inf, np.inf))
#            print('pars', parsexp2)
#            print('cov', covexp2)
#            plt.plot(exponential2(z[:10], *parsexp2),z[:10], linestyle='--', linewidth=2, color='black', label = 'exp2')
#            parsexp, covexp = curve_fit(f=exponential, xdata=z[10:], ydata=windfine[10:], p0=[0, 0], bounds=(-np.inf, np.inf))
#                parslog, covlog = curve_fit(f=log, xdata=z[10:], ydata=windfine[10:], p0=[0.1, 0.1,0.1], bounds=(0, np.inf))
#                parspow, covpow = curve_fit(f=power_law,   xdata=z[10:], ydata=windfine[10:], p0=[0, 0], bounds=(-np.inf, np.inf))
#                reslog = np.mean(windfine[10:] - log(z[10:], *parslog))
#                respow = np.mean(windfine[10:] - power_law(z[10:], *parspow))
#                pol = np.polyfit(z[10:], windfine[10:], 5)
#                plt.plot(log(z[10:], *parslog),z[10:], linestyle='--', linewidth=2, color='black', label = 'log')
#            plt.plot(exponential(z[10:], *parsexp),z[10:], linestyle='--', linewidth=2, color='black', label = 'exp')
#            plt.plot(power_law(z[10:], *parspow),z[10:], linestyle='--', linewidth=2, color='blue', label='power law')
#                p=np.poly1d(pol)
#                respol = np.mean(windfine[10:] -p(z[10:]))     
#            plt.plot(p(z[10:]),z[10:], linestyle='--', linewidth=2, color='red', label='polinom')

        
#            plt.figure()
#            plt.plot(mac044[:,1],mac044[:,0],'*')
#            plt.plot(windmorphology, z/he, label='morphology')
#            plt.plot(upoint3,zpoint3/he)
#            plt.legend()
#            plt.show()
    
           
    
      
    
def profile4(files,itter,axisindex,axispos,learnfrom,areaname,startarea,endarea,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui, verbose=False, sub=False, parallel=True):
    
    if 5==4:
        itter=40000
        files = [ r'/data4bk/nirb/Simulations/michaelstadtfloor1ml/windAroundcaseE.foam',
        ]
        areaname='michaelstadt'
        textboxtime1value = 300000.0
        axisindex = 1
        axispos = -200.0
        learnfrom = ''
        currentitem =  u'U_x'
        xmingui = 0.0
        ymingui = 0.0
        zmingui = 0.0
        xmaxgui = 0.0
        ymaxgui = 0.0
        zmaxgui = 0.0
        currentRow=1
        xmin=625
        xmax=750
        ymin=0
        ymax=0
        zmin=0
        zmax=20
        startarea=0
        endarea=1
        verbose=True
        sub=False
        area=0
        parallel=True
        
        
    print(files,itter,axisindex,axispos,learnfrom,startarea,endarea,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui)

    hadas=np.asarray(choosearea(areaname))

    if startarea==endarea:
        startarea=0 #41
        endarea=len(hadas) # 42
    plots=float(endarea-startarea)
    plotsx=math.ceil(plots**0.5)
    plotsy=math.ceil(plots/plotsx)
    readers=[]
    for fileindex in range(len(files)):
        if usehera:
            # dirname=r'/data4bk/nirb/Simulations/Mala/mala2b/'
            bse = paraviewOpenFOAM(shrinkdir(files[fileindex]))
            # bse = paraviewOpenFOAM(dirname)
            reader=bse.reader
            readers.append(reader)
        else:
            if parallel:
                CaseType='Decomposed Case'
            else:
                CaseType='Reconstructed Case'
            reader = bse.ReadCase(files[fileindex], files[fileindex], CaseType=CaseType)  # ' Reconstructed Decomposed Case')
            readers.append(reader)
    if verbose:
        plt.figure()
    print('proflie4',startarea,endarea)
    file1 = open(areaname+"scores.txt", "a")  # append mode
    file2 = open(areaname+"scores2.txt", "a")  # append mode
    for i in range(len(files)):
        file1.writelines(files[i]+'\n')
        file2.writelines(files[i]+'\n')
    file1.writelines('len areas (without limits) is '+str(len(hadas))+'\n')
    file1.close()
    file2.writelines('len areas (without limits) is '+str(len(hadas))+'\n')
    file2.close()
    if sub is True:
        file3 = open(areaname+"scoressub.txt", "a")  # append mode
        file3.writelines(files[i]+'\n')
        file3.writelines('len areas (without limits) is '+str(len(hadas))+'\n')
        file3.close()
        scoressub=np.zeros(len(files)+1)
        file4 = open(areaname+"scoressub2.txt", "a")  # append mode
        file4.writelines(files[i]+'\n')
        file4.writelines('len areas (without limits) is '+str(len(hadas))+'\n')
        file4.close()
        scoressub2=np.zeros(len(files)+1)
    
    scores=np.zeros(len(files)+1)
    scoresvalue=np.zeros(len(files)+1)    
    
    endarea = min(endarea,len(hadas))
    bigarea=False
    if (endarea-startarea)>0:
        bigarea=True
        db = get_clip(files[fileindex], [itter],
          clipxmin=xmingui, clipxmax=xmaxgui,
          clipymin=ymingui, clipymax=ymaxgui,
          clipzmin=zmingui, clipzmax=zmaxgui,
          reader=readers[0])
        if usehera:
            ap3 = db[['x', 'y', 'z', currentitem]]
        else:
            ap3 = db[0][['x', 'y', 'z', currentitem]]  
        xmintotal=ap3['x'].min() #cordxmin
        ymintotal=ap3['y'].min() #cordxmin
        zmintotal=ap3['z'].min() #cordxmin
        xmaxtotal=ap3['x'].max() #cordxmin
        ymaxtotal=ap3['y'].max() #cordxmin
        zmaxtotal=ap3['z'].max() #cordxmin
        
    
    for area in range(startarea,endarea): #range(len(hadas)): #41-42
#        for area in range(len(hadas)): #range(len(hadas)):
        print('new area that is checked:', area, ': ',hadas[area][0], ' ',hadas[area][1], ' ',hadas[area][2], ' ',hadas[area][3],' ',hadas[area][4],' ',hadas[area][5],' ',hadas[area][6])
        cordxmin=hadas[area][0]
        cordxmax=hadas[area][1]
        cordymin=hadas[area][2]
        cordymax=hadas[area][3]
        lambdaf = hadas[area][4] #0.43
        lambdap = hadas[area][5] #0.33
        he = hadas[area][6] #14.8
        print('coordif',xmintotal,cordxmin,xmaxtotal,cordxmax,ymintotal,cordymin,ymaxtotal,cordymax)
        if xmintotal<=cordxmin and xmaxtotal>=cordxmax and ymintotal<=cordymin and ymaxtotal>=cordymax:
            print('if-in-coord')
            if verbose:
                plt.subplot(plotsx, plotsy, area-startarea+1)
            comparez = []
            compareu = []
            xi=np.asarray([0])
            for fileindex in range(len(files)):
                print('now working on file ', files[fileindex] )
                db = get_clip(files[fileindex], [itter],
                  clipxmin=xmingui, clipxmax=xmaxgui,
                  clipymin=ymingui, clipymax=ymaxgui,
                  clipzmin=zmingui, clipzmax=zmaxgui,
                  reader=readers[fileindex])
                    
                if usehera:
                    dbsmall=db.copy()                   
                else:
                    dbsmall=db[0].copy()
                
                print('0:', len(dbsmall))
                dbsmall=dbsmall[dbsmall['x']>cordxmin]
                print('1:', len(dbsmall))
                dbsmall=dbsmall[dbsmall['x']<cordxmax]
                print('2:', len(dbsmall))
                dbsmall=dbsmall[dbsmall['y']>cordymin]
                print('3:', len(dbsmall))
                dbsmall=dbsmall[dbsmall['y']<cordymax]
                print('4', len(dbsmall))
                ap3 = dbsmall[['x', 'y', 'z', currentitem]]  
                if len(ap3)==0:
                    print ('len(ap3)==0')
                    comparez.append([])
                    compareu.append([np.nan])                
                    continue
                points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
                valuesux = ap3[currentitem]
                if 5==5: #np.asarray(xi).sum()==0: # only for the first available file # if fileindex==0:
                    xmin=ap3['x'].min() #cordxmin
                    xmax=ap3['x'].max() #cordxmax
                    ymin=ap3['y'].min() #cordymin
                    ymax=ap3['y'].max() #cordymax
                    zmin=ap3['z'].min()
                    zmax=ap3['z'].max()
                    # ticks = 600j
                    grid_points = 42000000.  # add 00
                    meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)
                    resolution = (meters / grid_points)**(1./3)
                    if resolution==0:
                        print('no clip boundaries in one of the axis',meters,xmax,xmin,ymax,ymin,zmax,zmin)
                        xmin = cordxmin
                        xmax = cordxmax
                        ymin = cordymin
                        ymax = cordymax
                        zmin = 0
                        meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)
                        resolution = 50.
                        print('no clip2', xmin, xmax, ymin, ymax)
                    
                    # resolution = 0.1
                    print('boundaries',xmax,xmin,ymax,ymin,zmax,zmin)

                    ticksx = max(int((xmax-xmin) / resolution / 2) ,1) * 1j
                    ticksy = max(int((ymax-ymin) / resolution / 2) ,1) * 1j
                    ticksz = max(int((zmax-zmin) / resolution * 4) ,1) * 1j
                    print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)
        #            print('minmax',xmin,xmax,ymin,ymax,zmin,zmax)
        #            x = np.linspace(ap3['x'].min(), ap3['x'].max(), int(ticksx.imag))
        #            y = np.linspace(ap3['y'].min(), ap3['y'].max(), int(ticksy.imag))
                    z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticksz.imag))
    
                    xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
                gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                                   method='nearest')  # Nearest for keeping zeros that are buildings
                print ('debug01')
    
                windfine = []
                windfinestd = []
                for i in range(int(ticksz.imag)):
                    windfine.append(gridux[:,:,i][gridux[:,:,i]!=0.].mean().mean())
                    windfinestd.append(gridux[:,:,i][gridux[:,:,i]!=0.].std())
    
                comparez.append(np.asarray(z))
                compareu.append(np.asarray(windfine))
#                if verbose:                
#                    for i in range(40):
#                        print('>>>',i,comparez[-1][i],compareu[-1][i])
                
                
                print('windfine ' + shrinktitle(files[fileindex]) + ' p='+str(round(lambdap,2))+ ' f='+str(round(lambdaf,2))+ ' he='+str(round(he,2)))
#                if verbose:
#                    plt.plot(windfine,z-z[0],label='windfine ' + shrinktitle(files[fileindex]) + ' p='+str(round(lambdap,2))+ ' f='+str(round(lambdaf,2))+ ' he='+str(round(he,2)))
    
            print('profile stat')
            foundall=True
            for i in range(len(comparez)):
                if len(comparez[i])<2:
                    foundall=False
            if foundall is True:
                ground=[]
                for i in range(len(comparez)):
                    j=0
                    while (np.isnan(compareu[i][j])) and (j<len(compareu[i])-1):
                        j+=1
                    ground.append(j)

                zground=[]
                uground=[]
                ugroundstd=[]                
                uground0=[]  # interpolated to base simulation with relative to ground height
                for i in range(len(comparez)):
                    zground.append([])
                    uground.append([])
                    ugroundstd.append([])
                    uground0.append([])
                    for j in range(len(comparez[i])):
                        if j>=ground[i]:
                            if ground[i]==0:
                                trueground=comparez[i][ground[i]] #0.0
                            else:
                                trueground = comparez[i][ground[i]]
                            zground[i].append(comparez[i][j]-trueground)
                            uground[i].append(compareu[i][j])
                            ugroundstd[i].append(compareu[i][j])
#                    print('interplog',i, len(zground[i]), len(uground[i]),zground[i],uground[i])
                    uground0[i]=np.interp(zground[0],zground[i],uground[i])
                            
        
                if len(comparez)>0:
                    print('lencomparez',len(comparez),len(comparez[0]), he )
                    if len(comparez[0])>0:
    #                    heindex = np.abs(comparez[0]-he).argmin()
                        he3index = np.abs(zground[0]-(he*3.)).argmin()
                        if he==0:
                            he3index=int(len(uground[0])/2.)
                        he3index = np.abs(zground[0]-(meanheight(areaname)*3.)).argmin()    
                        he3index = int(len(uground[0])/10.)   
                for i in range(len(comparez)):
    #                print('dbgr2 grnd, he3index',i, ground[i], he3index, comparez[0][:he3index],dduz[i][:he3index])
                    hhr2=stat(uground0[0][:he3index],uground0[i],'r2', nantozero=False, verbose=str(i))
                    if sub:
                        hhr2sub=stat(uground0[1][:he3index],uground0[i][:he3index],'r2', nantozero=False, verbose=str(i))
                    fb=stat(uground0[0][:he3index],uground0[i][:he3index],'fb', nantozero=False)
                    nmse=stat(uground0[0][:he3index],uground0[i][:he3index],'nmse', nantozero=False)
                    rmse=stat(uground0[0][:he3index],uground0[i][:he3index],'rmse', nantozero=False)
                    std = np.std(uground0[0][:he3index])
                    print('hh r2 nmse fb rmse/std, rmse, std:', i,hhr2, nmse, fb, rmse/std, rmse, std)
                    scoresvalue[i]+=hhr2
                    if ((hhr2>0.3) & (nmse<4.) & (math.fabs(fb)<0.3) & (rmse<=std)):
                           scores[i]+=1
                    if sub:
                        if (hhr2sub>0.3):
                               scoressub[i]+=1
                        scoressub2[i]+=hhr2sub
                        
                    if verbose:
#                        for k in range(len(uground[i])):
#                            print('malai',i,k, zground[i][k],uground[i][k])
                        plt.plot(uground[i],zground[i],label='windfine ' + shrinktitle(files[i]) + ' p='+str(round(lambdap,2))+ ' f='+str(round(lambdaf,2))+ ' he='+str(round(he,2)))
                        plt.plot(uground[i]+np.asarray(windfinestd),zground[i],label='windfine++ ' + shrinktitle(files[i]) + ' p='+str(round(lambdap,2))+ ' f='+str(round(lambdaf,2))+ ' he='+str(round(he,2)))
                        plt.plot(uground[i]-np.asarray(windfinestd),zground[i],label='windfine-- ' + shrinktitle(files[i]) + ' p='+str(round(lambdap,2))+ ' f='+str(round(lambdaf,2))+ ' he='+str(round(he,2)))
#                        print('minmax',min(zground[i]),max(zground[i]),min(uground[i]),max(uground[i]))
                logexpu=np.zeros(he3index)
        #        print('logdebug0',ground[0],he3index)
                for j in range(0,he3index):
#                    logexpu[j]=windprofile(zground[0][j], uref=windfine[len(zi[0,0,:])-3], href=zi[0,0,len(zi[0,0,:])-3]-comparezground[0], he=he, lambdap=lambdap, lambdaf=lambdaf, beta=0.2, verbose=False)
                    logexpu[j]=windprofile(zground[0][j], uref=uground[0][-3], href=zground[0][-3],he=he, lambdap=lambdap, lambdaf=lambdaf, beta=0.2, verbose=False)
                logexpr2=stat(uground0[0][:he3index],logexpu[:he3index],'r2', nantozero=False)
                if sub:
                    logexpr2sub=stat(uground0[1][:he3index],logexpu,'r2', nantozero=False)
                    if (logexpr2sub>0.3):
                           scoressub[-1]+=1
                    scoressub2[-1]+=logexpr2sub
                print('logexp r2=',logexpr2)
                scoresvalue[-1]+=logexpr2
                if (logexpr2>0.3):
                       scores[-1]+=1
                   
            if len(ap3)>0:
                z = np.asarray(z)
    
                if learnfrom==u'':
                    windmorphologytest0=[]
                    ztest = np.linspace(0,z[-1],int(z[-1]+1))
                    hearea = len(zi[0,0,:])-3
    
                    for i in range(len(ztest)):
                        windmorphologytest0.append(windprofile(ztest[i], uref=windfine[hearea], href=zi[0,0,hearea], he=he, lambdap=lambdap, lambdaf=lambdaf, beta=0.2, verbose=False))
                    if verbose:
                        plt.plot(windmorphologytest0, ztest, label='log-exp')
                        plt.plot(np.linspace(0,windfine[hearea],10), np.ones(10)*he)               
                
                if verbose:
                    plt.title(currentitem+', area '+str(area)+': '+str(hadas[area][0])+ ' '+str(hadas[area][1])+ ' '+str(hadas[area][2])+ ' '+str(hadas[area][3]))
    
                print('**scores=', area, scores, scoresvalue)
                
                file1 = open(areaname+"scores.txt", "a")  # append mode
                file1.writelines(str(area)+" "+ ",".join(str(scores))+"\n")
                file1.close()
                
                file2 = open(areaname+"scores2.txt", "a")  # append mode
                file2.writelines(str(area)+" "+ ",".join(str(scoresvalue))+"\n")
                file2.close()
                
                if sub:
                    file3 = open(areaname+"scoressub.txt", "a")  # append mode
                    file3.writelines(str(area)+" "+ ",".join(str(scoressub))+"\n")
                    file3.close()
                    file4 = open(areaname+"scoressub2.txt", "a")  # append mode
                    file4.writelines(str(area)+" "+ ",".join(str(scoressub2))+"\n")
                    file4.close()
                    
                
                if verbose:
                    plt.legend()
                    plt.xlabel('U [m/s]')
                    plt.ylabel('z [m]')
    #                plt.ylim(0,800)
                    plt.show()
    print('****scores final=', scores, scoresvalue)


def areaaccuracy(areaname):
    areaname='ashkelon'
#    areaname='natanya'
    areaname='bs'
#    areaname='yehuda200'
#    areaname='tlvbig250'
    
    with open(areaname+'scores.txt', 'r') as f:
        lines = f.read().splitlines()
    with open(areaname+'scores2.txt', 'r') as f:
        lines2 = f.read().splitlines()
    start = len(lines)
    start2 = len(lines2)
    while start>0:
        start -= 1
        if lines[start][:3]=='len'    :
            break
    print(lines[start])

    while start2>0:
        start2 -= 1
        if lines2[start2][:3]=='len'    :
            break
    print(lines2[start2])
    
    # count how many simulations we have - start
    lastline = lines[-1]
    lastline = lastline.replace(",", "")
    lastline = lastline.replace("[", "")
    lastline = lastline.replace("]", "")
    simulations = len(lastline.split())-1
    # count how many simulations we have - end
    # count how many simulations we have - start
    lastline2 = lines2[-1]
    lastline2 = lastline2.replace(",", "")
    lastline2 = lastline2.replace("[", "")
    lastline2 = lastline2.replace("]", "")
    simulations2 = len(lastline2.split())-1
    # count how many simulations we have - end
    
    # add columns for simulations - start
    table = choosearea(areaname)
    table = np.asarray(table)
    table1 = np.zeros((table.shape[0],table.shape[1]+simulations))
    table1[:,:-simulations] = table
    table1[:,table.shape[1]:]=np.nan
    table2 = np.zeros((table.shape[0],table.shape[1]+simulations))
    table2[:,:-simulations] = table
    table2[:,table.shape[1]:]=np.nan
    # add columns for simulations - end

    for i in range(start+1, len(lines)):
        if i>start+1:
            lastline = lines[i-1]
            lastline = lastline.replace(",", "")
            lastline = lastline.replace("[", "")
            lastline = lastline.replace("]", "")
            line1 = lastline.split()
        else:
            line1 = np.zeros(simulations+1)  # the +1 is for the area number
        lastline = lines[i]
        lastline = lastline.replace(",", "")
        lastline = lastline.replace("[", "")
        lastline = lastline.replace("]", "")
        line2 = lastline.split()
        for j in range(1, simulations+1):
            table1[int(float(line2[0]))-1,6+j] = float(line2[j])-float(line1[j])

    for i in range(start2+1, len(lines2)):
        if i>start2+1:
            lastline2 = lines2[i-1]
            lastline2 = lastline2.replace(",", "")
            lastline2 = lastline2.replace("[", "")
            lastline2 = lastline2.replace("]", "")
            line2 = lastline2.split()
        else:
            line2 = np.zeros(simulations2+1)  # the +1 is for the area number
            
        lastline2 = lines2[i]
        lastline2 = lastline2.replace(",", "")
        lastline2 = lastline2.replace("[", "")
        lastline2 = lastline2.replace("]", "")
        line22 = lastline2.split()
        for j in range(1, simulations2+1):
            table2[int(float(line22[0]))-1,6+j] = float(line22[j])-float(line2[j])

# 6 is height, 4,5 is lambda
    print(np.nansum(table1[table1[:,4]>=0], axis=0)[7:]) # height == 0
    print(np.nansum(table1[table1[:,4]>0.1], axis=0)[7:])
    print(np.nansum(table1[table1[:,4]==0], axis=0)[7:]) # height == 0
#    print(np.nansum(table1[table1[:,6]>2.1], axis=0)[7:])
#    print(np.nansum(table1[(table1[:,4]>0.0) & (table1[:,4]<0.1)], axis=0)[7:])

    print(np.nansum(table1[table1[:,4]>=0], axis=0)[7:]/np.nansum(table1[table1[:,4]>=0], axis=0)[7])
    print(np.nansum(table1[table1[:,4]>0.1], axis=0)[7:]/np.nansum(table1[table1[:,4]>0.1], axis=0)[7])
    print(np.nansum(table1[table1[:,4]==0], axis=0)[7:]/np.nansum(table1[table1[:,4]==0], axis=0)[7]) # height == 0
#    print(np.nansum(table1[table1[:,6]>2.1], axis=0)[7:]/np.nansum(table1[table1[:,6]>=2.1], axis=0)[7])
#    print(np.nansum(table1[(table1[:,4]>0.0) & (table1[:,4]<0.1)], axis=0)[7:]/np.nansum(table1[(table1[:,4]>0.0) & (table1[:,4]<0.1)], axis=0)[7])
    
    print(np.nansum(table2[table1[:,4]>=0], axis=0)[7:]) # height == 0
    print(np.nansum(table2[table1[:,4]>0.1], axis=0)[7:])
    print(np.nansum(table2[table1[:,4]==0], axis=0)[7:]) # height == 0

    print(np.nansum(table2[table1[:,4]>=0], axis=0)[7:]/np.nansum(table2[table1[:,4]>=0], axis=0)[7])
    print(np.nansum(table2[table1[:,4]>0.1], axis=0)[7:]/np.nansum(table2[table1[:,4]>0.1], axis=0)[7])
    print(np.nansum(table2[table1[:,4]==0], axis=0)[7:]/np.nansum(table2[table1[:,4]==0], axis=0)[7]) # height == 0
    
    
    import pandas as pd
    df = pd.DataFrame(table, columns = ['x','xmax','y','ymax','lambdaf','lambdap','he'])
    df['accuracy']=table1[:,10]
    xmin, xmax, ymin, ymax, zilf = makegrid(df, 'lambdaf', 2, ticks=5000, method='nearest')
    xmin, xmax, ymin, ymax, zilp = makegrid(df, 'lambdap', 2, ticks=5000, method='nearest')
    xmin, xmax, ymin, ymax, zilh = makegrid(df, 'he', 2, ticks=5000, method='nearest')
    xmin, xmax, ymin, ymax, zila = makegrid(df, 'accuracy', 2, ticks=5000, method='nearest')
    
    plt.figure()
    ax = plt.subplot(2,2,1)
    plt.imshow(zilp, origin='lower', interpolation='nearest')  # jet, Paired
    plt.colorbar()
    plt.title('$\lambda_f$')
    plt.subplot(2,2,2, sharex=ax, sharey=ax)
    zfshift=np.roll(zilp,125,axis=1)
    zfshift=zfshift-zilp
    zfshift[zfshift>3.]=3.
    plt.imshow(zfshift, origin='lower', interpolation='nearest')  # jet, Paired
    plt.colorbar()
    plt.title('$\lambda_p$')
    plt.subplot(2,2,3, sharex=ax, sharey=ax)
    plt.imshow(zilh, origin='lower', interpolation='nearest')  # jet, Paired ,extent=[xmin,xmax,ymin,ymax]
    plt.colorbar()
    plt.title('He')
    plt.subplot(2,2,4, sharex=ax, sharey=ax)
    plt.imshow(zila, origin='lower', interpolation='nearest')  # jet, Paired
    plt.colorbar()
    plt.title('Accuracy')
    
    plt.figure()
    plt.scatter(zilf.ravel(),zila.ravel())
    plt.xlabel('lambdaf')
    plt.ylabel('accuracy')
    plt.figure()
    plt.scatter(zilp.ravel(),zila.ravel())
    plt.xlabel('lambdap')
    plt.ylabel('accuracy')
    plt.figure()
    plt.scatter(zilh.ravel(),zila.ravel())
    plt.xlabel('zh')
    plt.ylabel('accuracy')
    
    
    return

def compare():
    files=[r'/data4bk/nirb/Simulations/michaelstadtfloor1ml/michaelstadtfloor1ml.foam']
    itter=[1,2300]
    files=[r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor1ml/michaelstadtfloor1ml.foam']
    itter=[1,6000]
    files=[r'/data4bk/nirb/Simulations/Dans/tlvsmallml/tlvs2.foam']
    itter=[1,5000]
    xmin=-600
    xmax=600
    ymin=-400
    ymax=400
    zmin=0
    zmax=60

    xmin=-300
    xmax=300
    ymin=-200
    ymax=200
    zmin=0
    zmax=30
    db = get_clip(files[0], [itter[0]],
                   clipxmin=xmin, clipxmax=xmax,
                   clipymin=ymin, clipymax=ymax,
                   clipzmin=zmin, clipzmax=zmax) 

    ap3 = db[['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)            
    points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
    # values = ap3[sel[0].data()]
    valuesux = ap3['U_x']
    valuesuy = ap3['U_y']
    valuesuz = ap3['U_z']
    if xmin==xmax:
        xmin=ap3['x'].min()
        xmax=ap3['x'].max()
    if ymin==ymax:
        ymin=ap3['y'].min()
        ymax=ap3['y'].max()
    if zmin==zmax:
        zmin=ap3['z'].min()
        zmax=ap3['z'].max()
    # ticks = 600j
    grid_points = 27000000
    meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)

    resolution = (meters / grid_points)**(1./3)
    if resolution==0:
        print('no clip boundaries in one of the axis')
    # resolution = 0.1
    ticksx = int((xmax-xmin) / resolution) * 1j
    ticksy = int((ymax-ymin) / resolution) * 1j
    ticksz = int((zmax-zmin) / resolution) * 1j
    print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)

    xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
    gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                       method='nearest')  # Nearest for keeping zeros that are buildings
    griduy = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuy), (xi, yi, zi),
                       method='nearest')  # Nearest for keeping zeros that are buildings
    griduz = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuz), (xi, yi, zi),
                       method='nearest')  # Nearest for keeping zeros that are buildings

    db2 = get_clip(files[0], [itter[1]],
                   clipxmin=xmin, clipxmax=xmax,
                   clipymin=ymin, clipymax=ymax,
                   clipzmin=zmin, clipzmax=zmax) 

    ap32 = db2[['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)            
    points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
    # values = ap3[sel[0].data()]
    valuesux2 = ap32['U_x']
    valuesuy2 = ap32['U_y']
    valuesuz2 = ap32['U_z']
    gridux2 = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux2), (xi, yi, zi),
                       method='nearest')  # Nearest for keeping zeros that are buildings
    griduy2 = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuy2), (xi, yi, zi),
                       method='nearest')  # Nearest for keeping zeros that are buildings
    griduz2 = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuz2), (xi, yi, zi),
                       method='nearest')  # Nearest for keeping zeros that are buildings

    gu=(gridux.ravel()**2.+griduy.ravel()**2.)**.5
    gu2=(gridux2.ravel()**2.+griduy2.ravel()**2.)**.5
    plt.figure()
    plt.scatter(gu2,gu, s=0.1)    
    plt.xlabel('RANS [m/s]')
    plt.ylabel('ML [m/s]')

    plt.figure()
    plt.scatter(griduy2.ravel(),griduy.ravel(), s=0.01)    
    plt.xlabel('RANS [m/s]')
    plt.ylabel('ML [m/s]')



###################################################3333
    print('start to read U file in the time directory')
    slash = files[0].rfind("/") + 1
    directory = files[0][:slash]+ str(int(float(itter[0]))) +"/"
    fileu = directory + "U"
    filex = directory + "Cx"
    filey = directory + "Cy"
    filez = directory + "Cz"    
    fhu = open(fileu, 'r')
    fhx = open(filex, 'r')
    fhy = open(filey, 'r')
    fhz = open(filez, 'r')    
    txtu = fhu.read()
    txtx = fhx.read()
    txty = fhy.read()
    txtz = fhz.read()    
    txtusplit = txtu.split()
    txtxsplit = txtx.split()  
    txtysplit = txty.split() 
    txtzsplit = txtz.split()    
    posu = txtusplit.index("internalField")
    posx = txtxsplit.index("internalField")
    posy = txtysplit.index("internalField") 
    posz = txtzsplit.index("internalField")    
    items = int(txtusplit[posu+3])
    us = []
    xs = []
    ys = []
    zs = []    
    for i in range(items):
        us.append([float(remove_left(txtusplit[posu+5+i*3])), float(txtusplit[posu+6+i*3]), float(remove_right(txtusplit[posu+7+i*3]))])
        xs.append(float(txtxsplit[posx+5+i])) 
        ys.append(float(txtysplit[posy + 5 + i])) 
        zs.append(float(txtzsplit[posz + 5 + i]))
    xs = np.asarray(xs)
    ys = np.asarray(ys)
    zs = np.asarray(zs) 
    clip = np.zeros_like(xs)  
    clip[clip==0]=True
    clip[xs<xmin]=False
    clip[xs>xmax]=False
    clip[ys<ymin]=False
    clip[ys>ymax]=False
    clip[zs<zmin]=False
    clip[zs>zmax]=False
        
    directory = files[0][:slash]+ str(int(float(itter[1]))) +"/"
    fileu = directory + "U"
    fhu = open(fileu, 'r')
    txtu = fhu.read()
    txtusplit = txtu.split()
    posu = txtusplit.index("internalField")
    items = int(txtusplit[posu+3])
    us2 = []
    for i in range(items):
        us2.append([float(remove_left(txtusplit[posu+5+i*3])), float(txtusplit[posu+6+i*3]), float(remove_right(txtusplit[posu+7+i*3]))])

    us=np.asarray(us)        
    us2=np.asarray(us2)            
    
    uss=(us[:,0]**2.+us[:,1]**2.)**.5
    us2s=(us2[:,0]**2.+us2[:,1]**2.)**.5
    
    
    plt.figure()
    # plt.scatter(us2s[::10],uss[::10], s=0.001)   
    axis = 0
    # plt.scatter(us2[::10,axis],us[::10,axis], s=0.001)    
    plt.scatter(us2[clip==1,axis][::10],us[clip==1,axis][::10], s=0.001)    
    # plt.scatter(us2s[clip==1][::10],uss[clip==1][::10], s=0.001)    
    plt.xlabel('RANS [m/s]')
    plt.ylabel('ML [m/s]')
    plt.title('axis:'+str(axis)+" time:"+str(itter[1])+' stat='+str(stat(us2[:,axis],us[:,axis],kind='r2')))


    X, Y, Z = density_estimation(us[:,0], us2[:,0])
    plt.figure()
    fig, ax = plt.subplots()

    # Show density
    ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
              extent=[min(us[:,0]), max(us[:,0]), min(us[:,0]), max(us[:,0])])

    # Add contour lines
    # CS = plt.contour(X, Y, Z, 10)
    # plt.clabel(CS, inline=1, fontsize=10)
    # plt.scatter(observation,model)
    plt.xlabel('observations')
    plt.ylabel('model')



class App(QWidget):

    slicomatic = False
    showfieldlist = {'U_x':0, 'U_y':1, 'U_z':2}

    def __init__(self):
        super(App, self).__init__()
        self.title = 'IIBR - OpenFoam Viewer'
        self.left = 10
        self.top = 10
        self.width = 820
        self.height = 700
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.legend = []

        button_para = QPushButton('ParaView', self)
        button_para.move(420, 20)
        button_para.clicked.connect(self.on_click_para)

        button_test = QPushButton('Error classification', self)
        button_test.move(420, 60)
        button_test.clicked.connect(self.on_click_error_classification)

        button_test = QPushButton('LSM', self)
        button_test.move(570, 60)
        button_test.clicked.connect(self.on_click_lsm)

        button_test = QPushButton('Porous', self)
        button_test.move(670, 60)
        button_test.clicked.connect(self.on_click_porous)

        button_save = QPushButton('Save', self)
        button_save.move(320, 20)
        button_save.clicked.connect(self.on_click_save_configuaration)

        button_density = QPushButton('Density', self)
        button_density.move(320, 60)
        button_density.clicked.connect(self.on_click_density)

        label3 = QLabel("Add to title:", self)
        label3.move(10, 30)
        self.addtotitle = QLineEdit(self)
        self.addtotitle.move(100, 30)
        self.addtotitle.resize(200, 20)
        self.addtotitle.setText('')
                     
        label3 = QLabel("Itteration #", self)
        label3.move(10, 350)

        self.itteration1 = QLineEdit(self)
        self.itteration1.move(100, 340)
        self.itteration1.resize(100, 40)
        self.itteration1.setText('300.0')

        self.itteration2 = QLineEdit(self)
        self.itteration2.move(200, 340)
        self.itteration2.resize(100, 40)
        self.itteration2.setText('4000.0')

        self.listfields = QListWidget(self)
        self.listfields.move(180, 420)
        self.listfields.resize(100, 140)
        self.listfields.addItem('U_x')
        self.listfields.setSelectionMode(QListWidget.MultiSelection)
        self.listfields.setCurrentRow(0)
        
        self.listfiles = QListWidget(self)
        self.listfiles.move(30, 140)
        self.listfiles.resize(280, 140)
        self.listfiles.setSelectionMode(QListWidget.MultiSelection)
        self.listfiles.setCurrentRow(0)

        button_profile = QPushButton('add file', self)
        button_profile.move(30, 110)
        button_profile.clicked.connect(self.on_click_add_file)

        button_profile = QPushButton('remove file', self)
        button_profile.move(120, 110)
        button_profile.clicked.connect(self.on_click_remove_file)        

        button_corr = QPushButton('correlation', self)
        button_corr.move(320, 150)
        button_corr.clicked.connect(self.on_click_corr)

        self.listaxis = QListWidget(self)
        self.listaxis.move(40, 420)
        self.listaxis.resize(20, 60)
        self.listaxis.addItem('x')
        self.listaxis.addItem('y')
        self.listaxis.addItem('z')
        self.listaxis.setCurrentRow(2)

        self.axispos = QLineEdit(self)
        self.axispos.move(60, 420)
        self.axispos.resize(50, 40)
        self.axispos.setText('1.0')

        button_itter = QPushButton('itterations', self)
        button_itter.move(320, 190)
        button_itter.clicked.connect(self.on_click_itter)

        label4 = QLabel("xmin, xmax", self)
        label4.move(10, 530)
        self.xmin = QLineEdit(self)
        self.xmin.move(40, 520)
        self.xmin.resize(50, 40)
        self.xmin.setText('0.0')
        self.xmax = QLineEdit(self)
        self.xmax.move(90, 520)
        self.xmax.resize(50, 40)
        self.xmax.setText('0.0')

        label5 = QLabel("ymin, ymax", self)
        label5.move(10, 570)
        self.ymin = QLineEdit(self)
        self.ymin.move(40, 560)
        self.ymin.resize(50, 40)
        self.ymin.setText('0.0')
        self.ymax = QLineEdit(self)
        self.ymax.move(90, 560)
        self.ymax.resize(50, 40)
        self.ymax.setText('0.0')

        label6 = QLabel("zmin, zmax", self)
        label6.move(10, 610)
        self.zmin = QLineEdit(self)
        self.zmin.move(40, 600)
        self.zmin.resize(50, 40)
        self.zmin.setText('0.0')
        self.zmax = QLineEdit(self)
        self.zmax.move(90, 600)
        self.zmax.resize(50, 40)
        self.zmax.setText('0.0')

        self.holdon = QCheckBox('hold on', self)
        self.holdon.move(10, 1)
        # self.holdon.toggle()
        self.holdon.stateChanged.connect(self.changeHoldOn)

        self.heightcontour = QCheckBox('height contour', self)
        self.heightcontour.move(100, 1)

        self.figsave = QCheckBox('save figure', self)
        self.figsave.move(200, 1)

        self.experiment = QComboBox(self)
        self.experiment.addItem("observationB")
        self.experiment.addItem("observationE")
        self.experiment.addItem("Michelstadt")
        self.experiment.addItem("None")
        self.experiment.move(320,320)
        button_experiment = QPushButton('corr with experiment', self)
        button_experiment.move(320, 360)
        button_experiment.clicked.connect(self.on_click_experiment)

        button_test = QPushButton('Stream Lines', self)
        button_test.move(420, 100)
        button_test.clicked.connect(self.on_click_stream_lines)

        button_test = QPushButton('eddies', self)
        button_test.move(535, 100)
        button_test.clicked.connect(self.on_click_eddies)

        button_test = QPushButton('statistics', self)
        button_test.move(420, 150)
        button_test.clicked.connect(self.on_click_statistics)

        button_learn = QPushButton('Learn', self)
        button_learn.move(420, 190)
        button_learn.clicked.connect(self.on_click_learn)

        self.learnfrom = QLineEdit(self)
        self.learnfrom.move(520, 180)
        self.learnfrom.resize(150, 40)
        self.learnfrom.setText('')      
        
        button_learn = QPushButton('sinks', self)
        button_learn.move(420, 230)
        button_learn.clicked.connect(self.on_click_sinks)

        button_learn = QPushButton('trajectory', self)
        button_learn.move(420, 270)
        button_learn.clicked.connect(self.on_click_trajectories)
       
        self.xcoord = QLineEdit(self)
        self.xcoord.move(520, 270)
        self.xcoord.resize(50, 40)
        self.xcoord.setText('0.0')

        self.ycoord = QLineEdit(self)
        self.ycoord.move(580, 270)
        self.ycoord.resize(50, 40)
        self.ycoord.setText('0.0')

        self.zcoord = QLineEdit(self)
        self.zcoord.move(640, 270)
        self.zcoord.resize(50, 40)
        self.zcoord.setText('0.0')

        self.use2 = QCheckBox('use 2 files', self)
        self.use2.move(70, 48)
                
        self.parallel = QCheckBox('parallel', self)
        self.parallel.move(170, 48)
                

        button1 = QPushButton('Plot files', self)
        button1.move(120, 320)
        button1.clicked.connect(self.on_click_plot)
        

        button = QPushButton('Diff files', self)
        button.move(20, 320)
        button.clicked.connect(self.on_click_diff_files)

        button_update_fields = QPushButton('Update fields', self)
        button_update_fields.move(180, 390)
        button_update_fields.clicked.connect(self.on_click_update_fields)
        
        button_profile = QPushButton('profile', self)
        button_profile.move(320, 120)
        button_profile.clicked.connect(self.on_click_profile)

        button_profile = QPushButton('profile2', self)
        button_profile.move(420, 120)
        button_profile.clicked.connect(self.on_click_profile2)

        button_profile = QPushButton('profile3', self)
        button_profile.move(520, 120)
        button_profile.clicked.connect(self.on_click_profile3)

        button2 = QPushButton('Diff times', self)
        button2.move(220, 320)
        button2.clicked.connect(self.on_click_diff_times)
        
        
        # https: // github.com / mmisono / pyqt5 - example / blob / master / plot.py
        # fig = plt.figure(figsize=(4, 3), dpi=100)
        # self.axes = fig.add_subplot(111)
        # self.axes.hold(False)
        # self.axes.move(590, 400)
        # self.axes.resize(250, 240)

        # Standard Matplotlib code to generate the plot

        # initialize the canvas where the Figure renders into
        # FigureCanvas.__init__(self,fig)
        #
        #
        # t = np.arange(0.0, 3.0, 0.01)
        # s = np.sin(2 * np.pi * t)
        # self.axes.plot(t, s)
        # self.canvas.draw()

        try:
            with open('ofplot.json') as data_file:
                data = json.load(data_file)
            self.itteration1.setText(data['itteration1'])
            self.itteration2.setText(data['itteration2'])
#            self.filename1.setText(data['file1'])
#            self.filename2.setText(data['file2'])
            self.axispos.setText(data['axisValue'])
            self.listaxis.setCurrentRow(data['axisIndex'])
            self.xmin.setText(data['clipxmin'])
            self.xmax.setText(data['clipxmax'])
            self.ymin.setText(data['clipymin'])
            self.ymax.setText(data['clipymax'])
            self.zmin.setText(data['clipzmin'])
            self.zmax.setText(data['clipzmax'])
            self.xcoord.setText(data['xcoord'])
            self.ycoord.setText(data['ycoord'])
            self.zcoord.setText(data['zcoord'])
            self.experiment.setCurrentText(data['experiment'])
            self.heightcontour.setCheckState(data['heightcontour'])
            self.use2.setCheckState(data['use2'])
            self.learnfrom.setText(data['learnfrom'])
            self.listfields.clear()
            for fieldindex in range(len(data['fields'])):
                self.listfields.addItem(data['fields'][fieldindex])
            self.listfields.setCurrentRow(data['fieldIndex'])
            for fieldindex in range(len(data['files'])):
                self.listfiles.addItem(data['files'][fieldindex])            
        except:
            print ('no valid configuration file')
        self.show()

#        try:
#            self.reader1 = bse.ReadCase('casename', data['files'][0], CaseType='Reconstructed Case')  # 'Reconstructed Case')
#            print('Initial timelist print:', self.reader1.TimestepValues)
#        except:
#            pass

    @pyqtSlot()

    def on_click_trajectories(self):
        print('start traj', datetime.datetime.now())
                        
 #self.filename1.text()
        if 5==7:
            textboxfile1value = r'/ibdata2/nirb/openFOAM/vortex/tlvs2/tlvs2.foam'
            textboxfile1value = r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomega/windAroundCube.foam'

            textboxtime1value = 999999999
            db = get_clip(textboxfile1value, [textboxtime1value],
                      clipxmin=float(0), clipxmax=float(0),
                      clipymin=float(0), clipymax=float(0),
                      clipzmin=float(0), clipzmax=float(0))
            xmin=0
            xmax=0
            ymin=0
            ymax=0
            zmin=0
            zmax=0
            r=0
        else:
            files=[]        
            listItems=self.listfiles.selectedItems()
            if listItems:
                for item in listItems:
                   files.append(item.text())
            else:
                for item in range(self.listfiles.count()):
                    files.append(self.listfiles.item(item).text())
    
            for filesindex in range(len(files)):
                textboxfile1value = files[filesindex] #self.filename1.text()
            
            textboxtime1value = float(self.itteration1.text())
            sel = self.listfields.selectedIndexes()
            print('1start traj', datetime.datetime.now())
            
            db = get_clip(textboxfile1value, [textboxtime1value],
                      clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                      clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                      clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()))
            xmin = float(self.xmin.text())
            xmax = float(self.xmax.text())
            ymin = float(self.ymin.text())
            ymax = float(self.ymax.text())
            zmin = float(self.zmin.text())
            zmax = float(self.zmax.text())

        # ap3 = db[['x', 'y', 'z', sel[0].data()]]  # at the ground (1 meter)
        ap3 = db[0][['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)
        points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
        # values = ap3[sel[0].data()]
        valuesux = ap3['U_x']
        valuesuy = ap3['U_y']
        valuesuz = ap3['U_z']

        if xmin==xmax:
            xmin=ap3['x'].min()
            xmax=ap3['x'].max()
        if ymin==ymax:
            ymin=ap3['y'].min()
            ymax=ap3['y'].max()
        if zmin==zmax:
            zmin=ap3['z'].min()
            zmax=ap3['z'].max()
        # ticks = 600j
        grid_points = 42000000  # add 00
        meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)

        resolution = (meters / grid_points)**(1./3)
        if resolution==0:
            print('no clip boundaries in one of the axis')
        # resolution = 0.1
        ticksx = int((xmax-xmin) / resolution) * 1j
        ticksy = int((ymax-ymin) / resolution) * 1j
        ticksz = int((zmax-zmin) / resolution) * 1j
        print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)
        print('3start traj', datetime.datetime.now())

        xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
        gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        griduy = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuy), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        griduz = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuz), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings

        x = np.linspace(ap3['x'].min(), ap3['x'].max(), int(ticksx.imag))
        y = np.linspace(ap3['y'].min(), ap3['y'].max(), int(ticksy.imag))
        z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticksz.imag))

        print('mid traj1', datetime.datetime.now())
        
        if self.use2.checkState():
            textboxfile2value = r'/ibdata2/nirb/openFOAM/vortex/tlvs2z20/tlvs2.foam'
 #self.filename2.text()
    
            db2 = get_clip(textboxfile2value, [textboxtime1value],
                      clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                      clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                      clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()))
    
    
            ap32 = db2[0][['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)
            points32 = ap32[['x', 'y', 'z']].values  # numpy[cells,3]
            # values = ap3[sel[0].data()]
            valuesux2 = ap32['U_x']
            valuesuy2 = ap32['U_y']
            valuesuz2 = ap32['U_z']
    
            gridux2 = griddata((points32[:, 0], points32[:, 1], points32[:, 2]), np.asarray(valuesux2), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            griduy2 = griddata((points32[:, 0], points32[:, 1], points32[:, 2]), np.asarray(valuesuy2), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            griduz2 = griddata((points32[:, 0], points32[:, 1], points32[:, 2]), np.asarray(valuesuz2), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
       
            print('mid traj2', datetime.datetime.now())


        dt = 0.1
        duration = 1000000
        dx = x[1]-x[0]
        dy = y[1]-y[0]
        dz = z[1]-z[0]
        realizations = 2
        sumerror = 0.
        sumdist = 0.
        for r in range(realizations):
            if r%10000==0:
                print('r=',r)
            xcoord = x[int(len(x)/2)]#float(self.xcoord.text())
            ycoord = y[int(len(y)/2)]+r*(y[1]-y[0])#float(self.ycoord.text())+r*(y[1]-y[0])
            zcoord = z[5]#float(self.zcoord.text())
    
            trajx, trajy, trajz, dur = traj(xcoord, ycoord, zcoord, x, y, z, gridux, griduy, griduz, duration=duration, dt=dt)
#            speed = gridux[trajx[-1],trajy[-1],trajz[-1]]
#            print(dur, 'traj-',trajx[0],'<->',trajx[-1],'|||',trajy[0],'<->',trajy[-1],'|||',trajz[0],'<->',trajz[-1], gridux[trajx[-1],trajy[-1],trajz[-1]])
#            if self.use2.checkState():
#                trajx2, trajy2, trajz2, dur = traj(xcoord, ycoord, zcoord, x, y, z, gridux2, griduy2, griduz2, duration=duration, dt=dt)
#                thiserror = ((dx*(trajx[-1]-trajx2[-1]))**2.+(dy*(trajy[-1]-trajy2[-1]))**2.+(dz*(trajz[-1]-trajz2[-1]))**2.)**.5
#                thisdist = ((dx*(trajx2[0]-trajx2[-1]))**2.+(dy*(trajy2[0]-trajy2[-1]))**2.+(dz*(trajz2[0]-trajz2[-1]))**2.)**.5
#                if thiserror>thisdist*0.4:
#                    trajx, trajy, trajz, dur = traj(xcoord, ycoord, zcoord, x, y, z, gridux, griduy, griduz, duration=duration, dt=dt, verbose=False)
#                    trajx2, trajy2, trajz2 = traj(xcoord, ycoord, zcoord, x, y, z, gridux2, griduy2, griduz2, duration=duration, dt=dt, verbose=False)
#                if trajx[0]!=trajx[-1]:
##                    print(dur, '>>>',trajx[0],'<->',trajx[-1],trajx2[-1],'|||',trajy[0],'<->',trajy[-1],trajy2[-1],'|||',trajz[0],'<->',trajz[-1],trajz2[-1])
##                print('traj-',trajx[0],'<->',trajx[-1],'|||',trajy[0],'<->',trajy[-1],'|||',trajz[0],'<->',trajz[-1])
##                print('traj2-',trajx2[0],'<->',trajx2[-1],'|||',trajy2[0],'<->',trajy2[-1],'|||',trajz2[0],'<->',trajz2[-1])
#                    print(r, 'diff:*************',thiserror,thisdist)
#                sumerror += thiserror
#                sumdist += thisdist
##                speed2 = gridux[trajx2[-1],trajy2[-1],trajz2[-1]]                  
                
#        print('sum error, dist, mean1, mean2', sumerror, sumdist, sumerror/sumdist, sumerror/realizations)                
        plt.figure()
        plt.imshow(gridux[:,:,int(trajz[0])], origin='lower')  # jet, Paired
        plt.colorbar()
        plt.plot(trajx, trajy, marker="+")
        plt.xlabel('x')
        plt.ylabel('y')
#        plt.plot(trajx2, trajy2, marker="*")
        plt.show()
        plt.figure()
        plt.plot(trajx, trajz, marker="+")
        plt.xlabel('x')
        plt.ylabel('z')


#    def filename1changed(self):
#        self.reader1 = None

    def on_click_stream_lines(self):
        cc=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
#            cc=['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        files=[]        

#        files.append(r'/ibdata2/nirb/openFOAM/angle/fvoptions/cylinder.foam')
#        files.append(r'/ibdata2/nirb/openFOAM/angle/fvoptionszero/cylinder.foam')        
#        files.append(r'/ibdata2/nirb/openFOAM/angle/kleombenchmark2/cylinder.foam')   
        files.append(r'/ibdata2/nirb/openFOAM/vortex/tlvs2/tlvs2.foam')
        files.append(r'/ibdata2/nirb/openFOAM/vortex/tlvs2z40/tlvs2.foam')
        files.append(r'/ibdata2/nirb/openFOAM/vortex/tlvs2z20/tlvs2.foam')
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
            for item in listItems:
               files.append(item.text())
        else:
            for item in range(self.listfiles.count()):
                files.append(self.listfiles.item(item).text())
       
        textboxtime1value = 100000.
        xcoord = float (155000) # float(218750.0)
        ycoord = float (564000)# float(631250.0)
        zcoord = float (1000) # float(909)      
        
        xcoord = 178950
        ycoord = 664200
        zcoord = 20


        textboxtime1value = float(self.itteration1.text())
        xcoord = float(self.xcoord.text())
        ycoord = float(self.ycoord.text())
        zcoord = float(self.zcoord.text())

        plt.figure()
        timelength=[]
        f=0
        onest=[]
        fsti=[]
        for f in range(len(files)):    
            print(files[f])
            reader = bse.ReadCase(files[f], files[f], CaseType='Decomposed Case')  # ' Reconstructed Decomposed Case')
            lines = 110+1 # 2+1
                
            streamTracer1 = pvsimple.StreamTracer(Input=reader,SeedType='Point Source')
            streamTracer1 = pvsimple.StreamTracer(Input=reader,SeedType='High Resolution Line Source')
            streamTracer1.Vectors = ['POINTS', 'U']
            streamTracer1.MaximumStreamlineLength = 300099.0  # initial length 850
            streamTracer1.InitialStepLength = 0.002
            streamTracer1.MaximumStepLength = 0.1
            streamTracer1.MinimumStepLength = 0.001
            streamTracer1.MaximumSteps = 99999999
    
            streamTracer1.SeedType.Point1 = [xcoord-1200, ycoord-1000, zcoord]#[220750.0, 631250.0, 909.0]
            streamTracer1.SeedType.Point2 = [xcoord-1200, ycoord+1000, zcoord]#[220750.0, 633250.0, 909.0]
            streamTracer1.SeedType.Resolution = lines
            streamTracer1.IntegrationDirection = 'FORWARD'
            streamTracer1.UpdatePipeline(textboxtime1value)
            itr = bse.to_pandas(reader, streamTracer1, timelist=textboxtime1value)
            dbb = itr.next()
            db = dbb[dbb.keys()[0]]
            if len(db) == 0:
                print('There are no records in stream lines, please check coordinates limits or slice location')
    #        print('AAC',db.iloc[0]['x']+db.iloc[0]['U_x']*db.iloc[1]['IntegrationTime'],db.iloc[1]['x'])
    
            sti = db.index[db['IntegrationTime']==0.].tolist()
            print('sti',len(sti), sti)
            fsti.append(len(sti)-1)
            dist=[]
            dist0=[]
            i=0
            for i in range(len(sti)-1):                
                onest.append(db.iloc[sti[i]:sti[i+1]])
                plt.subplot(2,2,1)
                plt.plot(onest[i]['x'],onest[i]['y'],'-+', color=cc[f])  #cc[f]        
                plt.text(onest[i]['x'].iloc[-1],onest[i]['z'].iloc[-1],int(db.iloc[sti[i+1]-1]['IntegrationTime']))
                plt.title('xy')
                plt.subplot(2,2,2)
                plt.plot(onest[i]['x'],onest[i]['z'],'-+', color=cc[f])  #cc[f]        
                plt.text(onest[i]['x'].iloc[-1],onest[i]['z'].iloc[-1],int(db.iloc[sti[i+1]-1]['IntegrationTime']))
                plt.title('xz')
                plt.subplot(2,2,3)
                plt.plot(onest[i]['x'],onest[i]['U_x']**2.+onest[i]['U_y']**2.+onest[i]['U_z']**2.,'-+', color=cc[f])  #cc[f]        
                plt.text(onest[i]['x'].iloc[-1],onest[i]['z'].iloc[-1],int(db.iloc[sti[i+1]-1]['IntegrationTime']))
                plt.title('xu')
#                print(i,onest[i]['x'])
                timelength.append(db.iloc[sti[i+1]-1]['IntegrationTime'])
                print('it',i,timelength[-1])
#            print ('@@@@@@', onest[0][['x','y','U_x','IntegrationTime']])
#            print ('!!!!!!', onest2[['x','y','U_x','IntegrationTime']])
        plt.axis('equal')    
        lgnd=[]
        for i in range(len(files)):
            lgnd.append(shrinktitle(files[i]))
            if len(files)>1:
                dist = np.asarray(dist)
                dist0 = np.asarray(dist0)
#                print ('dist diff mean+-std out of dist0 mean+-std, speed0, speed1 (last)')
#                print (round(dist.mean(),2), round(dist.std(),2), round(dist0.mean(),2), round(dist0.std(),2), round(dist.mean()/dist0.mean(),2)*100.,'%', onest[0].iloc[0]['U_x'], onest[0].iloc[0]['U_x'])
        plt.legend(lgnd)
        ax = plt.gca()       
        leg = ax.get_legend()
        for i in range(len(files)):
            leg.legendHandles[i].set_color(cc[i])
        print('TimeLength', timelength)
        
        tr=[]
        for l in range(lines-1):
            mint=9999999999999999999999999999
            for f in range(len(files)):
                currentt=onest[f*(lines-1)+l]['IntegrationTime'].max()
                mint=min(mint,currentt)
                print(mint, currentt)
            tr.append(mint)
        
        for l in range(lines-1):
            print('x')
            for f in range(len(files)):           
                startx=onest[f*(lines-1)+l]['x'][onest[f*(lines-1)+l]['IntegrationTime']==0].max()
                starty=onest[f*(lines-1)+l]['y'][onest[f*(lines-1)+l]['IntegrationTime']==0].max()
                startz=onest[f*(lines-1)+l]['z'][onest[f*(lines-1)+l]['IntegrationTime']==0].max()
                tempt=onest[f*(lines-1)+l]['IntegrationTime'][onest[f*(lines-1)+l]['IntegrationTime']<tr[l]].max()
                tempx=onest[f*(lines-1)+l]['x'][onest[f*(lines-1)+l]['IntegrationTime']==tempt].max()
                tempy=onest[f*(lines-1)+l]['y'][onest[f*(lines-1)+l]['IntegrationTime']==tempt].max()
                tempz=onest[f*(lines-1)+l]['z'][onest[f*(lines-1)+l]['IntegrationTime']==tempt].max()
                print(tempt,startx,tempx, starty, tempy, startz, tempz)

        plt.figure()
        fsti=np.asarray(fsti)
        for f in range(1,len(files)):           
            tx=[]
            ty=[]
            sfsti=fsti[:f].sum()
            for l in range(fsti[f]):
                tempt=onest[sfsti+l]['IntegrationTime'].max()
                tx.append(onest[sfsti+l]['x'][onest[sfsti+l]['IntegrationTime']==tempt].max())
                ty.append(onest[sfsti+l]['y'][onest[sfsti+l]['IntegrationTime']==tempt].max())
            plt.plot(tx,ty,'.'+cc[f])                
        
    def on_click_eddies(self):
        itter=[float(self.itteration1.text())]
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
            for item in listItems:
               files.append(item.text())
        else:
            for item in range(self.listfiles.count()):
                files.append(self.listfiles.item(item).text())
        axisindex = int(self.listaxis.currentRow())
        axispos = float(self.axispos.text())
        xmin=float(self.xmin.text())
        xmax=float(self.xmax.text())
        ymin=float(self.ymin.text())
        ymax=float(self.ymax.text())
        zmin=float(self.zmin.text())
        zmax=float(self.zmax.text())
        currentRow=self.listaxis.currentRow()
            
        eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow)
           
    def on_click_eddies2(self):
        textboxtime1value = float(self.itteration1.text())
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
            for item in listItems:
               files.append(item.text())
        else:
            for item in range(self.listfiles.count()):
                files.append(self.listfiles.item(item).text())

        for filesindex in range(len(files)):
            textboxfile1value = files[filesindex] #self.filename1.text()
            axisindex = int(self.listaxis.currentRow())
            axispos = float(self.axispos.text())
            db = get_slice(textboxfile1value, textboxtime1value, axisindex, axispos,
                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()))

            cma = plt.cm.OrRd
            cma.set_bad(alpha=0.0)

            ticks = 400
            xmin, xmax, ymin, ymax, u = makegrid(db, 'U_x', self.listaxis.currentRow(), ticks=ticks)
            xmin, xmax, ymin, ymax, w = makegrid(db, 'U_z', self.listaxis.currentRow(), ticks=ticks)
            u[u==0]= np.NaN
            w[w==0]= np.NaN
            xi = np.linspace(xmin,xmax,ticks)
            zi = np.linspace(ymin,ymax,ticks)
            xi = np.linspace(0,ticks,ticks)
            zi = np.linspace(0,ticks,ticks)
            dx = (xmax-xmin)/ticks
            dz = (ymax-ymin)/ticks
            xx,zz = np.meshgrid(xi,zi)
            
           
            u1 = np.sign(u)
            w1 = np.sign(w)
            
            ud = np.roll(u1, -1, axis=0)
            uu = np.roll(u1, 1, axis=0)
            ur = np.roll(u1, -1, axis=1)
            ul = np.roll(u1, 1, axis=1)
            wd = np.roll(w1, -1, axis=0)
            wu = np.roll(w1, 1, axis=0)
            wr = np.roll(w1, -1, axis=1)
            wl = np.roll(w1, 1, axis=1)
            
            test0 = np.zeros_like(u)
            test0[(u1==1)&(w1==1)]=1
            test0[(u1==1)&(w1==-1)]=2
            test0[(u1==-1)&(w1==1)]=3
            test0[(u1==-1)&(w1==-1)]=4
            test1=np.roll(test0,1,axis=0) #come from above
            test2=np.roll(test0,1,axis=1) #come from left
            test3=np.roll(test2,1,axis=0)
            
            test = np.zeros_like(u)+1
            # 2 directions are not good
            test[(test0==test1) & (test0==test2)] = 0
            test[(test0==test1) & (test0==test3)] = 0
            test[(test0==test3) & (test0==test2)] = 0 #  | np.isnan(test0) | np.isnan(test1) | np.isnan(test2) | np.isnan(test3)
            test[(test1==test3) & (test1==test2)] = 0
            # 2 direction parallel
            test[(test2==test3) & (test0==test1)] = 0 #for vertical border
            test[(test1==test3) & (test0==test2)] = 0 #for horizontal border
            
            test[np.isnan(u1)]=np.NaN
            test[test0==0]=np.NaN
            test[test1==0]=np.NaN
            test[test2==0]=np.NaN
            test[test3==0]=np.NaN
            test[:,0]=np.NaN
            test[:,-1]=np.NaN
            test[0,:]=np.NaN
            test[-1,:]=np.NaN
                 
            plt.figure()
            plt.imshow(test0,origin='lower')
            plt.colorbar()
            
            plt.figure()                      
            plt.streamplot(xx,zz,u, w, density=15, color=u, linewidth=2, cmap=plt.cm.autumn)
            plt.colorbar()    
            plt.title('stramlines')

            # make some data
            y=zi#numpy.linspace(0,10,40)
            x=xi#numpy.linspace(0,10,50)
            v=w
            v[np.isnan(v)]=0
            u[np.isnan(u)]=0
            X,Y=np.meshgrid(x,y)
            print('shapes', X.shape, v.shape)
            print(v)

            # integrate to make an intial guess
            intx=integrate.cumtrapz(v,X,axis=1,initial=0)[0]
            inty=integrate.cumtrapz(u,Y,axis=0,initial=0)
            psi1=intx-inty
#            plt.figure()
#            plt.imshow(inty,origin='lower')
#            plt.title('vvxxxv')
            print('shae intx',intx.shape)
            print('shae inty',inty.shape)
#            plt.figure()
#            plt.imshow(inty,origin='lower')
#            plt.title('vvyyyv')
            
            intx=integrate.cumtrapz(v,X,axis=1,initial=0)
            inty=integrate.cumtrapz(u,Y,axis=0,initial=0)[:,0][:,None]
            psi2=intx-inty
            print('shae2 intx',intx.shape)
            print('shae2 inty',inty.shape)
            
            psi=0.5*(psi1+psi2)
            
            intx=integrate.cumtrapz(u,X,axis=1,initial=0)[0]
            inty=integrate.cumtrapz(v,Y,axis=0,initial=0)
            chi1=intx+inty
            
            intx=integrate.cumtrapz(u,X,axis=1,initial=0)
            inty=integrate.cumtrapz(v,Y,axis=0,initial=0)[:,0][:,None]
            chi2=intx+inty
            
            chi=0.5*(chi1+chi2)
            print('shapes2', psi.shape, chi.shape)
            
            # pad to add 1 row/col
            psi[psi<-30]=-30
            plt.figure()
            plt.imshow(psi,origin='lower')
            plt.colorbar()
            plt.title('psi')
            
            psi1[psi1<-30]=-30
            plt.figure()
            plt.imshow(psi1,origin='lower')
            plt.colorbar()
            plt.title('psi1')
            
            psi2[psi2<-30]=-30
            plt.figure()
            plt.imshow(psi2,origin='lower')
            plt.colorbar()
            plt.title('psi2')

#            plt.figure()
#            plt.imshow(psi1*psi1+psi2*psi2,origin='lower')
#            plt.colorbar()
#            plt.title('psi1122')
            
#            plt.figure()
#            plt.imshow(chi,origin='lower')
#            plt.colorbar()
#            plt.title('chi')
            
            
            
            plt.figure()
            #plt.imshow(w1, origin='lower')  # jet, Paired
#            s = np.power(np.add(np.power(u1,2), np.power(w1,2)),0.5) #could do (u**2+v**2)**0.5, other ways...
            if (self.learnfrom.text()==u''):
                scl = 50
            else:
                scl = float(self.learnfrom.text()) #100
            plt.quiver(xx,zz,u,w, scale=scl)
            
            
#            u2 = np.roll(u1, 1, axis=0)
#            u0 = np.roll(u1, 1, axis=1)
#            w2 = np.roll(w1, 1, axis=0)
#            w0 = np.roll(w1, 1, axis=1)
#                      
#            s = np.zeros_like(u1)
#            s[(uu==1) & (ud==-1) & (ul==1) & (ur==-1)
#              & (wu==-1) & (wd==1) & (wl==1) & (wr==-1)] = 1
#            s[(uu==-1) & (ud==1) & (ul==1) & (ur==-1) 
#              & (wu==-1) & (wd==1) & (wl==-1) & (wr==1)] = 2

            #s[(u1!=u2) & (w1!=w2)] = 1
#            s[((((u1==-1) & (u2==1)) | (u1==1) & (u2==-1)) & (((w1==-1)&(w0==1)) | ((w1==1)&(w0==-1))))] = 1

#            s[(((u1==-1) & (u0==-1) & (u2== 1) & (w1==-1) & (w0== 1) & (w2==1)) | 
#               ((u1==1)  & (u0==1)  & (u2==-1) & (w1==-1) & (w0==1) & (w2==-1)) | # by experiment !!!
#               ((u1== 1) & (u0== 1) & (u2==-1) & (w1== 1) & (w0==-1) & (w2==-1)))] = 1
               
#            print(u1[s==1], u0[s==1], u2[s==1], w1[s==1], w0[s==1], w2[s==-1])
            s = test.copy()
            mask = np.zeros_like(test)
            s = np.zeros_like(test)
            s[np.isnan(test)]=np.NaN
            centers=np.where(test==1)
            centerslen = len(centers[0])
            print('centers', centers, centerslen)
            point = 6 # please don't use -1 because it's negative and we assume positive
#            centers[0][point]=3
#            centers[1][point]=180
            print('detect veortex #',point)
            exitreason='loop'
            radius=0
            radiusjump = 0.2
#            while exitreason=='loop':
#                sbackup = s.copy()
#                s = np.zeros_like(test)
#                s[np.isnan(test)]=np.NaN
#                if (radius>88888):
#                    exitreason='max'
#                    print('exit reason=',exitreason)
#                    break
#                    print('exit reason2=',exitreason)
#                radius += radiusjump
#                t1=[centers[0][point],centers[1][point]]
#                exitreason='end'
#                t1[1] += radius*(1) # we will start contour one grid above
#                t1[0] += radius*(1) # we will start contour one grid above
#                print('radius: ', radius, exitreason,t1)
#                counter=2
#                s[int(t1[0]),(t1[1])] = -counter  # assign first cell as a cell that was visited
#                for step in range(94000000):
#                    dt = min(abs(dx/u[t1[0],t1[1]]),dz/abs(w[t1[0],t1[1]])) / (-40.) # negative time because we go backwords
#                    velu, velw = interpvel(u,w,t1[0],t1[1])
#                    stepx=dt*u[int(t1[0]),int(t1[1])]/dx
#                    stepz=dt*w[int(t1[0]),int(t1[1])]/dz
#                    stepx=dt*velu/dx
#                    stepz=dt*velw/dz
##                    print('first t1',step, dt,t1, stepx,stepz,u[int(t1[0]),int(t1[1])], w[int(t1[0]),int(t1[1])])
#                    prevlocation=[t1[0],t1[1]]
#                    t1[1]+=stepx
#                    t1[0]+=stepz
#                    if s[int(t1[0]),int(t1[1])]==0:  # we haven't been here before
##                        print(step, 't debug',dt,t1,prevlocation,s[t1[0],t1[1]], stepx , stepz)
#                        s[int(t1[0]),int(t1[1])] = -counter
#                        counter+=1
#                    else:
#                        if np.isnan(s[int(t1[0]),int(t1[1])]):
#                            exitreason='boundary'
#                            print ('nan in ',t1[0],t1[1])
#                            break
#                        if (((int(t1[0])!=int(prevlocation[0])) | (int(t1[1])!=int(prevlocation[1]))) & (s[int(t1[0]),int(t1[1])]!=centerslen + radius) & (abs(s[int(t1[0]),int(t1[1])]+counter)>2)):
#                            exitreason='loop'
#                            print('loop after ', counter, t1[0],t1[1],prevlocation,s[int(t1[0]),int(t1[1])],centerslen + radius)
#                            break
##                        print(step, counter, 't2 debug',dt,t1,prevlocation,s[t1[0],t1[1]], stepx , stepz, int(t1[0]),int(prevlocation[0]),int(t1[1]),int(prevlocation[1]))
#            if exitreason=='end':
#                #s[s<0]=point*100+radius
#                s = sbackup.copy()
#                print('end loop', step, counter, point, radius,centerslen,radiusjump)
#            elif exitreason=='boundary':
#                print('boundary', step, counter, point, radius,centerslen,radiusjump)
#                s = sbackup.copy()
#                s[s<0]=radius + centerslen -radiusjump               
#            elif exitreason=='max':
#                print('max', step, counter, point, radius,centerslen,radiusjump,t1)
#                s = sbackup.copy()
#                s[s<0]=radius + centerslen -radiusjump
#            else: # loop
#                s[s<0]=radius + centerslen
#                print(exitreason)
#                print('finish loop after:',step ,' steps',radius ,centerslen)
#            if radius>radiusjump:
#                s[int(centers[0][point]),int(centers[1][point])] = 0 # remove center point
#                paintinside(s,[int(centers[0][point]),int(centers[1][point])],point, radius + centerslen-radiusjump)
#                s[s==(radius + centerslen-radiusjump)]=0
#            else:
#                print ('no vortex at this point')
#            
#            s[np.isnan(u1)]=np.NaN
#            plt.imshow(s, origin='lower')
#            plt.plot(centers[1][:],centers[0][:],'*')
#            print('vortex size is ', np.sum(s==point),' grid points')
##            print('U (min, mean, max):', u[s==point].min(),  u[s==point].mean(), u[s==point].max())
##            print('W (min, mean, max):', w[s==point].min(),  w[s==point].mean(), w[s==point].max())            
##            print('U (sum,abs, frac):', u[s==point].sum(),  np.abs(u[s==point]).sum(),u[s==point].sum()/np.abs(u[s==point]).sum())
##            print('W (sum,abs, frac):', w[s==point].sum(),  np.abs(w[s==point]).sum(),w[s==point].sum()/np.abs(w[s==point]).sum())
#            plt.colorbar()
#            plt.show()
#            print(centers)
            
            mask = np.zeros_like(test)
#            mask = np.zeros([15,15])
#            point=0
#            centers=np.asarray([[8],[7]])

            point = 0 
            centers=np.asarray([[8],[337]])
            
            mask[int(centers[0][point]),int(centers[1][point])]=5
            print(centers)
            for rad in range(0,1):
                # build the trajectories
                print('                    rad=',rad)
                for north in range(-rad,rad+1):
                    print('north=',north)
                    if not(np.isnan(u[int(centers[0][point])+rad,int(centers[1][point])+north])):
                        mask[int(centers[0][point])+rad,int(centers[1][point])+north] = self.trajectory(u, w, dx,dz, [int(centers[0][point])+rad,int(centers[1][point])+north], mask)
                    if not(np.isnan(u[int(centers[0][point])-rad,int(centers[1][point])+north])):  
                        mask[int(centers[0][point])-rad,int(centers[1][point])+north] = self.trajectory(u, w, dx,dz, [int(centers[0][point])-rad,int(centers[1][point])+north], mask)
                for east in range(-rad+1,rad):
                    print('east=',east)                    
                    if not(np.isnan(u[int(centers[0][point])+east,int(centers[1][point])+rad])):
                        mask[int(centers[0][point])+east,int(centers[1][point])+rad] =  self.trajectory(u, w, dx,dz, [int(centers[0][point])+east,int(centers[1][point])+rad], mask)
                    if not(np.isnan(u[int(centers[0][point])+east,int(centers[1][point])-rad])):
                        mask[int(centers[0][point])+east,int(centers[1][point])-rad] =  self.trajectory(u, w, dx,dz, [int(centers[0][point])+east,int(centers[1][point])-rad], mask)
                # check if it is the last radius
                last = False
                for north in range(-rad,rad+1):
                    if mask[int(centers[0][point])+rad,int(centers[1][point])+north]: last = True
                    if mask[int(centers[0][point])-rad,int(centers[1][point])+north]: last = True
                for east in range(-rad+1,rad):
                    if mask[int(centers[0][point])+east,int(centers[1][point])+rad]: last = True
                    if mask[int(centers[0][point])+east,int(centers[1][point])-rad]: last = True
                         
            plt.figure()
            plt.imshow(mask, origin='lower')
            plt.title('mask')
            plt.colorbar()
            plt.show()         
            
            
    def trajectory(self,u,w,dx,dz,location, mask):
        s = np.zeros_like(u)
        s[np.isnan(u)]=np.NaN
        t1=[location[0], location[1]]
        exitreason='end'
        counter=2
        s[int(t1[0]),(t1[1])] = -counter  # assign first cell as a cell that was visited
#        print ('debug1')
        distance=[]
        side2=[]
        xx=[0,0]
        yy=[0,0]
#        print ('debug1')
        for step in range(940000):
#            print('t1',t1,step)
            dt = min(abs(dx/u[int(t1[0]),int(t1[1])]),dz/abs(w[int(t1[0]),int(t1[1])])) / (-40.) # negative time because we go backwords
            velu, velw = interpvel(u,w,t1[0],t1[1])
#            print('t1',t1,step,velu,velw)
            if np.isnan(velu) | np.isnan(velw):
                exitreason='boundary'
#                    print ('nan in ',t1[0],t1[1])
                break
            
            stepx=dt*velu/dx
            stepz=dt*velw/dz
#                    print('first t1',step, dt,t1, stepx,stepz,u[int(t1[0]),int(t1[1])], w[int(t1[0]),int(t1[1])])

            prevlocation=[t1[0],t1[1]]
            t1[1]+=stepx
            t1[0]+=stepz
#            print('t3',t1,step,velu,velw)
            distance.append((t1[0]-location[0])**2.+(t1[1]-location[1])**2.)
            yy.append(t1[0])
            xx.append(t1[1])
            ax=xx[-2]
            ay=yy[-2]
            bx=xx[-3]
            by=yy[-3]
            cx=xx[-1]
            cy=yy[-1]
            dxl=xx[-2]-xx[-3] # dx of the vector of the previous line
            dxp=xx[-1]-xx[-2] # dx of the vector for the new point
            dyl=yy[-2]-yy[-3] # dx of the vector of the previous line
            dyp=yy[-1]-yy[-2] # dx of the vector for the new 
            if ((dxl>=0) and (dyl>=0)):
                quarterl=1
            elif ((dxl>=0) and (dyl<0)):
                quarterl=2
            elif ((dxl>0) and (dyl>=0)):
                quarterl=4
            else:
                quarterl=3
            if ((dxp>=0) and (dyp>=0)):
                quarterp=1
            elif ((dxp>=0) and (dyp<0)):
                quarterp=2
            elif ((dxp>0) and (dyp>=0)):
                quarterp=4
            else:
                quarterp=3
            left=0
            right=1
            problem=2
            side=problem
            
            if ((quarterl==1) and (quarterp==4)):
                side=left
            if ((quarterl==4) and (quarterp==1)):
                side=right
            if ((quarterl==2) and (quarterp==1)):
                side=left
            if ((quarterl==1) and (quarterp==2)):
                side=right
            if ((quarterl==3) and (quarterp==2)):
                side=left
            if ((quarterl==2) and (quarterp==3)):
                side=right
            if ((quarterl==4) and (quarterp==3)):
                side=left
            if ((quarterl==3) and (quarterp==4)):
                side=right
            
            if dxp==0:
                angp=0.
            else:
                angp=dyp/dxp
            if dxl==0:
                angl=0.
            else:                   
                angl=dyl/dxl

            if ((quarterl==quarterp) and (quarterl==1 or quarterl==4)):
                if angp>angl:
                    side=right
                else:
                    side=left
                    
            if ((quarterl==quarterp) and (quarterl==2 or quarterl==3)):
                if angp<angl:
                    side=right
                else:
                    side=left

            if ((quarterl==1 and quarterp==3) or (quarterl==3 and quarterp==1)):
                if angp>angl:
                    side=right
                else:
                    side=left                

            if ((quarterl==2 and quarterp==4) or (quarterl==4 and quarterp==2)):
                if angp<angl:
                    side=right
                else:
                    side=left                

                    
#            if side==2:
#                print ('side21---', dxp,dyp, dxl,dyl,len(xx))
#                
#                print ('side2---', dxp,dyp, dxl,dyl,quarterp,quarterl,angp,angl, len(xx))
            side2.append(side)
#            side2.append(((bx - ax)*(cy - by) - (by - ay)*(cx - ax)))
#            if len(xx)%500==0:
#                print(xx[-3],xx[-2],xx[-1],yy[-3],yy[-2],yy[-1], side2[-1], len(xx))            
#                print(bx-ax,cy-by,by-ay,cx-ax, side2[-1])
#                print ('side1---', dxp,dyp, dxl,dyl,quarterp,quarterl,angp,angl, len(xx))
            if ((t1[0]<1) | (t1[1]<1) | (t1[0]>s.shape[0]-1) | (t1[1]>s.shape[1]-1)):
#                    print('t22',t1,step,velu,velw)

                    exitreason='boundary'
#                    print ('nan in ',t1[0],t1[1])
                    break

            if np.isnan(s[int(t1[0]),int(t1[1])]):
#                    print('t23',t1,step,velu,velw)

                    exitreason='boundary'
#                    print ('nan in ',t1[0],t1[1])
                    break
                
                
#            print('t2',t1,step,velu,velw)
#            if mask[int(t1[0]),int(t1[1])]>0:  # we haven't been here before
#                    exitreason='mask'
#                    break

            if s[int(t1[0]),int(t1[1])]==0:  # we haven't been here before
#                        print(step, 't debug',dt,t1,prevlocation,s[t1[0],t1[1]], stepx , stepz)
                s[int(t1[0]),int(t1[1])] = -counter
                counter+=1
            else:
                if (((int(t1[0])!=int(prevlocation[0])) | (int(t1[1])!=int(prevlocation[1]))) & (abs(s[int(t1[0]),int(t1[1])]+counter)>2)):
                    exitreason='loop'
#                    print('loop after ', counter, t1[0],t1[1],prevlocation,s[int(t1[0]),int(t1[1])],centerslen + radius)
                    break
#                        print(step, counter, 't2 debug',dt,t1,prevlocation,s[t1[0],t1[1]], stepx , stepz, int(t1[0]),int(prevlocation[0]),int(t1[1]),int(prevlocation[1]))
        print('trajdebug1')
        if len(distance)>0:
            print('trajdebug2')
            distance = np.asarray(distance)
            dx = distance*.0+1
            dydx = np.gradient(distance, dx)
            dy2dx = np.gradient(dydx, dx)
#        plt.figure()
##        plt.plot(dydx,label='dxdy')
##        plt.plot(dy2dx,label='dxdy2')
#        plt.plot(side2,label='side2')
#        plt.legend()
#        plt.show()
        print('trajdebug3')
        plt.figure()
        plt.plot(xx[2:],yy[2:])
        plt.show()
        print('trajdebug4')
        import math
        a=[5,1]
        a=[0,2]
        b=[2,3]
        c=[3,5]
        a=[0,-2]
        b=[0,2]
        c=[2,0]


        # https://www.wikihow.com/Find-the-Angle-Between-Two-Vectors
        angle=math.degrees(math.acos(((a[0]-b[0])*(b[0]-c[0])+(a[1]-b[1])*(b[1]-c[1]))/((((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5)*(((b[0]-c[0])**2+(b[1]-c[1])**2)**.5))))
        # https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line        
        angle2 = ((b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0]))
#        print('angle is ',angle,angle2)

        a=[0,0]
        b=[2,-1]
        c=[-1,2]

        a=[0,0]
        b=[-1,-2]
        c=[2,1]

        
        ax=a[0]
        ay=a[1]
        bx=b[0]
        by=b[1]
        cx=c[0]
        cy=c[1]
#        print(((bx - ax)*(cy - by) - (by - ay)*(cx - ax)))
        side2=np.asarray(side2)
#        print('side2 max', side2.max())

        if exitreason=='end':
            #s[s<0]=point*100+radius
#            print('end loop', step, counter,t1)
            return 2
        elif exitreason=='boundary':
#            print('boundary', step, counter)
            return 3
        elif exitreason=='mask':
#            print('mask', step, counter)
            return mask[int(t1[0]),int(t1[1])]
        else: # loop
#            print(exitreason)
#            print('finish loop after:',step ,' steps')
            return 1

#        if np.isnan(u[location[0],location[1]]):
#            return np.NaN
#        else:
#            return location[0]+location[1]
        
    def on_click_eddiesold(self):
        textboxtime1value = float(self.itteration1.text())
        textboxfile1value = self.filename1.text()
        
        textboxtime1value=1000
        textboxfile1value = r'/ibdata2/nirb/openFOAM/angle/kleombenchmarkwest/cylinder.foam'

        textboxtime1value=3467000
        textboxfile1value = r'/ibdata2/nirb/openFOAM/ml/windAroundurbanMichelstadt2zomegablock1/windAroundCube.foam'
        textboxtime1value=3467000
        textboxfile1value = r'/ibdata2/nirb/openFOAM/vortex/tlvs2/tlvs2.foam'
        

#        if self.heightcontour.checkState():
#            db = get_slice_height(textboxfile1value, textboxtime1value, int(self.listaxis.currentRow()),
#                                  float(self.axispos.text()),
#                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
#                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
#                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
#                                  reader=self.reader1)
#        else:
#            db = get_slice(textboxfile1value, textboxtime1value, int(self.listaxis.currentRow()),
#                           float(self.axispos.text()),
#                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
#                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
#                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=self.reader1)


        slicey=632000
        clipx1=220500
        clipx2=221000
        clipy1=0        
        clipy2=0       
        clipz1=700       
        clipz2=830     
        
        slicey=-250000
        clipx1=0
        clipx2=0
        clipy1=0        
        clipy2=0       
        clipz1=0       
        clipz2=0     
        
        slicey=664200
        clipx1=178750
        clipx2=179000
        clipy1=0        
        clipy2=0       
        clipz1=0       
        clipz2=50     
        
        db = get_slice(textboxfile1value, textboxtime1value, 1,
                       slicey,
                       clipxmin=clipx1, clipxmax=clipx2,
                       clipymin=clipy1, clipymax=clipy2,
                       clipzmin=clipz1, clipzmax=clipz2)
        
        ticks=400
#        xmin, xmax, ymin, ymax, zi1 = makegrid(db, self.listfields.currentItem().text(), self.listaxis.currentRow())
        xmin, xmax, ymin, ymax, u = makegrid(db, 'U_x', 1, ticks=ticks)
        xmin, xmax, ymin, ymax, w = makegrid(db, 'U_z', 1, ticks=ticks)
#        xmin, xmax, ymin, ymax, ziy = makegrid(db, 'U_y', 1, ticks=400)

##################
        xi = np.linspace(0,ticks,ticks)
        zi = np.linspace(0,ticks,ticks)
        X,Y=np.meshgrid(xi,zi)
        w[np.isnan(w)]=0
        u[np.isnan(u)]=0
        u[np.abs(u)<0.000001]=0
        w[np.abs(w)<0.000001]=0
        intx=integrate.cumtrapz(w,X,axis=1,initial=0)[0]
        inty=integrate.cumtrapz(u,Y,axis=0,initial=0)
        psi1=intx-inty
        w[w==0.]=np.nan
        u[u==0]=np.nan
        
###################
          
#        zi1=((u)**2+(ziy)**2+(w)**2)**0.5
#        zi1=((u)**2+(w)**2)**0.5
#        zi1=np.round(zi1,decimals=6)
        zi1=psi1.copy()
        zi1[u == 0] = np.NaN

# multiple neibour find negative (sign change)
        psi1u=np.roll(psi1,1,axis=0) #comes from above   
        psi1d=np.roll(psi1,-1,axis=0) #comes from below    
        psi1l=np.roll(psi1,1,axis=1) #comes from left
        psi1lu=np.roll(psi1l,1,axis=0) #comes from left and above
        psi1ld=np.roll(psi1l,-1,axis=0) #comes from left and below
        psi1r=np.roll(psi1,-1,axis=1) #comes from right
        psi1ru=np.roll(psi1r,1,axis=0) #comes from right and above
        psi1rd=np.roll(psi1r,-1,axis=0) #comes from right and below
        psi1m=np.zeros_like(psi1)-1
#        psi1m[(psi1>=psi1u) & (psi1>=psi1d) & (psi1>=psi1l) & (psi1>=psi1r) & (psi1>=psi1lu) & (psi1>=psi1ld) & (psi1>=psi1ru) & (psi1>=psi1rd)]=-2 # all max
#        psi1m[(psi1<=psi1u) & (psi1<=psi1d) & (psi1<=psi1l) & (psi1<=psi1r) & (psi1<=psi1lu) & (psi1<=psi1ld) & (psi1<=psi1ru) & (psi1<=psi1rd)]=-3 # all min

#        psi1[np.isnan(u)]=-9999
        psi1m[(psi1>psi1u) & (psi1>psi1d) & (psi1>psi1l) & (psi1>psi1r) & (psi1>psi1lu) & (psi1>psi1ld) & (psi1>psi1ru) & (psi1>psi1rd)]=-2 # all max
#        psi1[np.isnan(u)]=9999
        psi1m[(psi1<psi1u) & (psi1<psi1d) & (psi1<psi1l) & (psi1<psi1r) & (psi1<psi1lu) & (psi1<psi1ld) & (psi1<psi1ru) & (psi1<psi1rd)]=-3 # all min
    
#            psi1m[np.isnan(unan)]==np.NaN                              
#            psi1m[u==0]=np.NaN
        
        print('1-',np.sum(psi1m==-2))
        print(np.where(psi1m==-2))
        print('2-',np.sum(psi1m==-3))
        print(np.where(psi1m==-3))
        suspect1 = np.where((psi1m==-2) | (psi1m==-3))
        print('suspected points: ', len(suspect1[0]), suspect1)
        suspect=[]
        for i in range(len(suspect1[0])):
            suspect.append([suspect1[0][i],suspect1[1][i],psi1[suspect1[0][i],suspect1[1][i]]])
        suspect=np.asarray(suspect)
        suspect = suspect[(suspect[:, 2]).argsort()]
        countindex = suspect
#        for i in range(len(countindex[0])):
#            remove = False
#            for j in range(-1,2):
#                for k in range(-1,2):
#                    if ((countindex[0][i]+j)<0) or ((countindex[1][i]+k)<0) or ((countindex[0][i]+j)>=mulz.shape[0]) or ((countindex[1][i]+k)>=mulz.shape[1]):
#                        remove = True
#                    else:
#                        if mulz[countindex[0][i]+j,countindex[1][i]+k]==0:
#                            remove = True
#            print(remove,countindex[0][i],countindex[1][i])
#            if remove:
#                mulz[countindex[0][i],countindex[1][i]]=1
                
        
        zi2=np.zeros_like(zi1)   
        counter=0
#        zi2[zi1==0]==counter
#        zi2[np.isnan(zi1)]=counter
        zi1[np.isnan(u)]=np.nan
            
########### option 2 start
        for k in range(len(suspect)):
            counter=k+1
            zi2[int(suspect[k][0]),int(suspect[k][1])]=counter
            print (k, 'CCOOUUNNTTEERR',counter, suspect[k][0], suspect[k][1])
            loop = 0 # debug print variable     
            while True:
                loop +=1
                found = False       
                countindex = np.where(zi2==counter)
                print(counter, loop, len(np.where(zi2==-counter)[0]), len(np.where(zi2==counter)[0]),countindex)
                zi2[int(suspect[k][0]),int(suspect[k][1])]=-counter
                for i in range(len(countindex[0])):
                    if countindex[0][i]>0:
                        if zi2[countindex[0][i]-1,countindex[1][i]]==0:
                            if zi1[countindex[0][i]-1,countindex[1][i]]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i]-1,countindex[1][i]]=counter
                                found=True               
                        if countindex[1][i]>0:
                            if zi2[countindex[0][i]-1,countindex[1][i]-1]==0:
                                if zi1[countindex[0][i]-1,countindex[1][i]-1]<=zi1[countindex[0][i],countindex[1][i]]:
                                    zi2[countindex[0][i]-1,countindex[1][i]-1]=counter
                                    found=True
                        if countindex[1][i]<zi1.shape[1]-1:
                            if zi2[countindex[0][i]-1,countindex[1][i]+1]==0:
                                if zi1[countindex[0][i]-1,countindex[1][i]+1]<=zi1[countindex[0][i],countindex[1][i]]:
                                    zi2[countindex[0][i]-1,countindex[1][i]+1]=counter
                                    found=True
                    if countindex[0][i]<zi1.shape[0]-1:
                        if zi2[countindex[0][i]+1,countindex[1][i]]==0:
                            if zi1[countindex[0][i]+1,countindex[1][i]]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i]+1,countindex[1][i]]=counter
                                found=True
                        if countindex[1][i]>0:
                            if zi2[countindex[0][i]+1,countindex[1][i]-1]==0:
                                if zi1[countindex[0][i]+1,countindex[1][i]-1]<=zi1[countindex[0][i],countindex[1][i]]:
                                    zi2[countindex[0][i]+1,countindex[1][i]-1]=counter
                                    found=True
                        if countindex[1][i]<zi1.shape[1]-1:
                            if zi2[countindex[0][i]+1,countindex[1][i]+1]==0:
                                if zi1[countindex[0][i]+1,countindex[1][i]+1]<=zi1[countindex[0][i],countindex[1][i]]:
                                    zi2[countindex[0][i]+1,countindex[1][i]+1]=counter
                                    found=True
                    if countindex[1][i]>0:
                        if zi2[countindex[0][i],countindex[1][i]-1]==0:
                            if zi1[countindex[0][i],countindex[1][i]-1]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i],countindex[1][i]-1]=counter
                                found=True
                    if countindex[1][i]<zi1.shape[1]-1:
                        if zi2[countindex[0][i],countindex[1][i]+1]==0:
                            if zi1[countindex[0][i],countindex[1][i]+1]<=zi1[countindex[0][i],countindex[1][i]]:
                                zi2[countindex[0][i],countindex[1][i]+1]=counter
                                found=True
                                    
                                    
                if found is False:
                    break

########### option 2 end            
            
            
########### option 1 start
#            
#        zi3=zi1.copy()
#        zi3[zi3 == 0] = np.nanmax(zi3)
#        zi3[np.isnan(zi3)] = np.nanmax(zi3)
#        # I didn't do it recursive due to heap stack problem
#        while (len(np.where(zi2==0)[0])>0):        
#            counter += 1
#            zi3[zi2!=0]=np.max(zi3)
#            minindex = np.where((zi3==np.amin(zi3)) & (zi2==0))
#            print ('CCOOUUNNTTEERR',counter, minindex,np.sum(zi2==0))           
#            zi2[minindex[0][0],minindex[1][0]] = counter            
#            while True:
#                found = False       
#                countindex = np.where(zi2==counter)
#                for i in range(len(countindex[0])):
#                    if countindex[0][i]>0:
#                        if zi2[countindex[0][i]-1,countindex[1][i]]==0:
#                            if zi1[countindex[0][i]-1,countindex[1][i]]>=zi1[countindex[0][i],countindex[1][i]]:
#                                zi2[countindex[0][i]-1,countindex[1][i]]=counter
#                                found=True               
#                        if countindex[1][i]>0:
#                            if zi2[countindex[0][i]-1,countindex[1][i]-1]==0:
#                                if zi1[countindex[0][i]-1,countindex[1][i]-1]>=zi1[countindex[0][i],countindex[1][i]]:
#                                    zi2[countindex[0][i]-1,countindex[1][i]-1]=counter
#                                    found=True
#                        if countindex[1][i]<zi1.shape[1]-1:
#                            if zi2[countindex[0][i]-1,countindex[1][i]+1]==0:
#                                if zi1[countindex[0][i]-1,countindex[1][i]+1]>=zi1[countindex[0][i],countindex[1][i]]:
#                                    zi2[countindex[0][i]-1,countindex[1][i]+1]=counter
#                                    found=True
#                    if countindex[0][i]<zi1.shape[0]-1:
#                        if zi2[countindex[0][i]+1,countindex[1][i]]==0:
#                            if zi1[countindex[0][i]+1,countindex[1][i]]>=zi1[countindex[0][i],countindex[1][i]]:
#                                zi2[countindex[0][i]+1,countindex[1][i]]=counter
#                                found=True
#                        if countindex[1][i]>0:
#                            if zi2[countindex[0][i]+1,countindex[1][i]-1]==0:
#                                if zi1[countindex[0][i]+1,countindex[1][i]-1]>=zi1[countindex[0][i],countindex[1][i]]:
#                                    zi2[countindex[0][i]+1,countindex[1][i]-1]=counter
#                                    found=True
#                        if countindex[1][i]<zi1.shape[1]-1:
#                            if zi2[countindex[0][i]+1,countindex[1][i]+1]==0:
#                                if zi1[countindex[0][i]+1,countindex[1][i]+1]>=zi1[countindex[0][i],countindex[1][i]]:
#                                    zi2[countindex[0][i]+1,countindex[1][i]+1]=counter
#                                    found=True
#                if found is False:
#                    break
                
########### option 1 end

               
        plt.figure()
        ax = plt.subplot(1,2,1)
        plt.imshow(u,origin='lower')
        plt.colorbar()
        plt.title('u')

#        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)

        plt.subplot(1,2,2, sharex=ax, sharey=ax)
        plt.imshow(w,origin='lower')
        plt.colorbar()
        plt.title('w')

#        plt.subplot(1,3,3, sharex=ax, sharey=ax)
        
        plt.figure()
        zi1[np.isnan(w)]=np.nan       
        plt.imshow(zi1,origin='lower')
        plt.colorbar()
        plt.title('psi1')

        plt.figure()
        zi2[zi2 == 1] = np.NaN
        X1,Y1 = np.meshgrid(np.linspace(0,u.shape[0],u.shape[0]), np.linspace(0,u.shape[1],u.shape[1]))
        n=10
        zoomx1=0
        zoomx2=X1.shape[0]
        zoomz1=0
        zoomz2=X1.shape[1]
        
        zoomx1=700
        zoomx2=900
        zoomz1=0
        zoomz2=200
        
        zi2[np.isnan(w)]=np.nan
        plt.imshow(zi2,origin='lower')
        plt.colorbar()
        plt.quiver(X1[zoomz1:zoomz2:n,zoomx1:zoomx2:n],Y1[zoomz1:zoomz2:n,zoomx1:zoomx2:n],u[zoomz1:zoomz2:n,zoomx1:zoomx2:n],w[zoomz1:zoomz2:n,zoomx1:zoomx2:n], scale_units='xy', scale=0.1)
        plt.title('flood')

        
        y=83
        y1=43
        y2=123
        nx = zi1.shape[1]
        d2=np.zeros(nx)
        d21=np.zeros(nx)
        d22=np.zeros(nx)
        d23=np.zeros(nx)
        c2=np.zeros(nx)
        c21=np.zeros(nx)
        c22=np.zeros(nx)
        c23=np.zeros(nx)

        for i in range(926,1010): #nx
            for j in range(926,1010):
                if ~np.isnan(zi1[y,i]) & ~np.isnan(zi1[y,j]):
                    c2[math.fabs(i-j)]+=(zi1[y,i]*zi1[y,j])
                    d2[math.fabs(i-j)]+=(zi1[y,i]-zi1[y,j])**2
                if ~np.isnan(zi1[y1,i]) & ~np.isnan(zi1[y1,j]):
                    c21[math.fabs(i-j)]+=(zi1[y1,i]*zi1[y1,j])
                    d21[math.fabs(i-j)]+=(zi1[y1,i]-zi1[y1,j])**2                  
                if ~np.isnan(zi1[y2,i]) & ~np.isnan(zi1[y2,j]):
                    c22[math.fabs(i-j)]+=(zi1[y2,i]*zi1[y2,j])
                    d22[math.fabs(i-j)]+=(zi1[y2,i]-zi1[y2,j])**2                       
        for i in range(43,123): #nx
            for j in range(43,123):
                if ~np.isnan(zi1[i,970]) & ~np.isnan(zi1[i,970]):
                    c23[math.fabs(i-j)]+=(zi1[i,970]*zi1[j,970])
                    d23[math.fabs(i-j)]+=(zi1[i,970]-zi1[j,970])**2
                       
        plt.figure()
        plt.plot(d2,'*r')
        plt.plot(d21,'og')
        plt.plot(d22,'+b')
        plt.plot(d23,'.k')
        plt.title('d2')
        plt.figure()
        plt.plot(c2,'*r')
        plt.plot(c21,'og')
        plt.plot(c22,'+b')
        plt.plot(c23,'.k')
        plt.title('c2')
        plt.figure()
        plt.plot(zi1[y,:],'*r')
        plt.plot(u[y,:],'+g')
        plt.title('wind750')
        plt.figure()
        plt.plot(zi1[y1,:],'*r')
        plt.plot(u[y1,:],'+g')
        plt.title('wind700')
        plt.figure()
        plt.plot(zi1[y2,:],'*r')
        plt.plot(u[y2,:],'+g')
        plt.title('wind800')
        
        
# https://stackoverflow.com/questions/49557329/compute-stream-function-from-x-and-y-velocities-by-integration-in-python        
#        from scipy import integrate      
#
#        # integrate to make an intial guess
#        intx=integrate.cumtrapz(w,X1,axis=1,initial=0)[0]
#        inty=integrate.cumtrapz(u,Y1,axis=0,initial=0)
#        psi1=intx-inty
#        
#        intx=integrate.cumtrapz(w,X1,axis=1,initial=0)
#        inty=integrate.cumtrapz(u,Y1,axis=0,initial=0)[:,0][:,None]
#        psi2=intx-inty
#        
#        psi=0.5*(psi1+psi2)
#        psi[np.isnan(zi2)] = np.NaN
#
#        plt.figure()
#        plt.imshow(-psi,origin='lower') #[0:200,:]
#        plt.colorbar()
#        plt.quiver(X1[zoomz1:zoomz2:n,zoomx1:zoomx2:n],Y1[zoomz1:zoomz2:n,zoomx1:zoomx2:n],u[zoomz1:zoomz2:n,zoomx1:zoomx2:n],w[zoomz1:zoomz2:n,zoomx1:zoomx2:n], scale_units='xy', scale=0.1)
#        plt.title('stream function')
        
        
        plt.figure()  
        ax = plt.subplot(2,2,1)
        plt.imshow(zi4x,origin='lower')
        plt.colorbar()
        plt.title('Ux')
        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
        plt.subplot(2,2,2, sharex=ax, sharey=ax)
        plt.imshow(zi4z,origin='lower')
        plt.colorbar()
        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
        plt.title('Uz')
        plt.subplot(2,2,3, sharex=ax, sharey=ax)
#        plt.imshow(zi2,origin='lower')
        plt.imshow(zi5x,origin='lower')
        plt.colorbar()
        plt.title('zi2')
        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
        
        plt.subplot(2,2,4, sharex=ax, sharey=ax)
        plt.imshow(zi5z,origin='lower')
        plt.colorbar()
        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
        plt.title('mulz')        


        plt.figure()
  
        ax = plt.subplot(2,2,1)
        plt.imshow(u,origin='lower')
        plt.colorbar()
        plt.title('Ux')
        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
        
        plt.subplot(2,2,2, sharex=ax, sharey=ax)
        plt.imshow(w,origin='lower')
        plt.colorbar()
        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
        plt.title('Uz')

        plt.subplot(2,2,3, sharex=ax, sharey=ax)
#        plt.imshow(zi2,origin='lower')
        plt.imshow(mulx,origin='lower')
        plt.colorbar()
        plt.title('zi2')
#        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
        
        plt.subplot(2,2,4, sharex=ax, sharey=ax)
        plt.imshow(mulz,origin='lower')
        plt.colorbar()
#        plt.quiver(X1[::n,::n],Y1[::n,::n],u[::n,::n],w[::n,::n], scale_units='xy', scale=0.1)
        plt.title('mulz')

        
    def on_click_learn(self):

        full = True
        files=[]   
        inside=True
        if not(inside):
            listItems=self.listfiles.selectedItems()      
            if listItems:           
                for item in listItems:
                   files.append(item.text())
            else:
                for item in range(self.listfiles.count()):
                    files.append(self.listfiles.item(item).text())
            textboxtime1value = float(self.itteration1.text())
            db = get_clip(textboxfile1value, [textboxtime1value],
                      clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                      clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                      clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), parallel=False)
            xmin = float(self.xmin.text())
            xmax = float(self.xmax.text())
            ymin = float(self.ymin.text())
            ymax = float(self.ymax.text())
            zmin = float(self.zmin.text())
            zmax = float(self.zmax.text())
            learnfrom=self.learnfrom.text()
            
        else:
            # postProcess -func writeCellCentres -time 0
            files=[r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor1/windAroundcaseE.foam']
            files=[r'/data4bk/nirb/Simulations/michaelstadtfloor1ml/windAroundcaseE.foam']
            files=[r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloornoborders/windAroundcaseE.foam']
            textboxtime1value=1000
            textboxtime1value=2300

            files=[r'/data4bk/nirb/Simulations/Dans/tlvsmall/tlvs2.foam']
            textboxtime1value=3800

            files=[r'/data4bk/nirb/Simulations/michaelstadtfloornoborders/windAroundcaseE.foam']
            textboxtime1value=10000


            textboxfile1value = files[0]
            filename=files[0]                       
            db = get_clip(textboxfile1value, [textboxtime1value],
                  clipxmin=0, clipxmax=0,
                  clipymin=0, clipymax=0,
                  clipzmin=0, clipzmax=0, parallel=True)
            xmin=0
            xmax=0
            ymin=0
            ymax=0
            zmin=0
            zmax=0
            learnfrom=''


        textboxfile1value = files[0]

        learnfile = shrinktitle(textboxfile1value)
        print('start learn',textboxfile1value, datetime.datetime.now())


        if usehera:
        # ap3 = db[['x', 'y', 'z', sel[0].data()]]  # at the ground (1 meter)
            if full:
                ap3 = db[['x', 'y', 'z', 'U_x', 'U_y', 'U_z','p','k','nut','epsilon']]  # at the ground (1 meter)
            else:
                ap3 = db[['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)            
        else:
            # ap3 = db[['x', 'y', 'z', sel[0].data()]]  # at the ground (1 meter)
            if full:
                ap3 = db[0][['x', 'y', 'z', 'U_x', 'U_y', 'U_z','p','k','nut','epsilon']]  # at the ground (1 meter)
            else:
                ap3 = db[0][['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)
        points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
        # values = ap3[sel[0].data()]
        valuesux = ap3['U_x']
        valuesuy = ap3['U_y']
        valuesuz = ap3['U_z']
        if full:
            valuesp = ap3['p']
            valuesk = ap3['k']
            valueso = ap3['epsilon']
            valuesn = ap3['nut']

        if xmin==xmax:
            xmin=ap3['x'].min()
            xmax=ap3['x'].max()
        if ymin==ymax:
            ymin=ap3['y'].min()
            ymax=ap3['y'].max()
        if zmin==zmax:
            zmin=ap3['z'].min()
            zmax=ap3['z'].max()
        # ticks = 600j
        grid_points = 1800000000 # tlv 1.7 1800000000 # 27000000
        meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)

        resolution = (meters / grid_points)**(1./3)
        if resolution==0:
            print('no clip boundaries in one of the axis')
        # resolution = 0.1
        ticksx = int((xmax-xmin) / resolution) * 1j
        ticksy = int((ymax-ymin) / resolution) * 1j
        ticksz = int((zmax-zmin) / resolution) * 1j
        print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)

        xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
        # grid = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(values), (xi, yi, zi),
        #                    method='nearest')  # Nearest for keeping zeros that are buildings
        gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        print('grid1', datetime.datetime.now())
        griduy = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuy), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        print('grid2', datetime.datetime.now())
        griduz = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuz), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        print('grid3', datetime.datetime.now())
        if full:
            gridrp = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesp), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            print('gridp',datetime.datetime.now())
            gridrk = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesk), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            print('gridk',datetime.datetime.now())
            gridro = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valueso), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            print('grido',datetime.datetime.now())
            gridrn = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesn), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            print('gridn',datetime.datetime.now())

        x = np.linspace(ap3['x'].min(), ap3['x'].max(), int(ticksx.imag))
        y = np.linspace(ap3['y'].min(), ap3['y'].max(), int(ticksy.imag))
        z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticksz.imag))

        if 5==5:
            print('start to read U file in the time directory')
            slash = textboxfile1value.rfind("/") + 1
            directory = textboxfile1value[:slash]+ str(int(float(textboxtime1value))) +"/"
            directory = r'/data4bk/nirb/Simulations/Haifa/haifa32p1/58000/'
            fileu = directory + "U"
            filex = directory + "Cx"
            filey = directory + "Cy"
            filez = directory + "Cz"
            if full:
                filep = directory + "p"
                filek = directory + "k"
                fileo = directory + "epsilon"
                filen = directory + "nut"
            fhu = open(fileu, 'r')
            fhx = open(filex, 'r')
            fhy = open(filey, 'r')
            fhz = open(filez, 'r')
            if full:
                fhp = open(filep, 'r')
                fhk = open(filek, 'r')
                fho = open(fileo, 'r')
                fhn = open(filen, 'r')
            txtu = fhu.read()
            txtx = fhx.read()
            txty = fhy.read()
            txtz = fhz.read()
            if full:
                txtp = fhp.read()
                txtk = fhk.read()
                txto = fho.read()
                txtn = fhn.read()
            txtusplit = txtu.split()
            txtxsplit = txtx.split()
            txtysplit = txty.split()
            txtzsplit = txtz.split()
            if full:
                txtpsplit = txtp.split()
                txtksplit = txtk.split()
                txtosplit = txto.split()
                txtnsplit = txtn.split()
            posu = txtusplit.index("internalField")
            posx = txtxsplit.index("internalField")
            posy = txtysplit.index("internalField")
            posz = txtzsplit.index("internalField")
            if full:
                posp = txtpsplit.index("internalField")
                posk = txtksplit.index("internalField")
                poso = txtosplit.index("internalField")
                posn = txtnsplit.index("internalField")
            items = int(txtusplit[posu+3])
            # first item is posu+5

            us = []
            xs = []
            ys = []
            zs = []
            if full:
                ps = []
                ks = []
                os = []
                ns = []
            print('start reading U field', datetime.datetime.now())
            for i in range(items):
                us.append([float(remove_left(txtusplit[posu+5+i*3])), float(txtusplit[posu+6+i*3]), float(remove_right(txtusplit[posu+7+i*3]))])
                xs.append(float(txtxsplit[posx+5+i]))
                ys.append(float(txtysplit[posy + 5 + i]))
                zs.append(float(txtzsplit[posz + 5 + i]))
                if full:
                    ps.append(float(txtpsplit[posp + 5 + i]))
                    ks.append(float(txtksplit[posk + 5 + i]))
                    os.append(float(txtosplit[poso + 5 + i]))
                    ns.append(float(txtnsplit[posn + 5 + i]))

            xs = np.asarray(xs)
            ys = np.asarray(ys)
            zs = np.asarray(zs)
            if full:
                ps = np.asarray(ps)
                ks = np.asarray(ks)
                os = np.asarray(os)
                ns = np.asarray(ns)

            print('finish reading U field', datetime.datetime.now())
            gridu0 = griddata((xs, ys, zs), np.asarray(us)[:,0], (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
            gridu1 = griddata((xs, ys, zs), np.asarray(us)[:,1], (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
            gridu2 = griddata((xs, ys, zs), np.asarray(us)[:,2], (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
            gridu0[gridux==0]=0.0
            gridu1[griduy==0]=0.0
            gridu2[griduz==0]=0.0
            
            
            if full:
                gridp = griddata((xs, ys, zs), np.asarray(ps), (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
                gridk = griddata((xs, ys, zs), np.asarray(ks), (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
                grido = griddata((xs, ys, zs), np.asarray(os), (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
                gridn = griddata((xs, ys, zs), np.asarray(ns), (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings

            plt.figure()
            testxnp = gridu0.ravel()
            sc = plt.scatter(gridu0.ravel(), gridux.ravel(), c=(gridu0.ravel()-gridux.ravel()), s=1) #, cmap=cm)
            plt.colorbar(sc)
            plt.xlabel('time directory ordered file')
            plt.ylabel('grid vtk')
            plt.title('READING file in order ? ')
            plt.show()
            print('corr for reading is:', str(np.corrcoef(gridu0.ravel(), gridux.ravel())[0, 1]),str(sklearn.metrics.r2_score(gridu0.ravel(), gridux.ravel())))

        if 5==5:
            np.save('mgridu0',gridu0)
            np.save('mgridu1',gridu1)
            np.save('mgridu2',gridu2)

            gridu0=gridux
            gridu1=griduy
            gridu2=griduz
            grid_z = np.zeros_like(gridu0)  # above sea level
            grid_y = np.zeros_like(gridu0)  # lat
            grid_x = np.zeros_like(gridu0)  # lon
            grid_h = np.zeros_like(gridu0)  # above ground (including building)
            grid_t = np.zeros_like(gridu0)  # above topograpyh
            grid_bl = np.zeros_like(gridu0)  # distance from boundary level
            grid_dp = np.zeros_like(gridu0)  # distance from previous building
            grid_dn = np.zeros_like(gridu0)  # distance from next building
            grid_dl = np.zeros_like(gridu0)  # distance between left
            grid_dr = np.zeros_like(gridu0)  # distance between right
            grid_hp = np.zeros_like(gridu0)  # height of previous building
            grid_hn = np.zeros_like(gridu0)  # height of previous building
            grid_hl = np.zeros_like(gridu0)  # height of previous building
            grid_hr = np.zeros_like(gridu0)  # height of previous building

            np.save('nbgrid_z',grid_z)
            print('grid_z', datetime.datetime.now())
            np.save('nbgrid_y',grid_y)
            print('grid_y', datetime.datetime.now())
            np.save('nbgrid_x',grid_x)
            print('grid_x', datetime.datetime.now())
            np.save('nbgrid_h',grid_h)
            print('grid_h', datetime.datetime.now())
            np.save('nbgrid_t',grid_t)
            print('grid_t', datetime.datetime.now())
            np.save('nbgrid_dp',grid_dp)
            print('grid_dp', datetime.datetime.now())
            np.save('nbgrid_dn',grid_dn)
            print('grid_dn', datetime.datetime.now())
            np.save('nbgrid_dl',grid_dl)
            print('grid_dl', datetime.datetime.now())
            np.save('nbgrid_dr',grid_dr)
            print('grid_dr', datetime.datetime.now())
            np.save('nbgrid_hp',grid_hp)
            print('grid_hp', datetime.datetime.now())
            np.save('nbgrid_hn',grid_hn)
            print('grid_hn', datetime.datetime.now())

            np.save('nbx.npy',x)
            np.save('nby.npy',y)
            np.save('nbz.npy',z)
            np.save('nbgridux.npy',gridux)
            np.save('nbgriduy.npy',griduy)
            np.save('nbgriduz.npy',griduz)
            
            from ml import ml

            import numpy as np
            import math
            import datetime
            x = np.load('nbx.npy')
            y = np.load('nby.npy')
            z = np.load('nbz.npy')
            gridux = np.load('nbgridux.npy')
            griduy = np.load('nbgriduy.npy')
            griduz = np.load('nbgriduz.npy')
            gridp=   np.load('nbp.npy')
            gridk=   np.load('nbk.npy')
            gridu0=gridux
            gridu1=griduy
            gridu2=griduz
#            gridu0 = np.load('nbgridu0.npy')
#            gridu1 = np.load('nbgridu1.npy')
#            gridu2 = np.load('nbgridu2.npy')
            
            grid_z = np.load('nbgrid_z.npy')
            grid_y = np.load('nbgrid_y.npy')
            grid_x = np.load('nbgrid_x.npy')
            grid_t = np.load('nbgrid_t.npy')
            grid_h = np.load('nbgrid_h.npy')
            grid_dp = np.load('nbgrid_dp.npy')
            grid_dn = np.load('nbgrid_dn.npy')
            grid_dl = np.load('nbgrid_dl.npy')
            grid_dr = np.load('nbgrid_dr.npy')
            grid_hp = np.load('nbgrid_hp.npy')
            grid_hn = np.load('nbgrid_hn.npy')
            grid_hr = np.load('nbgrid_hr.npy')
            grid_hl = np.load('nbgrid_hl.npy')

            features = np.load('nbfeatures.npy')


            # Todo: work with vetors and not matrics (get rid of the interpolation...)
            print('start loop   ', datetime.datetime.now())
            print('to add features: z / last house height !!!')
#            print('to  :',gridux[5, 5, :100])
            dx = x[2]-x[1]
            dy = y[2]-y[1]
            dz = z[2]-z[1]

            for i in range(0, gridu0.shape[0]): # 445
                print('ij',i)
                if i % 20 == 0:
                    print ('                      i', datetime.datetime.now())                    
                    np.save('nbgrid_x',grid_x)
                    np.save('nbgrid_y',grid_y)
                    np.save('nbgrid_z',grid_z)
                    np.save('nbgrid_h',grid_h)
                    np.save('nbgrid_dp',grid_dp)
                    np.save('nbgrid_dn',grid_dn)
                    np.save('nbgrid_dl',grid_dl)
                    np.save('nbgrid_dr',grid_dr)
                    
                for j in range(gridu0.shape[1]):
                    for k in range(gridu0.shape[2]):
                        grid_x[i,j,k] = x[i]
                        grid_y[i,j,k] = y[j]
                        grid_z[i,j,k] = z[k]
                        path = gridu0[i, j, :k]
                        if len(np.where(path == 0)[0])==0:
                            grid_h[i, j, k] = k * dz
                        else:
#                            print('test_h',i,j, np.where(path == 0)[0], np.where(path == 0)[0][-1])
                            grid_h[i, j, k] = (k - np.where(path == 0)[0][-1])*dz
                        path = gridu0[i:, j, k]
                        if len(np.where(path == 0)[0])==0:
                            grid_dp[i, j, k] = (gridu0.shape[0])*dx  # -i
                        else:
                            grid_dp[i, j, k] = (np.where(path == 0)[0][0])*dx

                        path = np.flipud(gridu0[:i, j, k])
                        if len(np.where(path == 0)[0])==0:
                            grid_dn[i, j, k] = (gridu0.shape[0])*dx
                        else:
                            grid_dn[i, j, k] = (np.where(path == 0)[0][0])*dx

                        path = gridu0[i, j:, k]
                        if len(np.where(path == 0)[0])==0:
                            grid_dl[i, j, k] = (gridu0.shape[1])*dy
                        else:
                            grid_dl[i, j, k] = (np.where(path == 0)[0][0])*dy

                        path = np.flipud(gridu0[i, :j, k])
                        if len(np.where(path == 0)[0])==0:
                            grid_dr[i, j, k] = (gridu0.shape[1])*dy
                        else:
                            grid_dr[i, j, k] = (np.where(path == 0)[0][0])*dy

            print('end loop1     ', datetime.datetime.now())
            
            
#            def fillgrid(i,dx,, gridu0,grid_h)
            
            
            for i in range(0,gridu0.shape[0]): # 2230
                ipp = min(gridu0.shape[0],i+50)
                inn = max(0,i-50)
                if i % 10 == 0:
                    print ('                      i',inn,i,ipp, gridu0.shape[0], datetime.datetime.now())                    
                    np.save('nbgrid_t',grid_t)
                    np.save('nbgrid_hp',grid_hp)
                    np.save('nbgrid_hn',grid_hn)
                    np.save('nbgrid_hl',grid_hp)
                    np.save('nbgrid_hr',grid_hn)
                for j in range(gridu0.shape[1]):
                    jpp = min(gridu0.shape[1],j+50)
                    jnn = max(0,j-50)
#                    print ('j',jnn,j,jpp)
                    for k in range(gridu0.shape[2]):
                        grid_t[i,j,k] = np.min(grid_h[inn:ipp,jnn:jpp,k])

                        path = grid_h[:i, j, k]
#                        path = path - grid_h[i, j, k]
                        if len(np.where(path == 0.)[0])==0:
                            grid_hn[i, j, k] = 3.
                        else:
                            place = np.where(path == 0)[0][-1]
                            colmn = 15 # grid_h[place,j,grid_h.shape[2]-5]
                            grid_hn[i, j, k] = colmn
                        
                        path = np.flipud(grid_h[i:, j, k])
                        path = path - grid_h[i, j, k]
                        if len(np.where(path != 0)[0])==0:
                            grid_hp[i, j, k] = i*dx
                        else:
                            grid_hp[i, j, k] = (np.where(path != 0)[0][0])*dx
                        
                        grid_t[i,j,k] = np.min(grid_h[inn:ipp,jnn:jpp,k])
####
                        path = grid_h[i, :j, k]
#                        path = path - grid_h[i, j, k]
                        if len(np.where(path == 0.)[0])==0:
                            grid_hr[i, j, k] = 3.
                        else:
                            place = np.where(path == 0)[0][-1]
                            colmn = 15 # grid_h[place,j,grid_h.shape[2]-5]
                            grid_hr[i, j, k] = colmn
                        
                        path = np.flipud(grid_h[i, j:, k])
                        path = path - grid_h[i, j, k]
                        if len(np.where(path != 0)[0])==0:
                            grid_hl[i, j, k] = i*dx
                        else:
                            grid_hl[i, j, k] = (np.where(path != 0)[0][0])*dx
                       
                        
            print('end loop2     ', datetime.datetime.now())
            
            angle = np.abs(grid_dl+grid_dr).ravel() / np.abs(grid_dp+grid_dn).ravel()
#            angle[angle>1.]= 1./angle[angle>1.]
            angle = np.arctan(angle)
            angle[np.isinf(angle)] = math.pi/2
            angle[np.isnan(angle)] = 0
#            angle[angle>math.pi/2] = math.pi/2
#            angle[angle<-math.pi/2] = -math.pi/2
            # features = [grid_x.ravel(), grid_y.ravel(), grid_z.ravel(), grid_dp.ravel(), grid_dn.ravel(), np.abs(grid_dl-grid_dr).ravel(), grid_dl.ravel(), grid_dr.ravel(), angle]
            
                          
            plt.figure()
            plt.subplot(2,4,1)
            plt.imshow(grid_h[:,:,30])
            plt.title('h')
            plt.colorbar()
            plt.subplot(2,4,2)
            plt.imshow(grid_t[:,:,30])
            plt.colorbar()
            plt.title('t')
            plt.subplot(2,4,3)
            plt.imshow(grid_dn[:,:,4]+grid_dp[:,:,4])
            plt.colorbar()
            plt.title('dn+dp')
            plt.subplot(2,4,4)
            plt.imshow(grid_dl[:,:,4]+grid_dr[:,:,4])
            plt.colorbar()
            plt.title('dl+dr')
            plt.subplot(2,4,5)
            plt.imshow(grid_dr[:,:,4])
            plt.colorbar()
            plt.title('dr')
            plt.subplot(2,4,6)
            plt.imshow(grid_hp[:,:,4])
            plt.colorbar()
            plt.title('hp')
            plt.subplot(2,4,7) #, sharex=True, sharey=True
            plt.imshow(grid_hn[:,:,4])
            plt.colorbar()
            plt.title('hn')
            plt.subplot(2,4,8)
            plt.imshow(angle.reshape(grid_h.shape[0], grid_h.shape[1],grid_h.shape[2])[:,:,4])
            plt.colorbar()
            plt.title('angle')
            plt.show()                          
                          
            features = [grid_h.ravel(), 
#                        grid_t.ravel(),              new
                        grid_dp.ravel(),
                        grid_dn.ravel(), 
#                        np.abs(grid_dp+grid_dn).ravel(), new
#                        np.abs(grid_dl+grid_dr).ravel(),  new
                        grid_dl.ravel(), 
                        grid_dr.ravel(), 
#                        grid_hp.ravel(), 
#                        grid_hn.ravel(), 
                        angle
                        ]
            
            lim=100.
            dplim=grid_dp.ravel()
            dplim[dplim>dx*lim]=dx*lim
            dnlim=grid_dn.ravel()
            dnlim[dnlim>dx*lim]=dx*lim
            dllim=grid_dl.ravel()
            dllim[dllim>dy*lim]=dy*lim
            drlim=grid_dr.ravel()
            drlim[drlim>dy*lim]=dy*lim
            dzlim=grid_z.ravel()
            dzlim[dzlim>dz*lim]=dz*lim
            
            allfeatures = [
#                        grid_x.ravel(),  #0
#                        grid_y.ravel(),  #1
#                        grid_z.ravel(),      #1
                        dplim, #4
                        dnlim, #4
                        dllim, #4
                        drlim, #4
                        grid_h.ravel(),  #2
                        grid_dl.ravel()/grid_hl.ravel(),  #2
                        grid_dr.ravel()/grid_hr.ravel(),  #3
                        grid_dp.ravel(), #4
                        grid_dp.ravel()/grid_hp.ravel(), #4
                        grid_dn.ravel(), #5
                        np.abs(grid_dp+grid_dn).ravel(),  #6
#                        np.abs(grid_dp+grid_dn).ravel()/grid_hp.ravel(),  # W/H
                        grid_dl.ravel(),  #6
                        grid_dl.ravel()/grid_hl.ravel(),
                        grid_dr.ravel(),  #7
                        np.abs(grid_dl+grid_dr).ravel(),  #7
                        grid_hp.ravel(),  #10
                        grid_hl.ravel(),  #10
                        grid_hr.ravel(),  #10
#                        grid_hn.ravel(),  #11
                        angle] #12

            # normdp = preprocessing.normalize([dplim])[0]
            # normdn = preprocessing.normalize([dnlim])[0]
            # normdl = preprocessing.normalize([dllim])[0]
            # normdr = preprocessing.normalize([drlim])[0]
            # normdz = preprocessing.normalize([dzlim])[0]

            # allfeatures = [normdp,
            #                normdn,
            #                normdl,
            #                normdr,
            #                normdz
            #                ]
                           

#            allfeatures[np.isinf(allfeatures)]=0 # W/H==0

            allfeatures[5][np.isnan(allfeatures[5])]=0.
            allfeatures[8][np.isnan(allfeatures[8])]=0.
            allfeatures[12][np.isnan(allfeatures[12])]=0.
            allfeatures[5][np.isinf(allfeatures[5])]=0.
            allfeatures[8][np.isinf(allfeatures[8])]=0.
            allfeatures[12][np.isinf(allfeatures[12])]=0.
            for i in range(len(allfeatures)):
                print(i,allfeatures[i].min(),allfeatures[i].max())
            
            
            af=np.zeros([len(dplim), 15])
            af[:,0]=dplim
            for i in range(len(af)//10000):
                print(i,len(af)//10000)
                af[i*10000:(i+1)*10000,0]=angle[i*10000:(i+1)*10000]
            af[-10000:-1,0]=angle[-100000:(i+1)*10000]
            
            features = np.asarray(features).T
            features = np.asarray(allfeatures).T
            
            # scaler = preprocessing.StandardScaler().fit(features)            
            # features0 = scaler.transform(features[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]) #
            # features1 = scaler.transform(features[(grid_y.ravel() <= part1y) | (grid_y.ravel() >= part2y)])

            # mlx.features_coeff(scaler.transform(features0[:,:]), labelsux0, featurestest = scaler.transform(features1), labeltest = labelsux1) #
            # mlx.features_coeff(features0[:,:], labelsux0, featurestest = features1, labeltest = labelsux1) #
            

            allfeatures = np.asarray(allfeatures).T
            np.save('allfeatures',allfeatures)
            print ('grid_x1',grid_x.min(), grid_x.max())
            np.save(learnfile+'x', grid_x)
            np.save(learnfile+'y', grid_y)
            np.save(learnfile+'z', grid_z)
            np.save(learnfile+'ux', gridu0)
            np.save(learnfile+'allfeatures', allfeatures)
            np.save(learnfile+'features', features)
                        
            indexprint = 5
           
            featureslevel = features[grid_z.ravel()==z[indexprint]]
            print('indexprint', len(grid_h.ravel()), ticksz, len(featureslevel))
            
            np.save('nbaf.npy',features)
            
            miny = grid_y.ravel().min()
            maxy = grid_y.ravel().max()
#            maxz = grid_z.ravel().max()
            part1y = miny + (maxy-miny) * 0.25# 25
            part2y = miny + (maxy-miny) * 0.5
            part1y1 = miny + (maxy-miny) * 0.5# 25
            part2y1 = miny + (maxy-miny) * 0.75
#            part1y = miny + (maxy-miny) * 0.4
#            part2y = miny + (maxy-miny) * 0.6

            
            print('features0', datetime.datetime.now(), miny, maxy, part1y, part2y)
#            print('features03', datetime.datetime.now(), grid_y.ravel() )
#            print('features05', datetime.datetime.now(), grid_y.ravel() > part1y)
#            print('features06', datetime.datetime.now(), (grid_y.ravel() > part1y) & (grid_y.ravel() < part2y))
#            print('features065', len(features), len(grid_y.ravel()), len(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y), (grid_y.ravel() > part1y) & (grid_y.ravel() < part2y))
            
            features0 = features[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
#            print('features07', datetime.datetime.now(), (grid_y.ravel() > part1y) & (grid_y.ravel() < part2y))

                                  
#            features1 = features[(grid_y.ravel() <= part1y) | (grid_y.ravel() >= part2y)]
            features1 = features[(grid_y.ravel() > part1y1) & (grid_y.ravel() < part2y1)]

            # features0 = scaler.transform(features[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]) #
            # features1 = scaler.transform(features[(grid_y.ravel() <= part1y) | (grid_y.ravel() >= part2y)])
                                  
            print('features1', len(features), len(features0), len(features1))
                                  
            labelsux = gridu0.ravel()
            labelsuy = gridu1.ravel()
            labelsuz = gridu2.ravel()
            labelsk  = gridk.ravel()
            labelsp  = gridp.ravel()
            labelsu2 = (gridu0.ravel()**2+gridu1.ravel()**2)**0.5
            labelsu2 = (gridu0.ravel()**2+gridu1.ravel()**2+gridu2.ravel()**2)**0.5

            labelsu20 = labelsu2[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelsu21 = labelsu2[(grid_y.ravel() <= part1y) | (grid_y.ravel() >= part2y)]
            labelsu21 = labelsu2[(grid_y.ravel() > part1y1) & (grid_y.ravel() < part2y1)]
            
            labelsk0 = labelsk[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelsp0 = labelsp[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelsk1 = labelsk[(grid_y.ravel() > part1y1) & (grid_y.ravel() < part2y1)]
            labelsp1 = labelsp[(grid_y.ravel() > part1y1) & (grid_y.ravel() < part2y1)]

            labelsux0 = labelsux[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelsux1 = labelsux[(grid_y.ravel() <= part1y) | (grid_y.ravel() >= part2y)]
            labelsuy0 = labelsuy[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelsuy1 = labelsuy[(grid_y.ravel() <= part1y) | (grid_y.ravel() >= part2y)]
            labelsuz0 = labelsuz[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
            labelsuz1 = labelsuz[(grid_y.ravel() <= part1y) | (grid_y.ravel() >= part2y)]

            if full:
                labelsp = gridp.ravel()
                labelsk = gridk.ravel()
                labelso = grido.ravel()
                labelsn = gridn.ravel()
                
#            plt.figure()
#            plt.scatter(grid_z[5,5,:],grid_h[5,5,:], s=1, c='b')
#            plt.xlabel('z')
#            plt.ylabel('h')
#            plt.title('H vs Z features in flat city should be the same')
#            plt.show()

#                
            features0[np.isnan(features0)]=0.0
            features1[np.isnan(features1)]=0.0
            features0[np.isinf(features0)]=0.0
            features1[np.isinf(features1)]=0.0
            
            
            if 5==6:
                from ml import ml
                import numpy as np 
                allfeatures=np.load('allfeatures.npy') 
                gridux = np.load('gridux.npy')
                griduy = np.load('griduy.npy')
                griduz = np.load('griduz.npy')
                gridu0=gridux
                gridu1=griduy
                gridu2=griduz
                grid_z = np.load('grid_z.npy')
                grid_y = np.load('grid_y.npy')
                grid_x = np.load('grid_x.npy')

                miny = grid_y.ravel().min()
                maxy = grid_y.ravel().max()
#                maxz = grid_z.ravel().max()
                part1y = miny + (maxy-miny) * 0.25# 25
                part2y = miny + (maxy-miny) * 0.5
                part1y1 = miny + (maxy-miny) * 0.5# 25
                part2y1 = miny + (maxy-miny) * 0.75
    
                features = np.asarray(allfeatures).T

                features0 = features[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
                features1 = features[(grid_y.ravel() > part1y1) & (grid_y.ravel() < part2y1)]

                labelsux = gridu0.ravel()
                labelsuy = gridu1.ravel()
                labelsuz = gridu2.ravel()
                labelsu2 = (gridu0.ravel()**2+gridu1.ravel()**2+gridu2.ravel()**2)**0.5
    
                labelsu20 = labelsu2[(grid_y.ravel() > part1y) & (grid_y.ravel() < part2y)]
                labelsu21 = labelsu2[(grid_y.ravel() > part1y1) & (grid_y.ravel() < part2y1)]
               
                features0[np.isnan(features0)]=0.0
                features1[np.isnan(features1)]=0.0
                features0[np.isinf(features0)]=0.0
                features1[np.isinf(features1)]=0.0
                features=features0
                labels=labelsu20
                show='tU3'
                featurestest = features1
                labeltest = labelsu21                     
            if (learnfrom==''): 
                ml2 = ml()
                clf2, scaler, score = ml2.fit(features0, labelsu20, show='U3', featurestest = features1, labeltest = labelsu21)
                ml2.save('ml3-'+learnfile)

                u2ml = ml2.predict(features1)







                jumps=len(features)//2000000
                labelsux=labelsux[::jumps]
                labelsuy=labelsuy[::jumps]
                labelsuz=labelsuz[::jumps]
                features=features[::jumps]


                mlx = ml()

#                print ('start feature selection', datetime.datetime.now())
#                mlx.features_coeff(features0, labelsux0, featurestest = features1, labeltest = labelsux1)
#                mlx.features_coeff(features0, labelsuy0, featurestest = features1, labeltest = labelsuy1)
#                mlx.features_coeff(features0, labelsuz0, featurestest = features1, labeltest = labelsuz1)
#                print ('end feature selection', datetime.datetime.now())
    
                print('mlx1', datetime.datetime.now())           
                learnfile=''
                clfx, scaler, score = mlx.fit(features, labelsux, show='Ux', featurestest = None, labeltest = None)
                learnfile='nb'
#                clfx, scaler, score = mlx.fit(features0, labelsux0, show='Ux')#, featurestest = features1, labeltest = labelsux1)
                print('mlx2', datetime.datetime.now())
#                print('mlx3', learnfile, datetime.datetime.now())
                mlx.save('mlx-'+learnfile)
                print('>>>UX', score, clfx)
                mly = ml()
#                clfy, scaler, score = mly.fit(features, labelsuy, show='Uy')
                learnfile=''
                clfy, scaler, score = mly.fit(features, labelsuy, show='Uy')#, featurestest = features1, labeltest = labelsuy1)
                learnfile='nb'
                mly.save('mly-'+learnfile)
                print('>>>UY', score, clfy)
                mlz = ml()
#                clfz, scaler, score = mlz.fit(features, labelsuz, show='Uz')
                learnfile=''
                clfz, scaler, score = mlz.fit(features, labelsuz, show='Uz')#, featurestest = features1, labeltest = labelsuz1)
                learnfile='nb'
                mlz.save('mlz-'+learnfile)
                print('>>>UZ', score, clfz)
                if full:
                    mlp = ml()
                    clfp, scaler, score = mlp.fit(features, labelsp)
                    print('>>>P', score, clfp)
                    mlk = ml()
                    clfk, scaler, score = mlk.fit(features, labelsk)
                    print('>>>K', score, clfk)
                    mlo = ml()
                    clfo, scaler, score = mlo.fit(features, labelso)
                    print('>>>Omega', score, clfo)
                    mln = ml()
                    clfn, scaler, score = mln.fit(features, labelsn)
                    print('>>>nu', score, clfn)
                
            else:    
                 print ('before 2 cases score test')
#            if (self.learnfrom.text()!=''): 
                 print('loading:',self.learnfrom.text())
                 mlx = ml()
                 mlx.load('mlx-'+self.learnfrom.text())
                 print('x loaded')
                 mly = ml()
                 mly.load('mly-'+self.learnfrom.text())
                 print('y loaded')
                 mlz = ml()
                 mlz.load('mlz-'+self.learnfrom.text())
                 print('z loaded', len(labelsux))
                 rndmsk = np.random.rand(len(labelsux)) < 0.01
                 print('f010:',features[rndmsk])
                 u = mlx.predict(features[rndmsk])
#                 u = mlx.predict(features)
                 print('upredict',u)

                 corr = np.corrcoef(labelsux[rndmsk], u)[0, 1]
                 print('score corr ux:', corr, len(labelsux), len(labelsux[rndmsk]),np.min(labelsux),np.min(u),np.max(labelsux),np.max(u))
                 v = mly.predict(features[rndmsk])
                 corr = np.corrcoef(labelsuy[rndmsk], v)[0, 1]
                 print('score corr uy:', corr)
                 w = mlz.predict(features[rndmsk])
                 corr = np.corrcoef(labelsuz[rndmsk], w)[0, 1]
                 print('score corr uz:', corr)
                 
                 
                 testmin = labelsux[rndmsk].min()
                 testmax = labelsux[rndmsk].max()
                 test11 = np.linspace(testmin, testmax, 100)
                 plt.figure()
                 plt.scatter(labelsux[rndmsk], u, s=1, c='b')
                 plt.scatter(test11, test11, s=1, c='r')
                 plt.xlabel('testx')
                 plt.ylabel('predictvalues from file')
                 plt.title('ux learn from a friend')
                 plt.show()

                 
                 testmin = labelsuy[rndmsk].min()
                 testmax = labelsuy[rndmsk].max()
                 test11 = np.linspace(testmin, testmax, 100)
                 plt.figure()
                 plt.scatter(labelsuy[rndmsk], v, s=1, c='b')
                 plt.scatter(test11, test11, s=1, c='r')
                 plt.xlabel('testy')
                 plt.ylabel('predictvalues from file')
                 plt.title('uy learn from a friend')
                 plt.show()

                 
                 testmin = labelsuz[rndmsk].min()
                 testmax = labelsuz[rndmsk].max()
                 test11 = np.linspace(testmin, testmax, 100)
                 plt.figure()
                 plt.scatter(labelsuz[rndmsk], w, s=1, c='b')
                 plt.scatter(test11, test11, s=1, c='r')
                 plt.xlabel('testz')
                 plt.ylabel('predictvalues from file')
                 plt.title('uz learn from a friend')
                 plt.show()                     
                 

        if 5==5: # plot score plots
            a=features1[features1[:,3]<500].copy()
            a=features1[features1[:,0]<30].copy()
            b=a[a[:,3]>0].copy()
            b0=labelsux1[features1[:,0]<30].copy()
            b1=b0[a[:,3]>0].copy()
            
            b=features1.copy()
            # b1=labelsux1.copy()
            u2ml = ml2.predict(b)
            uml = mlx.predict(b)
            vml = mly.predict(b)
            wml = mlz.predict(b)
            
            uvml=(uml*uml+vml*vml)**0.5
            labelsuxy1=(labelsux1*labelsux1+labelsuy1*labelsuy1)**.5
            labelsu2=(labelsux1*labelsux1+labelsuy1*labelsuy1)**.5
            labelsu2=(labelsux1*labelsux1+labelsuy1*labelsuy1+labelsuz1*labelsuz1)**.5
            plt.figure()
            # plt.scatter(labelsuxy1,uvml,s=0.01)
            plt.scatter(labelsu21[::20],u2ml[::20],s=0.009)
            plt.xlabel('RANS')
            plt.ylabel('ML')
            plt.title('tlv u2:'+str(stat(labelsu21,u2ml, kind='r')))

            speedlimit=14.5
            X, Y, Z = density_estimation(labelsu21[((u2ml<speedlimit) & (labelsu21<speedlimit) & (u2ml!=0.) & (labelsu21!=0.))],
                                         u2ml[((u2ml<speedlimit) & (labelsu21<speedlimit) & (u2ml!=0.) & (labelsu21!=0.))])
            # plt.figure()
            fig, ax = plt.subplots()   
            # Show density
            # ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
                      # extent=[min(labelsu21), max(labelsu21), min(u2ml), max(u2ml)]) 
            # Add contour lines
            CS = plt.contour(X, Y, Z, 500)
            plt.clabel(CS, inline=1, fontsize=10)
            plt.scatter(labelsu21[((u2ml<speedlimit) & (labelsu21<speedlimit) & (u2ml!=0.) & (labelsu21!=0.))],u2ml[((u2ml<speedlimit) & (labelsu21<speedlimit) & (u2ml!=0.) & (labelsu21!=0.))],s=0.1)
            plt.xlabel('observations')
            plt.ylabel('model')

            steps=np.zeros(20)
            bins =np.zeros(20)
            for i in range(len(steps)):
                steps[i]=np.abs(x-y)[(x>=0.25*i) & (x<0.25*(i+1))].std()
                bins[i] = 0.25*i+0.125
            plt.figure()
            plt.plot(bins,steps)


            plt.figure()
            plt.scatter(x,(y-x),s=0.01)

            plt.figure()
            plt.scatter(x,np.abs(x-y),s=0.01)
            X1, Y1, Z1 = density_estimation(x,np.abs(x-y))
            fig, ax = plt.subplots()   
            CS = plt.contour(X1, Y1, Z1, 50)
            plt.clabel(CS, inline=1, fontsize=10)
            plt.scatter(x, np.abs(x-y), s=0.1)
            plt.xlabel('observations')
            plt.ylabel('model')
        
        if 5==5:
            
            streetwide=30
            for i in range(streetwide):
                street = [
#                            3., #grid_z.ravel(),  #1
                            i, #dplim, #4
                            streetwide - i, #dnlim, #4
                            dllim.max(), # dllim, #4
                            drlim.max(), #drlim, #4
                            3, #grid_h.ravel(),  #2
                            1., #grid_h.ravel()/grid_z.ravel(),  #2
                            i, #grid_dp.ravel(), #4
                            i/3., #grid_dp.ravel()/grid_z.ravel(), #4
                            streetwide - i, #grid_dn.ravel(), #5
                            streetwide, #np.abs(grid_dp+grid_dn).ravel(),  #6
                            200, #grid_dl.ravel(),  #6
                            200, #grid_dr.ravel(),  #7
                            400, #np.abs(grid_dl+grid_dr).ravel(),  #7
                            20.,# grid_hp.ravel(),  #10
    #                        grid_hn.ravel(),  #11
                            1. #angle #12
                            ]
                
                print(i, ml2.predict([np.asarray(street)]))

            streetwide  = 20
            parcelheight = 5
            streetlength = 220
            for i in range(streetwide):
                street = [
#                            parcelheight, #grid_z.ravel(),  #1
                            dplim.max(), #dplim, #4
                            dnlim.max(), #dnlim, #4
                            i, # dllim, #4
                            streetwide - i, #drlim, #4
                            parcelheight, #grid_h.ravel(),  #2
                            1., #grid_h.ravel()/grid_z.ravel(),  #2
                            streetlength, #grid_dp.ravel(), #4
                            streetlength/parcelheight, #grid_dp.ravel()/grid_z.ravel(), #4
                            streetlength, #grid_dn.ravel(), #5
                            streetlength+streetlength, #np.abs(grid_dp+grid_dn).ravel(),  #6
                            i, #grid_dl.ravel(),  #6
                            streetwide-i, #grid_dr.ravel(),  #7
                            streetwide, #np.abs(grid_dl+grid_dr).ravel(),  #7
                            20.,# grid_hp.ravel(),  #10
    #                        grid_hn.ravel(),  #11
                            0.4 #angle #12
                            ]
                
                print(i, ml2.predict([np.asarray(street)]))


        if 5==5:
            print ('start saving ml results')
            casedir = textboxfile1value
            casedir = casedir[:casedir.rfind('/')+1]+str(int(textboxtime1value))+'/U'
            casedir = '/data4bk/nirb/Simulations/michaelstadtfloornoborders/10000/U'
            fhur = open(casedir, 'r')
#            fhur = open(r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor/2000/U', 'r')
            fhuw = open(r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor1ml/0/U', 'w')
            casedir = '/data4bk/nirb/Simulations/michaelstadtfloornobordersml/0/U'
            fhuw = open(r'/data4bk/nirb/Simulations/michaelstadtfloornobordersml/0/U', 'w')
            utxt = fhur.readlines()
            posu = utxt.index('internalField   nonuniform List<vector> \n')
            if full:
                fhpr = open(r'/data2/nirb/openFOAM/windAroundurbanMichelstadt2zomegascaleeps/bkup2000/p', 'r')
                fhkr = open(r'/data2/nirb/openFOAM/windAroundurbanMichelstadt2zomegascaleeps/bkup2000/k', 'r')
                fhor = open(r'/data2/nirb/openFOAM/windAroundurbanMichelstadt2zomegascaleeps/bkup2000/omega', 'r')
                fhnr = open(r'/data2/nirb/openFOAM/windAroundurbanMichelstadt2zomegascaleeps/bkup2000/nut', 'r')
                fhpw = open(r'/data2/nirb/openFOAM/windAroundurbanMichelstadt2zomegascaleepsml/0/p', 'w')
                fhkw = open(r'/data2/nirb/openFOAM/windAroundurbanMichelstadt2zomegascaleepsml/0/k', 'w')
                fhow = open(r'/data2/nirb/openFOAM/windAroundurbanMichelstadt2zomegascaleepsml/0/omega', 'w')
                fhnw = open(r'/data2/nirb/openFOAM/windAroundurbanMichelstadt2zomegascaleepsml/0/nut', 'w')
                ptxt = fhpr.readlines()
                ktxt = fhkr.readlines()
                otxt = fhor.readlines()
                ntxt = fhnr.readlines()
                posp = ptxt.index('internalField   nonuniform List<scalar> \n')
                posk = ktxt.index('internalField   nonuniform List<scalar> \n')
                poso = otxt.index('internalField   nonuniform List<scalar> \n')
                posn = ntxt.index('internalField   nonuniform List<scalar> \n')
            print('before items = ', posu)
            print('before items2 = ', utxt[posu+1])
            items = int(utxt[posu+1])
            print('items = ', items)
            predictfeatures=[]
            print('axis:',x.min(),(x[1]-x[0]),y.min(),(y[1]-y[0]),z.min(),(z[1]-z[0]))

            print('minmaxlabel', datetime.datetime.now(), min(labelsux), max(labelsux),min(gridu0.ravel()),max(gridu0.ravel()))
            
            ################3 nir start
            xg, yg ,zg = np.meshgrid(x, y, z, indexing='ij', sparse=True)
            from scipy.interpolate import RegularGridInterpolator
            pts1=np.asarray([xs,ys,zs]).T
            interpgrid_h = RegularGridInterpolator((x, y, z), grid_h, bounds_error=False,fill_value=0.)
            interpgrid_dp = RegularGridInterpolator((x, y, z), grid_dp, bounds_error=False,fill_value=0.)
            interpgrid_dn = RegularGridInterpolator((x, y, z), grid_dn, bounds_error=False,fill_value=0.)
            interpgrid_dpn = RegularGridInterpolator((x, y, z), np.abs(grid_dp+grid_dn), bounds_error=False,fill_value=0.)
            interpgrid_dl = RegularGridInterpolator((x, y, z), grid_dl, bounds_error=False,fill_value=0.)
            interpgrid_dr = RegularGridInterpolator((x, y, z), grid_dr, bounds_error=False,fill_value=0.)
            interpgrid_dlr = RegularGridInterpolator((x, y, z), np.abs(grid_dl+grid_dr), bounds_error=False,fill_value=0.)
            interpgrid_hp = RegularGridInterpolator((x, y, z), grid_hp, bounds_error=False,fill_value=0.)
            interpgrid_hr = RegularGridInterpolator((x, y, z), grid_hr, bounds_error=False,fill_value=0.)
            interpgrid_hl = RegularGridInterpolator((x, y, z), grid_hl, bounds_error=False,fill_value=0.)
            interpgrid_z = RegularGridInterpolator((x, y, z), grid_z, bounds_error=False,fill_value=0.)
            interpgrid_angle = RegularGridInterpolator((x, y, z), angle.reshape(grid_h.shape[0], grid_h.shape[1],grid_h.shape[2]), bounds_error=False,fill_value=0.)
            igrid_h = interpgrid_h(pts1)
            igrid_dp = interpgrid_dp(pts1)
            igrid_dn = interpgrid_dn(pts1)
            igrid_dpn = interpgrid_dpn(pts1)
            igrid_dl = interpgrid_dl(pts1)
            igrid_dr = interpgrid_dr(pts1)
            igrid_dlr = interpgrid_dlr(pts1)
            igrid_hp = interpgrid_hp(pts1)
            igrid_hr = interpgrid_hr(pts1)
            igrid_hl = interpgrid_hl(pts1)
            igrid_z = interpgrid_z(pts1)
            igrid_angle = interpgrid_angle(pts1)

            allfeatures = [
#                        grid_x.ravel(),  #0
#                        grid_y.ravel(),  #1
                        igrid_h,  #2
#                        grid_t.ravel(),  #3
                        igrid_dp, #4
                        igrid_dn, #5
                        igrid_dpn,  #6
#                        np.abs(grid_dp+grid_dn).ravel()/grid_hp.ravel(),  # W/H
                        igrid_dl,  #6
                        igrid_dr,  #7
                        igrid_dlr,  #7
                        igrid_hp,  #10
#                        grid_hn.ravel(),  #11
                        igrid_angle] #12
            
            lim=100.
            dplim=igrid_dp.ravel()
            dplim[dplim>dx*lim]=dx*lim
            dnlim=igrid_dn.ravel()
            dnlim[dnlim>dx*lim]=dx*lim
            dllim=igrid_dl.ravel()
            dllim[dllim>dy*lim]=dy*lim
            drlim=igrid_dr.ravel()
            drlim[drlim>dy*lim]=dy*lim
            dzlim=igrid_z.ravel()
            dzlim[dzlim>dz*lim]=dz*lim
            
            allfeatures = [
#                        grid_x.ravel(),  #0
#                        grid_y.ravel(),  #1
#                        grid_z.ravel(),      #1
                        dplim, #4
                        dnlim, #4
                        dllim, #4
                        drlim, #4
                        igrid_h.ravel(),  #2
                        igrid_dl.ravel()/igrid_hl.ravel(),  #2
                        igrid_dr.ravel()/igrid_hr.ravel(),  #3
                        igrid_dp.ravel(), #4
                        igrid_dp.ravel()/igrid_hp.ravel(), #4
                        igrid_dn.ravel(), #5
                        np.abs(igrid_dp+igrid_dn).ravel(),  #6
#                        np.abs(grid_dp+grid_dn).ravel()/grid_hp.ravel(),  # W/H
                        igrid_dl.ravel(),  #6
                        igrid_dl.ravel()/igrid_hl.ravel(),
                        igrid_dr.ravel(),  #7
                        np.abs(igrid_dl+igrid_dr).ravel(),  #7
                        igrid_hp.ravel(),  #10
                        igrid_hl.ravel(),  #10
                        igrid_hr.ravel(),  #10
#                        grid_hn.ravel(),  #11
                        igrid_angle] #12
            

            features = np.asarray(allfeatures)




            ################# nir end
            
            
            dx = (x[1]-x[0])
            dy = (y[1]-y[0])
            dz = (z[1]-z[0])
            print('testmin0', dx, dy, dz, x.min())
            print('testmin02', xs)
            xsi = np.round((xs - x.min())/dx)
            ysi = np.round((ys - y.min())/dy)
            zsi = np.round((zs - z.min())/dz)
            print('testmin1', ysi)
            xsi = xsi.astype(int)
            ysi = ysi.astype(int)
            zsi = zsi.astype(int)
            # pidx = np.zeros_like(xsi)
            print('testmin2', zsi)
            pidx = zsi+ysi*len(z)+xsi*len(z)*len(y)
            # pidx = np.zeros(len(x)*len(z)*len(y))
            print('testmin3')
            predictfeatures = features.T#[pidx]
            predictfeatures[np.isnan(predictfeatures)]=0.
            predictfeatures[np.isinf(predictfeatures)]=0.

            print('fin buildint test', datetime.datetime.now(), len(predictfeatures))
            labelspredictx = mlx.predict(predictfeatures)
            print('fin predict U_x', datetime.datetime.now())
            labelspredicty = mly.predict(predictfeatures)
            print('fin predict U_y', datetime.datetime.now())
            labelspredictz = mlz.predict(predictfeatures)
            print('fin predict U_z', datetime.datetime.now())
            if full:
                labelspredictp = mlp.predict(predictfeatures)
                print('fin predict p', datetime.datetime.now())
                labelspredictk = mlk.predict(predictfeatures)
                print('fin predict k', datetime.datetime.now())
                labelspredicto = mlo.predict(predictfeatures)
                print('fin predict o', datetime.datetime.now())
                labelspredictn = mln.predict(predictfeatures)
                print('fin predict n', datetime.datetime.now())
            # labelspredict[msklabel] = 0.0
            for i in range(len(labelspredictx)):  # todo items
                u0 = str(labelspredictx[i])
                u1 = str(labelspredicty[i])
                u2 = str(labelspredictz[i])
                utxt[posu+3+i]='('+u0+' '+u1+' '+u2+')\n'
                if full:
                    p = str(labelspredictp[i])
                    k = str(labelspredictk[i])
                    o = str(labelspredicto[i])
                    n = str(labelspredictn[i])
                    ptxt[posp + 3 + i] = p + '\n'
                    ktxt[posk + 3 + i] = k + '\n'
                    otxt[poso + 3 + i] = o + '\n'
                    ntxt[posn + 3 + i] = n + '\n'
            fhuw.writelines(utxt)
            fhuw.close()
            if full:
                fhpw.writelines(ptxt)
                fhkw.writelines(ktxt)
                fhow.writelines(otxt)
                fhnw.writelines(ntxt)
                fhpw.close()
                fhkw.close()
                fhow.close()
                fhnw.close()
            print('fin new U, p, k, omega, nut', datetime.datetime.now())

        if 5==5:       
            print('start test ml')
          
            grid_x = np.load(learnfile+'x.npy')
            grid_y = np.load(learnfile+'y.npy')
            grid_z1 = np.load(learnfile+'z.npy')  
#            print('grid_z', grid_z1)
            gridux = np.load(learnfile+'ux.npy')
            allfeatures = np.load(learnfile+'allfeatures.npy')
            features = np.load(learnfile+'features.npy')
            print('mlx3', learnfile, datetime.datetime.now())
            mlx = ml()  
            mlx.load('mlx-'+learnfile)
            indexprint = 5         
            z = np.unique(allfeatures[:,2])
            print ('grid_x2',grid_x.min(), grid_x.max())
            print ('mlx:  ', mlx)
#            print ('before print predict1', z[indexprint])
#            print ('before print predict11', grid_z.ravel())
#            print ('before print predict12', len(allfeatures))
            featureslevel = features[grid_z.ravel()==z[indexprint]]     
#            print ('before print predict2', featureslevel)
#            print ('before print predict21', z)
            labelspredictx = mlx.predict(featureslevel)
#            print ('after print predict')
            uxml = labelspredictx.reshape(gridux.shape[0],gridux.shape[1])
#            print ('after print predict', uxml)
            uxml[gridux[:,:,indexprint]==0]=0
#            print ('after predict reshape', uxml)
            plt.figure()
#            print ('after predict reshape11', uxml)
            plt.subplot(2,2,1) # , sharex=True, sharey=True
            print ('after predict reshape11')
#            print('a', gridux)
#            print('b', gridux[:,:,indexprint])
#            plt.imshow(uxml)
            plt.imshow(gridux[:,:,indexprint])
#            plt.title ('original')
            plt.colorbar()
#            print('sofarsogood1')
            plt.subplot(2,2,2)
            plt.imshow(uxml)
#            print('sofarsogood2')
            plt.title ('ml')
            plt.colorbar()
            plt.show()            
            plt.subplot(2,2,3)
            plt.imshow(gridux[:,:,indexprint] - uxml)
#            print('sofarsogood3')
            
            plt.title ('original - ml')
            plt.colorbar()
            plt.show()            
            plt.subplot(2,2,4)
            mainsimulation = self.experiment.currentText()
#            print('sofarsogood4')
            observation = chooseobservation(mainsimulation)  
#            print('sofarsogood5')
#            print('sofarsogood5',observation)
            showfieldtext = self.listfields.currentItem().text()    
            print ('con;t',len(observation))
            model, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation,
                                                  self.showfieldlist[showfieldtext])
            model = np.asarray(model)
            observation = np.asarray(observation)[:,0]
            print ('con;t2',len(model), len(features))
            dx = (x[1]-x[0])
            dy = (y[1]-y[0])
            dz = (z[1]-z[0])            
#            print('predict start4',x.min(), dx)
#            print('predict start5',y.min(), dy)
#            print('predict start6',z.min(), dz)
            for i in range(len(model)):
#                print('predict start')
#                print('predict start1',grid_x.ravel()[i])
#                print('predict start2',grid_y.ravel()[i])
#                print('predict start3',grid_z.ravel()[i])
#                indexx = np.where(x==grid_x.ravel()[i])[0][0]
#                print('predict startx', indexx)
#                indexy = np.where(y==grid_y.ravel()[i])[0][0]
#                print('predict starty', indexy)
#                indexz = np.where(z==grid_z.ravel()[i])[0][0]
#                print('predict startz', indexz)
#                pidx = indexz+indexy*len(z)+indexx*len(z)*len(y)
#                print(i,indexx, len(x),indexy, len(y),indexz, len(z), len(features))
#                print('features[pidx]:',features[pidx])
#                model[i]=mlx.predict(features[pidx])
#                print('predict;',model[i])
                
                xsi = round((coordinates[i,0] - x.min())/dx)
#                print('predict starta1',xsi)
                ysi = round((coordinates[i,1] - y.min())/dy)
#                print('predict starta1',ysi)
                zsi = round((coordinates[i,2] - z.min())/dz)
#                print('predict starta1',zsi)
                xsi = int(xsi)
#                print('predict starta2',xsi)
                ysi = int(ysi)
#                print('predict startb2',ysi)
                zsi = int(zsi)
#                print('predict startc2',zsi)
                # pidx = np.zeros_like(xsi)
                pidx = zsi+ysi*len(z)+xsi*len(z)*len(y)
#                print( '<<<<<<<<<<', coordinates[i,:])
#                print ('>>>>>>>>>>>>>>>>try', i, pidx, grid_x.ravel()[pidx], xsi, grid_y.ravel()[pidx], ysi, grid_z.ravel()[pidx],  zsi)
                
                
#            fea =  np.asarray([3.92786880e+00,    3.92786880e+00,   2.00e+00,   2.00000000e+00, 
#               4.00,   800.0e+00,   600.e+00,           8.52369517e-02])
#            modfea = mlx.predict(fea)
#            print ('>>',modfea)
               
            print ('con;t3',observation)
            print ('con;t4',model)
            plt.scatter(model, observation)
            print ('con;t5')
            plt.xlabel('ml')
            plt.ylabel('observation')
            corr = np.corrcoef(model, observation)[0, 1]
            plt.title ('obs vs XXXmlXX '+ str(corr) + ', '+ str(len(model)))
            plt.show()     

#            predictfeatures = features[1234]            
#            test2 = mlx.predict(predictfeatures)
#            print ('21<<',features[1234])
#            print ('2<<',allfeatures[1234])
#            print ('2>>',test2)
        print ('fin ml', datetime.datetime.now())

    def on_click_sinks(self):

        print('start sink', datetime.datetime.now())

        textboxfile1value = self.filename1.text()
        textboxtime1value = float(self.itteration1.text())

        sel = self.listfields.selectedIndexes()
        db = get_clip(textboxfile1value, [textboxtime1value],
                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=self.reader1)


        # ap3 = db[['x', 'y', 'z', sel[0].data()]]  # at the ground (1 meter)
        ap3 = db[['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)
        points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
        # values = ap3[sel[0].data()]
        valuesux = ap3['U_x']
        valuesuy = ap3['U_y']
        valuesuz = ap3['U_z']

        xmin = float(self.xmin.text())
        xmax = float(self.xmax.text())
        ymin = float(self.ymin.text())
        ymax = float(self.ymax.text())
        zmin = float(self.zmin.text())
        zmax = float(self.zmax.text())
        if xmin==xmax:
            xmin=ap3['x'].min()
            xmax=ap3['x'].max()
        if ymin==ymax:
            ymin=ap3['y'].min()
            ymax=ap3['y'].max()
        if zmin==zmax:
            zmin=ap3['z'].min()
            zmax=ap3['z'].max()
        # ticks = 600j
        grid_points = 10000000 # 7000000
        meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)

        resolution = (meters / grid_points)**(1./3)
        if resolution==0:
            print('no clip boundaries in one of the axis')
        # resolution = 0.1
        ticksx = int((xmax-xmin) / resolution) * 1j
        ticksy = int((ymax-ymin) / resolution) * 1j
        ticksz = int((zmax-zmin) / resolution) * 1j
        print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)

        xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
        # grid = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(values), (xi, yi, zi),
        #                    method='nearest')  # Nearest for keeping zeros that are buildings
        # gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
        #                    method='linear')  # Nearest for keeping zeros that are buildings
        # griduy = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuy), (xi, yi, zi),
        #                    method='linear')  # Nearest for keeping zeros that are buildings
        # griduz = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuz), (xi, yi, zi),
        #                    method='linear')  # Nearest for keeping zeros that are buildings

        def gd(result, x, y, z, v, xi, yi, zi, m):
            result[0] = griddata((x, y, z), v, (xi, yi, zi), method=m)  # Nearest for keeping zeros that are buildings

            return

        threads = []

        gridux = [np.zeros([int(ticksx.imag), int(ticksy.imag), int(ticksz.imag)]) + 1.]
        griduy = [np.zeros([int(ticksx.imag), int(ticksy.imag), int(ticksz.imag)])]
        griduz = [np.zeros([int(ticksx.imag), int(ticksy.imag), int(ticksz.imag)])]

        t = Thread(target=gd, args=(
        gridux, points3[:, 0], points3[:, 1], points3[:, 2], np.asarray(valuesux), xi, yi, zi, 'nearest'))
        threads.append(t)
        t = Thread(target=gd, args=(
        griduy, points3[:, 0], points3[:, 1], points3[:, 2], np.asarray(valuesuy), xi, yi, zi, 'nearest'))
        threads.append(t)
        t = Thread(target=gd, args=(
        griduz, points3[:, 0], points3[:, 1], points3[:, 2], np.asarray(valuesuz), xi, yi, zi, 'nearest'))
        threads.append(t)

        # Start all threads
        for x in threads:
            print('AAA', x, datetime.datetime.now())
            x.start()

        # Wait for all of them to finish
        for x in threads:
            print('BBB', x, datetime.datetime.now())
            x.join()

#        x = np.linspace(ap3['x'].min(), ap3['x'].max(), int(ticksx.imag))
#        y = np.linspace(ap3['y'].min(), ap3['y'].max(), int(ticksy.imag))
#        z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticksz.imag))

        rollxp = np.roll(gridux,1,axis =0)
        rollyp = np.roll(griduy,1,axis =1)
        rollzp = np.roll(griduz,1,axis =2)
        rollxm = np.roll(gridux,-1,axis =0)
        rollym = np.roll(griduy,-1,axis =1)
        rollzm = np.roll(griduz,-1,axis =2)
        rollxyz = rollxp + rollyp + rollzp - rollxm - rollym - rollzm
        gridxyz = gridux + griduy + griduz
        plt.figure()
        plt.subplot(1,2,1)
        plt.imshow(rollxyz[:,:,10])
        plt.colorbar()
        plt.subplot(1,2,2)
        plt.imshow(gridxyz[:,:,10])
        plt.colorbar()
        plt.show()




        print('fin sink', datetime.datetime.now())

    def on_click_statistics(self):
        gui=False
        rtheta=True
        mainsimulation = self.experiment.currentText()
        observation = chooseobservation(mainsimulation)
        
        if gui:
            files=[]        
            listItems=self.listfiles.selectedItems()
           
            if listItems:           
                for item in listItems:
                   files.append(item.text())
            else:
                for item in range(self.listfiles.count()):
                    files.append(self.listfiles.item(item).text())
    
            textboxfile1value = files[0]
        else:
            textboxfile1value = r'/data4bk/nirb/Simulations/michaelstadtfloor1ml/michaelstadtfloor1ml.foam'

        textboxtime1value = float(self.itteration1.text())
        
        items = self.listfields.selectedItems()
        fields = []
        for i in range(len(items)):
            fields.append(str(self.listfields.selectedItems()[i].text()))

        showfieldtext = self.listfields.currentItem().text()
        print('debugfield',showfieldtext , self.showfieldlist[showfieldtext])
        model, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation,
                                              self.showfieldlist[showfieldtext])
        
        if len(fields)>1:
            print('debugfield',showfieldtext , self.showfieldlist[showfieldtext])
            model2, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation,
                                                  1)
            

        if self.use2.checkState():
            textboxfile2value = self.filename2.text()
            textboxtime2value = float(self.itteration2.text())
            
            model2, coordinates2 = read_measurement(textboxfile2value, textboxtime2value, mainsimulation,
                                              self.showfieldlist[showfieldtext])
            rmseval = rmse(model, model2)
            print('rmse 2 files:', rmseval)        
        
        if float(self.xmin.text()) == float(self.xmax.text()):
            boolcoordx = (coordinates[:, 0] == coordinates[:, 0])
        else:
            boolcoordx = (coordinates[:, 0] > float(self.xmin.text())) & (coordinates[:, 0] < float(self.xmax.text()))
        if float(self.ymin.text()) == float(self.ymax.text()):
            boolcoordy = (coordinates[:, 1] == coordinates[:, 1])
        else:
            boolcoordy = (coordinates[:, 1] > float(self.ymin.text())) & (coordinates[:, 1] < float(self.ymax.text()))
        if float(self.zmin.text()) == float(self.zmax.text()):
            boolcoordz = (coordinates[:, 2] == coordinates[:, 2])
        else:
            boolcoordz = (coordinates[:, 2] > float(self.zmin.text())) & (coordinates[:, 2] < float(self.zmax.text()))
        boolcoord = boolcoordx & boolcoordy & boolcoordz
        newcoord = coordinates[boolcoord]
        print('lennewcoord', len(newcoord), len(coordinates))

        if mainsimulation == 'observationB':
            observation = observation[:, self.showfieldlist[showfieldtext]]
            print('mainsimulation-observationB:', mainsimulation)
        elif mainsimulation == 'Michelstadt':
            if len(fields)>1:
                observation = (observation[:,0]**2.+observation[:,1]**2.)**.5
                # observation[:,1][observation[:,1]==0]=0.0001
                # observationbkup=observation.copy()
                # observation = np.sin(observation[:,0]/observation[:,1])              
            else:
                observation = observation[:, self.showfieldlist[showfieldtext]]
                
            print('mainsimulation-Michel:', mainsimulation)

        observation = observation[boolcoord]
        if len(fields)>1:
            model= (model[boolcoord]**2.+model2[boolcoord]**2.)**.5
            # model2[boolcoord][model2[boolcoord]==0.0]=0.0001
            # model= np.sin(model[boolcoord]/model2[boolcoord])
        else:
            model = model[boolcoord]

        if 5==6:
            lowwinds=observationbkup[:,0]>.5
            observation=observation[lowwinds]
            model=model[lowwinds]
            print('cut of by wind gives len of ', len(observation))


        fac2 = 0
        for i in range(len(observation)):
            if abs(observation[i] / 2) <= abs(model[i]) <= abs(observation[i] * 2):
                fac2 = fac2 + 1
        fac2 = float(fac2) / float(len(model))

        print(textboxfile1value)
        print(textboxtime1value)
        print(showfieldtext)
        print('newcoord',
              float(self.xmin.text()), float(self.xmax.text()), float(self.ymin.text()), float(self.ymax.text()))
        print('observations=' + str(np.mean(observation))+'+-'+str(np.std(observation)))
        print('simulation  =' + str(np.mean(model))+'+-'+str(np.std(model)))
        print('fac2=' + str(fac2))
        print('median_absolute_error=' + str(sklearn.metrics.median_absolute_error(observation, model)))
        print('mae=' + str(sklearn.metrics.mean_absolute_error(observation, model)))
        print('rmse=' + str(rmse(observation, model)))
        print('std=' + str(np.std(observation-model)))
        print('r=' + str(np.corrcoef(observation, model)[0][1]))
        print('explained_variance_score=' + str(sklearn.metrics.explained_variance_score(observation, model)))
        print('r2_score=' + str(sklearn.metrics.r2_score(observation, model)))
        print('Frictional Bias=' + str((np.mean(observation) - np.mean(model))/(np.mean(observation) + np.mean(model))/0.5))
        fbindex = model > observation
        print('Frictional Bias fp=' + str((np.abs(np.mean(observation[fbindex] - model[fbindex]))+(np.mean(model[fbindex]) - np.mean(observation[fbindex])))/(np.mean(observation[fbindex]) + np.mean(model[fbindex]))))  # nopep8
        print('Frictional Bias fn=' + str((np.abs(np.mean(observation[~fbindex] - model[~fbindex]))+(np.mean(observation[~fbindex]) - np.mean(model[~fbindex])))/(np.mean(observation[~fbindex]) + np.mean(model[~fbindex]))))  # nopep8
        
        plt.figure()
        plt.scatter(observation,model)
        plt.xlabel('wind tunnel [m/s]')
        plt.ylabel('model [m/s]')
        
        

        
    def on_click_experiment(self):
        gui=False
        cca = True
        mainsimulation = self.experiment.currentText()
        observation = chooseobservation(mainsimulation)

        if gui:
            files=[]        
            listItems=self.listfiles.selectedItems()
           
            if listItems:           
                for item in listItems:
                   files.append(item.text())
            else:
                for item in range(self.listfiles.count()):
                    files.append(self.listfiles.item(item).text())
    
            textboxfile1value = files[0]
        else:
            textboxfile1value = r'/data4bk/nirb/Simulations/michaelstadtfloor1ml/michaelstadtfloor1ml.foam'
        textboxtime1value = float(self.itteration1.text())

        showfieldtext = self.listfields.currentItem().text()

        if cca:
            modelx, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation, 0)
            modely, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation, 1)
            modelz, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation, 2)

        model, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation,
                                              self.showfieldlist[showfieldtext])
           
        
#        probefile = 'OFTutorials/'+traj.v_root[i].z+'/postProcessing/probes/0/U'

        # if mainsimulation == 'observationE':
        #     model *= 3.12


        if mainsimulation == 'observationB':
            observation = observation[:, self.showfieldlist[showfieldtext]]
        elif mainsimulation == 'Michelstadt':
            if cca:
                observationx = observation[:, 0]
                observationy = observation[:, 1]
            observation = observation[:, self.showfieldlist[showfieldtext]]

        if float(self.xmin.text()) == float(self.xmax.text()):
            boolcoordx = (coordinates[:, 0] == coordinates[:, 0])
        else:
            boolcoordx = (coordinates[:, 0] > float(self.xmin.text())) & (coordinates[:, 0] < float(self.xmax.text()))
        if float(self.ymin.text()) == float(self.ymax.text()):
            boolcoordy = (coordinates[:, 1] == coordinates[:, 1])
        else:
            boolcoordy = (coordinates[:, 1] > float(self.ymin.text())) & (coordinates[:, 1] < float(self.ymax.text()))
        if float(self.zmin.text()) == float(self.zmax.text()):
            boolcoordz = (coordinates[:, 2] == coordinates[:, 2])
        else:
            boolcoordz = (coordinates[:, 2] > float(self.zmin.text())) & (coordinates[:, 2] < float(self.zmax.text()))
        boolcoord = boolcoordx & boolcoordy & boolcoordz
        newcoord = coordinates[boolcoord]
        print('lennewcoord', len(newcoord), len(coordinates))

        print('debug', type(observation), type(boolcoord), len(observation), len(boolcoord[1:]), len(model), len(newcoord))
        if cca:
            observationx = observationx[boolcoord]
            modelx = modelx[boolcoord]
            observationy = observationy[boolcoord]
            modely = modely[boolcoord]
            modelz = modelz[boolcoord]
        observation = observation[boolcoord]
        model = model[boolcoord]
        print('debug3', type(observation), type(boolcoord), len(observation), len(boolcoord[1:]), len(model))

        fac2 = 0
        for i in range(len(observation)):
            if abs(model[i]) >= abs(observation[i]/2) and abs(model[i]) <= abs(observation[i]*2):
                fac2 = fac2+1
        fac2 = float(fac2) / float(len(model))

        X, Y, Z = density_estimation(observation, model)

        fig, ax = plt.subplots()

        # Show density
        ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
                  extent=[min(observation), max(observation), min(model), max(model)])

        # Add contour lines
        CS = plt.contour(X, Y, Z, 10)
        plt.clabel(CS, inline=1, fontsize=10)
        # plt.scatter(observation,model)
        plt.xlabel('observations')
        plt.ylabel('model')

        plt.title('r=' + str(np.corrcoef(observation, model)[0][1]) + ', fac2=' + str(fac2)
                  + ',rmse=' + str(rmse(observation, model)))

        line11 = np.linspace(min(observation), max(observation))
        ax.plot(line11, line11, 'b-', markersize=2)
        ax.plot(observation, model, 'k.', markersize=2)
        plt.show()

    
        plt.figure()
        plt.scatter(observation,model,s=5)  # 2 for michaelstadt, 5 for case E
        plt.plot(model,model)       
        plt.xlabel('observation')
        plt.ylabel('model')
        plt.xlim(min(min(observation), min(model)), max(max(observation), max(model)))
        plt.ylim(min(min(observation), min(model)), max(max(observation), max(model)))
        plt.title('')
        plt.show()              
        
        plt.figure()
        n1, bins1, patches1 = plt.hist(model, 50, alpha=0.5, label='model')
        n2, bins2, patches2 = plt.hist(observation, 50, alpha=0.5, label='observations')
        # plt.hist(observation-model,alpha=0.5, label='observation-model')
        plt.legend(loc='upper right')
        # plt.plot(bins1, n1, '--')
        plt.show()

        plt.figure()
        plt.scatter(observation,model-observation)
        plt.xlabel('observation')
        plt.ylabel('model-observation')
        plt.title('error distribution')
        plt.show()
        if gui:
            textboxfile1value = self.filename1.text()
        textboxtime1value = float(self.itteration1.text())
        axisindex = int(self.listaxis.currentRow())
        axispos = float(self.axispos.text())
        if self.heightcontour.checkState():
            db = get_slice_height(textboxfile1value, textboxtime1value, axisindex, axispos,
                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                  reader=self.reader1)
        else:
            db = get_slice(textboxfile1value, textboxtime1value, axisindex, axispos,
                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()))
        cma = plt.cm.OrRd
        cma.set_bad(alpha=0.0)

        if len(db) == 0:
            print ('len(db) is 0 !!! ')
        else:
            plt.figure()
            if self.listaxis.currentRow() == 0:
                plt.tricontourf(db['y'], db['z'], db[self.listfields.currentItem().text()], 100)
            if self.listaxis.currentRow() == 1:
                plt.tricontourf(db['x'], db['z'], db[self.listfields.currentItem().text()], 100)
            if self.listaxis.currentRow() == 2:
                plt.tricontourf(db['x'], db['y'], db[self.listfields.currentItem().text()], 100)  # , cmap=cma

            plt.colorbar()
            # plt.scatter(newcoord[:,0],newcoord[:,1],c=observation)

            print('debug4', type(observation), type(boolcoord), len(observation), len(boolcoord[1:]), len(model), len(newcoord))

            if self.listaxis.currentRow() == 0:
                plt.scatter(newcoord[:, 1], newcoord[:, 2], c=observation, s=1)
            if self.listaxis.currentRow() == 1:
                plt.scatter(newcoord[:, 0], newcoord[:, 2], c=observation, s=1)
            if self.listaxis.currentRow() == 2:
                plt.scatter(newcoord[:, 0], newcoord[:, 1], c=observation, s=1)
            plt.title(textboxfile1value+'-'+self.listfields.currentItem().text())
            plt.show()
            
            

            # X = [[0., 0., 1.], [1., 0., 0.], [2., 2., 2.], [3., 5., 4.]]
            # Y = [[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]]

            if self.figsave.checkState():
                plt.savefig('/home/nirb/test.png')
                print('savefig')

            print('max observations:', np.max(observation))
            if cca:
                X = np.vstack((modely,modelx, modelz))
                Y = np.vstack((observationx,observationy))
                cca = CCA(n_components=1)
                cca.fit(X.T, Y.T)
                U_c, V_c = cca.fit_transform(X.T, Y.T)
                result = np.corrcoef(U_c.T, V_c.T)[0, 1]
                resultx = np.corrcoef(modelx, observationx)[0, 1]
                resulty = np.corrcoef(modely, observationy)[0, 1]
                print('cca:', result, resultx, resulty)

                
            plt.figure()
            plt.plot(model[115:127],coordinates[115:127][:,2],'r', linewidth=3, label='before building - model')
            plt.plot(observation[115:127],coordinates[115:127][:,2],'--r', linewidth=3, label='before building - observation')
            plt.xlabel('wind [m/s]')
            plt.ylabel('H [m]')
            plt.plot(model[159:171],coordinates[159:171][:,2],'b', linewidth=3, label='after building - model')
            plt.plot(observation[159:171],coordinates[159:171][:,2],'--b', linewidth=3, label='after building - observation')
            plt.legend(loc=2)
            plt.show()
            
    def on_click_porous(self):
        angle = 0 # 0 is west wind
        print ('!!!! we used angle:', angle, '!!!!')
        if self.learnfrom.text()!="":
            resolution = float(self.learnfrom.text())
        else:
            resolution = 150
        textboxfile1value = r'/ibdata2/nirb/openFOAM/angle/kleombenchmark/cylinder.foam'
        textboxfile1value = r'/ibdata2/nirb/openFOAM/porous/hadasfine/windAroundCube.foam'
        textboxfile1value = r'/ibdata2/nirb/openFOAM/porous/lambda08s2/windAroundCube.foam'
        print('start porous', datetime.datetime.now())
#        textboxfile1value = self.filename1.text()
        textboxtime1value = float(self.itteration1.text())
        
#        textboxfile1value = r'/ibdata2/nirb/openFOAM/ml/MALA2a120/windAroundCube.foam'
#        textboxtime1value = 1000.0

        try:
            print ('try')
            grid_z1n = np.load('grid_z1.npy')
            xi = np.load('xi.npy')
            yi = np.load('yi.npy')
            zi = np.load('zi.npy')
        except IOError:
            print ('except')
            if self.heightcontour.checkState():
                db = get_slice_height(textboxfile1value, textboxtime1value, 1, float(self.axispos.text()),
                                      offset=range(0, 200, 1),
                                      clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                      clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                      clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                      reader=self.reader1)
            else:
                db = get_clip(textboxfile1value, [textboxtime1value],
                               clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                               clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                               clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()))
#                db = get_clip(textboxfile1value, [textboxtime1value],
#                               clipxmin=0., clipxmax=0.,
#                               clipymin=0., clipymax=0.,
#                               clipzmin=0., clipzmax=0.)
    
            ticksh = 1000j  # 1000
            ticksv = 750j  # 600
    
            ap3 = db[['x', 'y', 'z', 'U_x']]  # at the ground (1 meter)
            points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
            values = ap3['U_x']
            print('x:',ap3['x'].min(),ap3['x'].max(),(ap3['x'].max()-ap3['x'].min())/ticksh)
            print('y:',ap3['y'].min(),ap3['y'].max(),(ap3['y'].max()-ap3['y'].min())/ticksh)
            print('z:',ap3['z'].min(),ap3['z'].max(),(ap3['z'].max()-ap3['z'].min())/ticksv)
            xi, yi, zi = np.mgrid[ap3['x'].min():ap3['x'].max():ticksh, ap3['y'].min():ap3['y'].max():ticksh,
                         ap3['z'].min():ap3['z'].max():ticksv]
            print('start griddata cubic', datetime.datetime.now())
            grid_z1n = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(values), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            print('end griddata cubic', datetime.datetime.now())

            np.save('xi',xi)           
            np.save('yi',yi)           
            np.save('zi',zi)           
            np.save('grid_z1',grid_z1n)
            
#        print ('return')
        
        grid_z1 = grid_z1n.copy()
        grid_z1[grid_z1n==0]=0

        print ('lc start')
        print(xi.shape)
        xi2 = xi[:,0,0]
        xi3 = np.zeros(int((xi[-1,0,0]-xi[0,0,0])/resolution))
        for i in range(len(xi3)):
            i1=int((float(len(xi2))/len(xi3)/2)+i*(float(len(xi2))/len(xi3)))
#            print(i,i1,xi2[i1])
            xi3[i]=xi2[i1]
        yi2 = yi[0,:,0]
        yi3 = np.zeros(int((yi[0,-1,0]-yi[0,0,0])/resolution))
        for i in range(len(yi3)):
            i1=int((float(len(yi2))/len(yi3)/2)+i*(float(len(yi2))/len(yi3)))
#            print(i,i1,yi2[i1])
            yi3[i]=yi2[i1]
        zi3 = zi[0,0,:]

        print('dx, dy, dz:',xi2[1]-xi2[0],yi2[1]-yi2[0],zi[0,0,1]-zi[0,0,0])
        grid_z2 = np.zeros([len(xi3), len(yi3), xi.shape[2]])
        
        testi=11  # 11 is 850
        testj=17  # 17 is 1350
        testi=19  # 17 is 1300
        testj=11  # 15 is 1150
        testi=40  # 17 is 1300
        testj=60  # 15 is 1150        
        k= 80  #80

        i1=testi
        j1=testj
        
        grid_z1[:,:,:][grid_z1[:,:,:]!=0.0]=1.
#        grid_z1[i,:,:][grid_z1[i,:,:]!=0.0]=1.

#        playbox=np.zeros([len(xi3),len(yi2)])

        di=float(len(xi2))/len(xi3)/2
        dj=float(len(yi2))/len(yi3)/2
        
        print('shapes:',grid_z1.shape, len(xi2),len(xi3))
        for k1 in range(len(zi3)):
            for i1 in range(len(xi3)):
    #        for i1 in range(testi,testi+1):
                print('i1 out of:',i1,len(xi3))
                for j1 in range(len(yi3)):
    #            for j1 in range(testj,testj):
    #            print(i1,j1,di,dj,di*(2*i1),di*(2*i1+2)-1,dj*(2*j1),dj*(2*j1+2)-1,k1)
                    playbox = grid_z1[int(di*(2*i1)):int(di*(2*i1+2)-1),int(dj*(2*j1)):int(dj*(2*j1+2)-1),k1]
                    nebe=0
                    for i in range(playbox.shape[0]):
                        for j in range(playbox.shape[1]):
                            if (playbox[i,j]==0):
                                if i==0:
                                    nebe+=1
                                elif (playbox[i-1,j]!=0):
                                    nebe+=1
                    nebe/=playbox.shape[0]
                    space=playbox.shape[0]*playbox.shape[1]
                    area=space-playbox.sum()  # because buildings are 0 and outdoor is 0
                    print('test1',area)
                    if area == 0:
                        lc1 = 0
                        lc = 888.
                    elif area == space:
                        lc1 = 999.
                        lc = 0
                    else:
                        lc1 = (nebe/((space-area)*(area/space)))*(yi2[2]-yi2[1])
                        lc=1./lc1
        #                else:
        #                    if np.isfinite(lc):
        #                        lc1=0
        #                    else:
        #                        lc1 = 99999999.
                    print('lc=',i1,j1,k1,nebe, space, area, lc, lc1, 'NEED TO MUL/DEVIDE by DX and not only dy')
                    grid_z2[i1,j1,k1]=lc1
        
        np.save('grid_z2',grid_z2)

        print ('lc end')
        maxlc = grid_z2.max()
        maxlc1 = grid_z2[grid_z2!=maxlc].max()
        print('maxlc:',maxlc,maxlc1)
        grid_z2[grid_z2==maxlc]=maxlc1
        
#        plt.figure()
##        plt.imshow(grid_z1[:,:,z])
#        plt.imshow(grid_z1[i,:,:])
#        plt.colorbar()
#        plt.show()
        
        print ('k=',k)
        plt.figure()
        plt.imshow(grid_z1[:,:,k])
        plt.colorbar()
        plt.show()
        
        plt.figure()
#        plt.imshow(playbox)
        plt.imshow(grid_z2[:,:,k])
        plt.colorbar()
        plt.show()
                     
        np.savetxt('porousx',xi3, fmt='%s')
        file1 = open("porousxi","w")
        print('start ')
        lines=[str(len(xi3))+" 1\n","(\n"]
        print(str(len(xi3))+" 1")
        for i in range(len(xi3)):
#                 print("("+xi3[i]+")")
                 lines.append("("+str(xi3[i])+")\n")
        print(")")
        lines.append(")\n")
        file1.writelines(lines)
        file1.close() #to change file access modes 
        print("VVVVVVVVVVVVVVVVVVV")
        np.savetxt('porousy',yi3, fmt='%s')
        file1 = open("porousyi","w")
        lines=[str(len(yi3))+" 1\n","(\n"]
        for i in range(len(yi3)):
                 lines.append("("+str(yi3[i])+")\n")
        lines.append(")\n")
        file1.writelines(lines)
        file1.close() #to change file access modes         
        
        np.savetxt('porousz',zi3, fmt='%s')
        file1 = open("porouszi","w")
        lines=[str(len(zi3))+" 1\n","(\n"]
        for i in range(len(zi3)):
                 lines.append("("+str(zi3[i])+")\n")
        lines.append(")\n")
        file1.writelines(lines)
        file1.close() #to change file access modes         
        
#        for i in range(len(zi)):
#            np.savetxt('porouslevel'+str(i).zfill(4),grid_z2[:,:,i], fmt='%s')
                
        file1 = open("porousleveli","w")
        lines=[str(len(xi3)*len(yi3)*len(zi3))+" 4\n","(\n"]
        for i in range(len(xi3)):
            for j in range(len(yi3)):
                    for k in range(len(zi3)):
                        lines.append("("+str(xi3[i])+" "+str(yi3[j])+" "+str(zi3[k])+" "+str(grid_z2[i,j,k])+")\n")
        lines.append(")\n")
        file1.writelines(lines)
        file1.close() #to change file access modes         


       # stacy clearing start
        file1 = open("porouslevellambda08s2","w")
        lines=[str(480*40*160)+" 4\n","(\n"]
        for i in range(480): # x
           for j in range(40): # y
                 for k in range(160): # z
                       if (k<12.1 and (i<200 or i>1000)): 
                           lines.append("("+str(i/1.)+" "+str(j/1.)+" "+str(k/1.)+" "+str(1/9.)+")\n")
                       else:
                           lines.append("("+str(i/1.)+" "+str(j/1.)+" "+str(k/1.)+" "+str(0.)+")\n")
        lines.append(")\n")
        file1.writelines(lines)
        file1.close() #to change file access modes         
       # stacy clearing end 
          
        ####################
        ## to add a loop on the playbox in order to do it on all particles
        ####################           
        
    def on_click_lsm(self):
        textboxfile1value = self.filename1.text()
        textboxtime1value = float(self.itteration1.text())

        db = get_clip(textboxfile1value, [textboxtime1value],
                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=self.reader1)


        # ap3 = db[['x', 'y', 'z', sel[0].data()]]  # at the ground (1 meter)
        ap3 = db[['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)
        points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
        # values = ap3[sel[0].data()]
        valuesux = ap3['U_x']
        valuesuy = ap3['U_y']
        valuesuz = ap3['U_z']

        xmin = float(self.xmin.text())
        xmax = float(self.xmax.text())
        ymin = float(self.ymin.text())
        ymax = float(self.ymax.text())
        zmin = float(self.zmin.text())
        zmax = float(self.zmax.text())
        if xmin==xmax:
            xmin=ap3['x'].min()
            xmax=ap3['x'].max()
        if ymin==ymax:
            ymin=ap3['y'].min()
            ymax=ap3['y'].max()
        if zmin==zmax:
            zmin=ap3['z'].min()
            zmax=ap3['z'].max()
        # ticks = 600j
        grid_points = 17000  # add 00
        meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)

        resolution = (meters / grid_points)**(1./3)
        if resolution==0:
            print('no clip boundaries in one of the axis')
        # resolution = 0.1
        ticksx = int((xmax-xmin) / resolution) * 1j
        ticksy = int((ymax-ymin) / resolution) * 1j
        ticksz = int((zmax-zmin) / resolution) * 1j
        print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)
        print(xmin,xmax,ticksx, ymin,ymax,ticksy, zmin,zmax,ticksz)

        print(np.mgrid[xmin:xmax:ticksx])

        xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
        print('ofwinda')

        gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        griduy = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuy), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        griduz = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuz), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        print('ofwindb')

        x = np.linspace(ap3['x'].min(), ap3['x'].max(), int(ticksx.imag))
        y = np.linspace(ap3['y'].min(), ap3['y'].max(), int(ticksy.imag))
        z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticksz.imag))
        print('ofwindc')

#        print(valuesuy.shape, datetime.datetime.now())
#        print('x',x)
#        with open("xcoord.txt", "w") as f:
#            for s in x:
#                f.write(str(s) +"\n")
#        with open("ycoord.txt", "w") as f:
#            for s in y:
#                f.write(str(s) +"\n")
#        with open("zcoord.txt", "w") as f:
#            for s in z:
#                f.write(str(s) +"\n")
#        
#        with open("u.txt", "w") as f:
#            for s in gridux:
#                f.write(str(s) +"\n")
#        with open("v.txt", "w") as f:
#            for s in griduy:
#                f.write(str(s) +"\n")
#        with open("w.txt", "w") as f:
#            for s in griduz:
#                f.write(str(s) +"\n")
                
###########################################################################
        
        from scipy.io.netcdf import netcdf_file
        
        # https://stackoverflow.com/questions/35609684/how-to-write-variable-attributes-in-netcdf-file-using-python
        print('ofwind0')
        
        # Define some dummy data
        time_arr = range(1)
        x_arr = x #np.array([30.5, 40., 40.5, 50])
        y_arr = y #np.array([200., 203., 206.])
        z_arr = z #np.array([20., 30., 60.])
        ntim, nx, ny, nz = len(time_arr), len(x_arr), len(y_arr), len(z_arr)
        print('ofwind1')
        
#        u_arr = np.random.randn(ntim, nx, ny, nz)
#        v_arr = np.random.randn(ntim, nx, ny, nz)
#        w_arr = np.random.randn(ntim, nx, ny, nz)
        u_arr = np.zeros([ntim, nx, ny, nz])
        v_arr = np.zeros([ntim, nx, ny, nz])
        w_arr = np.zeros([ntim, nx, ny, nz])
        print('ofwind2')

        u_arr[0,:,:,:] = gridux
        v_arr[0,:,:,:] = griduy
        w_arr[0,:,:,:] = griduz
        print('ofwind3')
        
        # Write out data to a new netCDF file with some attributes
        filename = netcdf_file('./ofwind.nc', 'w')
        
        # Dimensions
        filename.createDimension('time', ntim)
        filename.createDimension('x', nx)
        filename.createDimension('y', ny)
        filename.createDimension('z', nz)
        
        # Variables
        time = filename.createVariable('time', 'i', ('time',))
        x = filename.createVariable('x', 'f4', ('x',))
        y = filename.createVariable('y', 'f4', ('y',))
        z = filename.createVariable('z', 'f4', ('z',))
        u = filename.createVariable('u', 'f4', ('time', 'x', 'y', 'z',))
        v = filename.createVariable('v', 'f4', ('time', 'x', 'y', 'z',))
        w = filename.createVariable('w', 'f4', ('time', 'x', 'y', 'z',))
        
        # Attributes
        time.units = ''
        x.units = 'utm'
        y.units = 'utm'
        z.units = 'utm'
        u.units = 'm/s'
        u.missing_val = 1e20
        v.units = 'm/s'
        v.missing_val = 1e20
        w.units = 'm/s'
        w.missing_val = 1e20
        
        # Populate the variables with data
        time[:] = time_arr
        x[:] = x_arr
        y[:] = y_arr
        z[:] = z_arr
        u[:,:,:] = u_arr[:,:,:]
        v[:,:,:] = v_arr[:,:,:]
        w[:,:,:] = w_arr[:,:,:]
        
        filename.close()  
        print ('fin netcdf')
                                
    
    def on_click_error_classification(self):

        mainsimulation = self.experiment.currentText()
        observation = chooseobservation(mainsimulation)

        textboxfile1value = self.filename1.text()
        textboxtime1value = float(self.itteration1.text())

        showfieldtext = self.listfields.currentItem().text()

        model, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation,
                                              self.showfieldlist[showfieldtext])

        if mainsimulation == 'observationB':
            observation = observation[:, self.showfieldlist[showfieldtext]]
        elif mainsimulation == 'Michelstadt':
            observation = observation[:, self.showfieldlist[showfieldtext]]

        if float(self.xmin.text()) == float(self.xmax.text()):
            boolcoordx = (coordinates[:, 0] == coordinates[:, 0])
        else:
            boolcoordx = (coordinates[:, 0] > float(self.xmin.text())) & (coordinates[:, 0] < float(self.xmax.text()))
        if float(self.ymin.text()) == float(self.ymax.text()):
            boolcoordy = (coordinates[:, 1] == coordinates[:, 1])
        else:
            boolcoordy = (coordinates[:, 1] > float(self.ymin.text())) & (coordinates[:, 1] < float(self.ymax.text()))
        if float(self.zmin.text()) == float(self.zmax.text()):
            boolcoordz = (coordinates[:, 2] == coordinates[:, 2])
        else:
            boolcoordz = (coordinates[:, 2] > float(self.zmin.text())) & (coordinates[:, 2] < float(self.zmax.text()))
        boolcoord = boolcoordx & boolcoordy & boolcoordz
        newcoord = coordinates[boolcoord]
        print('lennewcoord', len(newcoord), len(coordinates))

        print('observations', len(observation), type(observation), observation)
        print('boolcoords', len(boolcoord), boolcoord, boolcoord)
        print('model', len(model), type(model), model)
        observation = observation[boolcoord]
        model = model[boolcoord]

        db = get_clip(textboxfile1value, [textboxtime1value],
                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=self.reader1)
        
        z80 = (max(db['z']) - db['z'].min()) * 0.8 + db['z'].min()
        print('z80', z80)

        # textboxfile1value = self.filename1.text()
        # textboxtime1value = float(self.itteration1.text())
        # axisindex = int(self.listaxis.currentRow())
        # axispos = float(self.axispos.text())
        # if self.heightcontour.checkState():
        #     db = get_slice_height(textboxfile1value, textboxtime1value, axisindex, axispos
        #               , clipxmin = float(self.xmin.text()), clipxmax = float(self.xmax.text())
        #               , clipymin = float(self.ymin.text()), clipymax = float(self.ymax.text())
        #               , clipzmin = float(self.zmin.text()), clipzmax = float(self.zmax.text()))
        # else:
        db = get_slice(textboxfile1value, textboxtime1value, 2, z80,
                       clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                       clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                       clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=self.reader1)
        ticks = 200
        xmin, xmax, ymin, ymax, zi80 = makegrid(db, 'U_x', ticks=ticks)
        # db = None

        ug = []
        ux = []
        x = []
        y = []
        z = []
        uy = []
        uxuy = []
        uz = []
        uxuz = []
        distwall = np.zeros(len(newcoord))
        distwall -= 999.0
        probe_location1 = pvsimple.ProbeLocation(Input=self.reader1,
                                                 ProbeType='Fixed Radius Point Source')
        # Properties modified on probe_location1
        probe_location1.Tolerance = 2.22044604925031e-16
        probe_location1.ProbeType.Radius = 1  # 21000.0
        print('len newcoord', len(newcoord))
        for i in range(len(newcoord)):  # len(newcoord)
            xx = newcoord[i][0]
            yy = newcoord[i][1]
            zx = int(((xx-xmin)/(xmax-xmin))*ticks)
            zy = int(((yy-ymin)/(ymax-ymin))*ticks)
            ug.append(zi80[zx][zy])
            # Properties modified on probe_location1.ProbeType
            probe_location1.ProbeType.Center = newcoord[i]  # [100000.0, 100000.0, 5000.0]
            poly_data = servermanager.Fetch(probe_location1)
            point_data = poly_data.GetPointData()
            rho_array = point_data.GetArray('U')
            value = rho_array.GetTuple(0)
            ux.append(value[0])
            uy.append(value[1])
            uz.append(value[2])
            if value[0] == 0:
                uxuy.append(0)
                uxuz.append(0)
            else:
                ux_uy_tmp = math.atan(value[1]/value[0])
                ux_uz_tmp = math.atan(value[2]/value[0])
                if value[0] < 0:
                    if value[1] < 0:
                        ux_uy_tmp -= math.pi
                    else:
                        ux_uy_tmp += math.pi
                    if value[2] < 0:
                        ux_uz_tmp -= math.pi
                    else:
                        ux_uz_tmp += math.pi
                uxuy.append(ux_uy_tmp)
                uxuz.append(ux_uz_tmp)

            # rho_array = point_data.GetArray('z')
            # value = rho_array.GetTuple(0)
            x.append(newcoord[i][0])
            y.append(newcoord[i][1])
            z.append(newcoord[i][2])
            if distwall[i] == -999.0:
                tickswall = 500
                Awall = get_slice(textboxfile1value, textboxtime1value, 2, newcoord[i][2],
                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                  reader=self.reader1)
                xminwall, xmaxwall, yminwall, ymaxwall, zwall = makegrid(Awall, 'U_x', ticks=tickswall)
                tickx = (xmaxwall-xminwall)/tickswall
                ticky = (ymaxwall-yminwall)/tickswall
                heightindex = np.where(newcoord[:,2]==z[-1])
                for hi in heightindex[0]:
                    currx = int(((newcoord[hi][0]-xminwall)/(xmaxwall-xminwall))*tickswall)
                    curry = int(((newcoord[hi][1]-yminwall)/(ymaxwall-yminwall))*tickswall)
                    mindist = min(xmaxwall-xminwall, ymaxwall-yminwall)
                    # print('hi',hi,newcoord[hi][0],max(0,currx-500),min(currx+500,tickswall),mindist)
                    for xi in range(max(0, currx-500), min(currx+500, tickswall)):
                        for yi in range(max(0, curry - 500), min(curry + 500, tickswall)):
                            if zwall[xi,yi] == 0:
                                currdist = ((float(xi-curry)*tickx)**2+(float(yi-curry)*ticky)**2)**0.5
                                # print(hi,mindist,currdist,currx,curry,xi,yi,tickx,ticky)
                                if mindist > currdist:
                                    mindist = currdist
                    # print('new mindist=',mindist)
                    distwall[hi] = mindist
                print('calculating isoheight', z[-1], np.where(newcoord[:, 2] == z[-1]))

            print('i', i, value[0], newcoord[i])

        # print('ux=',ux)
        print('newcoord', len(newcoord))
        # print('c[:, 0]',len(c),[float(s) for s in c[i:,0]])

        distwallz = np.zeros(len(distwall))
        for i in range(len(distwallz)):
            distwallz[i] = min(distwall[i], z[i])

        print('test2', value[0])
        print('colorscale', min(z), max(z), z[int(len(z)/2)], ((z[int(len(z)/2)]-min(z))/(max(z)-min(z)))*1.0)

        colorscale = [int(i) for i in ((z-min(z))/(max(z)-min(z)))*5.0]
        colorscale = [int(i) for i in ((np.asarray(z)-min(z))/(max(z)-min(z)))*5.0]
        ssize = 3
        colortitle = 'z'
        plt.figure()
        plt.scatter(ux, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('U_x')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(x, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('x')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(y, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('y')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(z, (model-observation)/observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('z')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(uxuy, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('uxuy')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(uy, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('uy')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(uxuz, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('uxuz')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(uz, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('uz')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(distwall, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('dist wall')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

        plt.figure()
        plt.scatter(distwallz, model-observation, s=ssize, c=colorscale, alpha=0.5)
        plt.ylabel('model-observations')
        plt.xlabel('dist wall+z')
        plt.title('error classification: ' + showfieldtext + ' color: ' + colortitle)
        plt.show()

#    def on_click_choose_file1(self):
#        options = QFileDialog.Options()
#        options |= QFileDialog.DontUseNativeDialog
#        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
#                                                  "All Files (*);;Python Files (*.py)", options=options)
#        if fileName:
#            self.filename1.setText(fileName)


    def on_click_add_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            self.listfiles.addItem(fileName)
        
    def on_click_remove_file(self):

        listItems=self.listfiles.selectedItems()
        if not listItems: return        
        for item in listItems:
           self.listfiles.takeItem(self.listfiles.row(item))
            
    def on_click_para(self):

        textboxtime1value = float(self.itteration1.text())
#        textboxfile1value = self.filename1.text()
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
           for item in listItems:
               files.append(item.text())
        else:
           for item in range(self.listfiles.count()):
              files.append(self.listfiles.item(item).text())
        textboxfile1value = files[0]
        db = get_slice(textboxfile1value, textboxtime1value, int(self.listaxis.currentRow()),
                       float(self.axispos.text()),
                       clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                       clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                       clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()))

        midx = (db['x'].min()+db['x'].max())/2.
        print('midx',midx)        
        midy = (db['y'].min()+db['y'].max())/2.
        print('midy',midy)        
        midz = (db['z'].min()+db['z'].max())/2. # db['z'].min()+(db['z'].max()+db['z'].min())*0.1
        print('midz',midz)        

#        pvsimple.FindSource(self.reader1)
        # pvsimple.LoadState("/ibdata2/nirb/openFOAM/MALA/test1.pvsm")
        # pvsimple.Show()
        # pvsimple.Render()
        # pvsimple.Interact(view=None)
        # pvsimple.WriteImage("testpvstate.png")

        print ('fin test0')

        reader1 = bse.ReadCase('casename', files[0], CaseType='Decomposed Case')  # 'Reconstructed Case')
        slice_1 = pvsimple.Slice(reader1)
        slice_1.SliceOffsetValues = [0.0]
        slice_1.SliceType.Origin = [midx, midy, midz]
        slice_1.SliceType.Normal = [0.0, 0.0, 1.0]

        print ('fin test1')

#        animationScene1 = pvsimple.GetAnimationScene()
#        animationScene1.GoToLast()
        # animationScene1.GoToNext()

        # slicer = servermanager.filters.Slice(Input=self.reader1)
        # slicer.SliceOffsetValues = [0.0]
        # slicer.SliceType.Origin = [155250.0, 563750.0, 180.75]
        # slicer.SliceType.Normal = [0.0, 0.0, 1.0]

        # repClip = servermanager.CreateRepresentation(slicer, view)

        print ('fin test3')
        # get active view
        renderView1 = pvsimple.GetActiveViewOrCreate('RenderView')
        # uncomment following to set a specific view size
        renderView1.ViewSize = [900, 808]
        # show data in view
        slice1Display = pvsimple.Show()
        slice1Display.ColorArrayName = [None, '']

        slice1Display.SetScalarBarVisibility(renderView1, True)
        print ('fin test5')

        # set scalar coloring
        pvsimple.ColorBy(slice1Display, ('POINTS', 'U'))
        # rescale color and/or opacity maps used to include current data range
        slice1Display.RescaleTransferFunctionToDataRange(True)
        # show color bar/color legend
        slice1Display.SetScalarBarVisibility(renderView1, True)
        #
        slice1Display.SetRepresentationType('Wireframe')
        #

        # current camera placement for renderView1
        # renderView1.CameraPosition = [156282.273525518, 563645.495340609, 629.022076989757]
        # renderView1.CameraFocalPoint = [155250.0, 563750.0, 299.749999999999]
        # renderView1.CameraViewUp = [-0.307939724010595, -0.0478789914526784, 0.950200362320365]
        # renderView1.CameraParallelScale = 1069.89254717472

        # pvsimple.Show()
        # pvsimple.Render()
        print ('fin test8')
        pvsimple.WriteImage("testpv.png")
        pvsimple.Interact(view=None)
        print ('fin test')
    def on_click_save_configuaration(self):
        items = []
        for index in range(self.listfields.count()): #xrange
            items.append(self.listfields.item(index).text())
        files = []
        for index in range(self.listfiles.count()):
            files.append(self.listfiles.item(index).text())

        data = {#'file1': self.filename1.text(),
                #'file2': self.filename2.text(),
                'itteration1': self.itteration1.text(),
                'itteration2': self.itteration2.text(),
                'axisValue': self.axispos.text(),
                'axisIndex': self.listaxis.currentRow(),
                'clipxmin': self.xmin.text(),
                'clipxmax': self.xmax.text(),
                'clipymin': self.ymin.text(),
                'clipymax': self.ymax.text(),
                'clipzmin': self.zmin.text(),
                'clipzmax': self.zmax.text(),
                'xcoord': self.xcoord.text(),
                'ycoord': self.ycoord.text(),
                'zcoord': self.zcoord.text(),
                'experiment': self.experiment.currentText(),
                'fieldIndex': self.listfields.currentRow(),
                'heightcontour': self.heightcontour.checkState(),
                'use2': self.use2.checkState(),
                'learnfrom': self.learnfrom.text(),
                'fields': items,
                'files': files
                }

        with open('ofplot.json', 'w') as outfile:
            json.dump(data, outfile, indent=4)

    def changeHoldOn(self):
        print ('clear legend', self.holdon.checkState())
        self.legend = []

    def on_click_density(self):
        textboxtime1value = float(self.itteration1.text())
        textboxfile1value = self.filename1.text()
        if self.heightcontour.checkState():
            db = get_slice_height(textboxfile1value, textboxtime1value, int(self.listaxis.currentRow()),
                                  float(self.axispos.text()),
                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                  reader=self.reader1)
        else:
            db = get_slice(textboxfile1value, textboxtime1value, int(self.listaxis.currentRow()),
                           float(self.axispos.text()),
                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=self.reader1)

        xmin, xmax, ymin, ymax, zi1 = makegrid(db, self.listfields.currentItem().text(), self.listaxis.currentRow())
        total = len(zi1.ravel())
        totalzeros = np.count_nonzero(zi1.ravel())
        QMessageBox.question(self, 'Density', str (totalzeros) + '/' + str(total) + '=' + str(float(totalzeros)/total),
                             QMessageBox.Ok, QMessageBox.Ok)

    def on_click_itter(self):

#        if self.holdon.checkState():
#            fig = plt.figure(110)
#        else:
#            fig = plt.figure()
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
           for item in listItems:
               files.append(item.text())
        else:
           for item in range(self.listfiles.count()):
              files.append(self.listfiles.item(item).text())    
        textboxfile1value = files[0] #self.filename1.text()
#        ritter, rtime = logall(textboxfile1value)
#        # plt.figure()
#        print('itter:',ritter[4],rtime[4],textboxfile1value)
#        plt.plot(ritter, rtime, '.')
#        plt.xlabel('itteration #')
#        plt.ylabel('time')
#        plt.title('itteration time - '+ textboxfile1value)
#        plt.legend()
#        plt.show()
#        print('itter7:',ritter[4],rtime[4],textboxfile1value)
#        
        sel = self.listfields.selectedIndexes()
        corrlast, corrnear, times = corrall(textboxfile1value, sel[0].data(), self.listaxis.currentRow(),
                                            float(self.axispos.text()),
                                            float(self.xmin.text()), float(self.xmax.text()),
                                            float(self.ymin.text()), float(self.ymax.text()),
                                            float(self.zmin.text()), float(self.zmax.text()),
                                            self.heightcontour.checkState())
        if self.holdon.checkState():
            plt.figure(109)
            self.legend.append(textboxfile1value)

        else:
            plt.figure()
        plt.plot(times, corrlast,'-*')
        if self.holdon.checkState():
            plt.legend(self.legend)
            print('self legend',self.legend)
        plt.plot(times, corrnear, 'r')
        plt.legend(['last', 'index'])
        plt.title('Convergence' + textboxfile1value)
        plt.xlabel('itterations')
        plt.ylabel('mean R^2')
        plt.show()
        print ('corrlast', corrlast)
        print ('corrnear', corrnear)
        print ('times', times)
        
###############################3
        mainsimulation = self.experiment.currentText()
        
        if mainsimulation!='None':
            corrobs = []
            observation = chooseobservation(mainsimulation)
    
            showfieldtext = self.listfields.currentItem().text()
    
            for t in times:
                model, coordinates = read_measurement(textboxfile1value, t, mainsimulation,
                                              self.showfieldlist[showfieldtext])
        
                if mainsimulation == 'observationB':
                    observation = observation[:, self.showfieldlist[showfieldtext]]
                elif mainsimulation == 'Michelstadt':
                    observation = observation[:, self.showfieldlist[showfieldtext]]
        
                if float(self.xmin.text()) == float(self.xmax.text()):
                    boolcoordx = (coordinates[:, 0] == coordinates[:, 0])
                else:
                    boolcoordx = (coordinates[:, 0] > float(self.xmin.text())) & (coordinates[:, 0] < float(self.xmax.text()))
                if float(self.ymin.text()) == float(self.ymax.text()):
                    boolcoordy = (coordinates[:, 1] == coordinates[:, 1])
                else:
                    boolcoordy = (coordinates[:, 1] > float(self.ymin.text())) & (coordinates[:, 1] < float(self.ymax.text()))
                if float(self.zmin.text()) == float(self.zmax.text()):
                    boolcoordz = (coordinates[:, 2] == coordinates[:, 2])
                else:
                    boolcoordz = (coordinates[:, 2] > float(self.zmin.text())) & (coordinates[:, 2] < float(self.zmax.text()))
                boolcoord = boolcoordx & boolcoordy & boolcoordz
                newcoord = coordinates[boolcoord]
                observation = observation[boolcoord]
                model = model[boolcoord]
                corrobs.append(rmse(model,observation))
                print ('rmse',t,rmse(model,observation))
            print(times, corrobs)
#            plt.figure()
#            plt.xlabel('itterations')
#            plt.ylabel('mean RMSE m/s')
            plt.plot(times, corrobs,'--*','r')
            plt.title('corr obs')
            plt.show()  

        
        
    def on_click_corr(self):
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
           for item in listItems:
               files.append(item.text())
        else:
           for item in range(self.listfiles.count()):
              files.append(self.listfiles.item(item).text())        
        textboxfile1value = files[0] #self.filename1.text()
        textboxtime1value = float(self.itteration1.text())
        sel = self.listfields.selectedIndexes()
        if self.heightcontour.checkState():
            db = get_slice_height(textboxfile1value, textboxtime1value, int(self.listaxis.currentRow()),
                                  float(self.axispos.text()),
                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                  reader=self.reader1)
        else:
            db = get_slice(textboxfile1value, textboxtime1value, int(self.listaxis.currentRow()),
                           float(self.axispos.text()),
                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()))
        plt.figure()
        ax1 = plt.subplot(311)
        if self.listaxis.currentRow() == 0:
            plt.tricontourf(db['y'], db['z'], db[sel[0].data()], 100)
        if self.listaxis.currentRow() == 1:
            plt.tricontourf(db['x'], db['z'], db[sel[0].data()], 100)
        if self.listaxis.currentRow() == 2:
            plt.tricontourf(db['x'], db['y'], db[sel[0].data()], 100)

        plt.colorbar()
        plt.title(textboxfile1value+'-'+sel[0].data())
        plt.show()

        # Plot density map.
        plt.subplot(312, sharex=ax1, sharey=ax1)

        if self.listaxis.currentRow() == 0:
            plt.tricontourf(db['y'], db['z'], db[sel[1].data()], 100)  # db[sel[1].data()], 100
        if self.listaxis.currentRow() == 1:
            plt.tricontourf(db['x'], db['z'], db[sel[1].data()], 100)
        if self.listaxis.currentRow() == 2:
            plt.tricontourf(db['x'], db['y'],db[sel[1].data()], 100)
        print('end BBB2')

        plt.colorbar()
        plt.title(textboxfile1value+'-'+sel[1].data())
        plt.show()
        print('end BBB3')

        plt.subplot(313)
        plt.scatter(db[sel[0].data()], db[sel[1].data()], s=1)
        plt.xlabel(sel[0].data())
        plt.ylabel(sel[1].data())
        plt.title('corr'+str(np.corrcoef(db[sel[0].data()], db[sel[1].data()])[1, 0]))
        plt.show()
       
        
    def on_click_profile(self):
        ##############################
        # calculating vertical profile
        ##############################
        textboxfile1value = self.filename1.text()
        textboxtime1value = float(self.itteration1.text())

        sel = self.listfields.selectedIndexes()
        if self.heightcontour.checkState():
            db = get_slice_height(textboxfile1value, textboxtime1value, 1, float(self.axispos.text()),
                                  offset=range(0, 200, 1),
                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                  reader=self.reader1)
        else:
            db = get_clip(textboxfile1value, [textboxtime1value],
                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=self.reader1)

        mainsimulation = self.experiment.currentText()
        showfieldtext = self.listfields.currentItem().text()

        if mainsimulation!='None':
            model, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation,
                                                  self.showfieldlist[showfieldtext])

        ticks = 150j  # 600
        z = np.zeros(int(ticks.imag))
        unz = np.zeros_like(z)
        roof = np.zeros_like(z)

        ap3 = db[['x', 'y', 'z', 'U_x']]  # at the ground (1 meter)
        points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
        values = ap3['U_x']
        xi, yi, zi = np.mgrid[ap3['x'].min():ap3['x'].max():ticks, ap3['y'].min():ap3['y'].max():ticks,
                     ap3['z'].min():ap3['z'].max():ticks]
        print('start griddata cubic', datetime.datetime.now())
        # grid_z1 = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(values), (xi, yi, zi))  # Nearest for keeping zeros that are buildings
        print('mid griddata cubic', datetime.datetime.now())
        grid_z1n = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(values), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        print('end griddata cubic', datetime.datetime.now())
        grid_z1 = grid_z1n

        grid_z1[grid_z1n==0]=0
        z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticks.imag))
        for i in range(len(z)):
            if z[i]<2.5:  # 1.5 meter
                ground = i
        ground += 1
        unz[0] = 0
        for i in range(1,int(ticks.imag)):
            grid_slice = grid_z1[:, :, i]
            # unz[i] = grid_slice[grid_slice != 0].mean()
            # print('unzlen i:', i)
            unz1 = db[db['z'] > z[i - 1]]
            unz2 = unz1[unz1 < z[i]]
            unz3 = unz2.groupby('z')[sel[0].data()].mean()
            if len(unz3)>0:
                unz[i] = np.average(unz3)
                # unz[i] = np.average(db[db[db['z']>z[i-1]]<z[i]].groupby('z')[sel[0].data()].mean())
                # print('unzlen0', unz[i], grid_slice[grid_slice != 0].mean())
            else:
                unz[i] = unz[i-1]
            roof[i] = (grid_z1[:, :, i] == 0).sum()
            # print ('3d:', i, z[i], roof[i], unz[i], (grid_z1[:, :, i] == 0).sum(), (grid_z1[:, :, i] != 0).mean())

        h = 0
        for i in range(len(z)-2,ground,-1):
            h += z[i] * (roof[i]-roof[i+1])
            # print('III',i, z[i], h, roof[i])

        h /= roof[ground]

        for i in range(len(z)-1,1,-1):
            if z[i] > h:
                hindex = i

        hindex += 1
        lambdap = float((grid_z1[:, :, ground] == 0).sum()) / (grid_z1[:, :, ground] == grid_z1[:, :, ground]).sum()


        # lamdaF start

        floor = grid_z1[:, :, ground].copy()
        floor = flood_main(floor)
        maxfloor = floor.max()
        lambdaf = 0
        for i in range(1, maxfloor):
            # index for height
            c = np.where(floor == i)
            height = np.nonzero(grid_z1[c[0][0], c[1][0], :])[0][0]
            position = np.where(floor == i)
            posymin = min(position[:][1])  # maybe 0 ???
            posymax = max(position[:][1])  # maybe 0 ???
            lambdaf += (posymax-posymin)*height
        lambdaf /= float((grid_z1[:, :, ground] == grid_z1[:, :, ground]).sum())
        print ('LambdaF =>', lambdaf)

        # plt.figure()
        # plt.imshow(floor)
        # plt.colorbar()
        # plt.show()


        for i in range(int(ticks.imag)):
            if i == 130:
                grid_slice = grid_z1[i, :, :]
                down = grid_slice[:, ground]
                downud = np.flipud(down)
                wall = []
                wallud = []
                for j in range(int(ticks.imag)-1):
                    if down[j] == 0 and down[j+1] !=0:
                        wall.append(j)
                    if downud[j] == 0 and downud[j + 1] != 0:
                        wallud.append(len(down)-j-1)
                # print('down:')
                # for j in range(len(down)):
                    # print(j, len(down)-j-1, down[j])
                # print('wall1:',wall)
                # print('wall2:',wallud)
                # for j in range(len(wall)):
                    # print(j,wallud[len(wall)-j-1],wall[j])
                    # print('upthere',grid_slice[int((wallud[len(wall)-j-1]+wall[j])/2),:])
                    # print('height',np.nonzero(grid_slice[int((wallud[len(wall)-j-1]+wall[j])/2),:])[0][0])

        # lambdaF end

        print ('hindex', hindex, h, lambdap, ground)

        if mainsimulation!='None':
            zunique = np.unique(coordinates[:, 2])
            zvalue = np.zeros_like(zunique)
            for i in range(len(zunique)):
                zvalue[i] = np.average(model[coordinates[:, 2] == zunique[i]])

        # lambdap = 0.35  # to calculate
        # lambdaf = 0.25  # to calculate
        lc = (h*axisscale)*(1-lambdap)/lambdaf
        uh = unz[hindex]   # 4.6

        # ustar = beta * uh  ==> beta = ustar /uh
        tmp = hindex - 1

        mu = 1.82E-5
        dudy = (unz[1]-unz[0])/((z[1]-z[0])*axisscale)
        # for j in range(10):
        #     print('jjjj',j,z[j],unz[j])
        tau = math.fabs(mu * dudy)
        # eq 2.28 at https://www.springer.com/cda/content/document/cda_downloaddocument/9783540882534-c1.pdf?SGWID=0-0-45-710813-p173848437
        tau = 1.225*0.4*0.4*z[tmp]**2*((unz[tmp+1]-unz[tmp-1])/(z[tmp+1]-z[tmp-1]))**2
        ro = 1.225
        ustar = math.sqrt(tau/ro)

        beta = 0.13
        beta = ustar / uh

        print ('beta',beta, ustar, uh, z[ground], z[hindex])
#        print ('tau',1.225*0.4*0.4*z[tmp]**2*((unz[tmp+1]-unz[tmp-1])/(z[tmp+1]-z[tmp-1]))**2)
#
#        lvar = 2.0*beta**3.0*lc * 0.3
#        k = 0.41
#        z0 = h*0.1
#        d = h*0.3
#        d = lvar/k
#        z0 = lvar/k*math.exp(-k/(beta))  # 1.0, 0.8, 0.2
##        level1 = len(z)*0.8
##        level2 = len(z)*0.6
#        # ustar1 = unz[level1] * k / math.log(z[level1/z0)
#        # ustar2 = unz[level2] * k / math.log(z[level2/z0)
#        # ustar = 0.5 * (ustar1 + ustar2)
#        # beta = ustar / uh
#        # ustar = beta * uh
#        # print('ustar',ustar1, ustar2, 'beta', beta)
#
#        # d = 14.52
#        # d = 2.52
#        # z0= 0.85
#        ## d = 1.E-8
#        ## z0 = 320
#        uprofile = np.zeros(len(z))
#        print('h',h,'uh',uh,'lvar',lvar,'k',k,'beta',beta,'d',d,'zo',z0)
#        for i in range(len(z)):
#            if z[i] <= h:
#                uprofile[i] = uh * math.exp(beta*((z[i]-h)*axisscale)/lvar)
#                # print(i, 'exp')
#            else:
#                uprofile[i] = uh*beta/k*math.log(((z[i]-h)*axisscale+d)/z0)
#            # print(i,z[i],z[i]-h,uprofile[i])
#            # print(i, z[i], h, u[i])
#        # params, params_covariance = curve_fit(test_profile, z[int(hindex):]-h, unz[int(hindex*1):],
#        #                                                p0=[uh*beta/k, 1/z0, d/z0],bounds=(0, [9993., 99991., 99990.5]))
#        # params, params_covariance = curve_fit(log_profile, z[int(hindex*1):]-h, unz[int(hindex*1):],
#        #                                                p0=[uh*beta/k, 1/z0, 0, d/z0])
#        # params1, params_covariance1 = curve_fit(exp_profile, z[ground:hindex]-h, unz[ground:hindex],
#        #                                                p0=[uh, -beta/lvar, 0])
#        # params2, params_covariance2 = curve_fit(logexp_profile, z[:]-h, unz[:],
#        #                                                p0=[uh*beta/k, 1/z0, d/z0, uh, -beta/lvar, 1])
#        # print('exp profile',params1[0], params1[1], params1[2])
#        # print('log profile',params[0], params[1], params[2], params[3])
#        # print('logexp profile',params2[0], params2[1], params2[2], params2[3],params2[4], params2[5])
#        # print(len(z),'guess',[uh*beta/k, 1/z0, d/z0, uh, -beta/lvar, 1])

        if self.holdon.checkState():
            fig = plt.figure(100)
        else:
            fig = plt.figure()
        if self.heightcontour.checkState():
            plt.plot(db.groupby('HeightFromTopo')[sel[0].data()].mean(),
                     db.groupby('HeightFromTopo')['HeightFromTopo'].mean(), label='db height')
        else:
            plt.plot(db.groupby('z')[sel[0].data()].mean(), db.groupby('z')['z'].mean(), 'b.', markersize=2, label='db')
#        print ('assuming building height is ',h,'->', hindex, 'lambdaP:', lambdap, 'lambdaF:', lambdaf, 'ground:', ground)

        # plt.plot(test_profile(z[hindex:]-h, params[0], params[1], params[2]),z[hindex:],  'c.',
        #          label='Fitted function')
        # plt.plot(exp_profile(z[:hindex]-h, params1[0], params1[1], params1[2]), z[:hindex], 'c--',
        #          label='Fitted function')

        # plt.plot(uh * beta / k * np.log(((z[hindex:] - h) + d) / z0), z[hindex:], 'g.', label='Fitted function - Eyal')
        # plt.plot(uh * np.exp(beta * (z[:hindex] - h) / lvar), z[:hindex], 'g--', label='Fitted function - Eyal')
        # plt.plot(uprofile[hindex:], z[hindex:], 'g.', label='Fitted function - Eyal')
        # plt.plot(uprofile[:hindex], z[:hindex], 'g--', label='Fitted function - Eyal')

        plt.plot(unz, z, 'r.', label='gridded')
        # plt.plot(unz, z, 'r*')
        if mainsimulation!='None':
            plt.plot(zvalue, zunique, 'k.', label='observations ?')
        title = shrinktitle(textboxfile1value) + '(' + sel[0].data() + '):' + self.xmin.text() + ',' + self.xmax.text() + ',' + self.ymin.text() + ',' + self.ymax.text()
        plt.title(title)
        if self.holdon.checkState():
            self.legend.append(title)
            plt.legend(self.legend)
        plt.legend()
        fig.canvas.draw()
        plt.figure()
        # plt.plot(z[hindex:], uprofile[hindex:], 'g.', label='Fitted function - Eyal')
        # plt.plot(z[:hindex], uprofile[:hindex], 'g--', label='Fitted function - Eyal')
        plt.plot(z,unz, 'r.', label='gridded')
        # plt.plot(z[hindex:], test_profile(z[hindex:]-h, params[0], params[1], params[2]), 'c.',
        #          label='Fitted function')
#        print(uh * beta / k * math.log(((z[i] - h) * axisscale + d) / z0 ) )
        # print('uh * beta / k', params[0], uh * beta / k)
        # print('1/z0', params[1], 1./z0, z0)
        # print('d/z0', params[2], d/z0, d)
        # print('c', params[3])

        # return a * np.log(b * x + c) [+ d]

        plt.show()

        
        if self.holdon.checkState():
            fig = plt.figure(160)
            plt.plot(unz,z, '.', label='gridded')
            plt.legend(self.legend, loc = 'upper left')
#            fig.canvas.draw()       
            plt.show()
        
        
    def on_click_profile2(self):
        ##############################
        # calculating vertical profile
        ##############################
        textboxfile1value = self.filename1.text()
        textboxtime1value = float(self.itteration1.text())
        textboxfile2value = self.filename2.text()

        sel = self.listfields.selectedIndexes()
        reader2 = bse.ReadCase(textboxfile2value, textboxfile2value, CaseType='Decomposed Case')  # ' Reconstructed Decomposed Case')

        if self.heightcontour.checkState():
            db = get_slice_height(textboxfile1value, textboxtime1value, 1, float(self.axispos.text()),
                                  offset=range(0, 200, 1),
                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                  reader=self.reader1)
            db2 = get_slice_height(textboxfile2value, textboxtime1value, 1, float(self.axispos.text()),
                                  offset=range(0, 200, 1),
                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                  reader=reader2)
        else:
            db = get_clip(textboxfile1value, [textboxtime1value],
                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=self.reader1)
            db2 = get_clip(textboxfile2value, [textboxtime1value],
                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()), reader=reader2)

        mainsimulation = self.experiment.currentText()
        showfieldtext = self.listfields.currentItem().text()

        if mainsimulation!='None':
            model, coordinates = read_measurement(textboxfile1value, textboxtime1value, mainsimulation,
                                                  self.showfieldlist[showfieldtext])

        ticks = 150j  # 600
        z = np.zeros(int(ticks.imag))


        ap3 = db[['x', 'y', 'z', 'U_x']]  # at the ground (1 meter)
        ap32 = db2[['x', 'y', 'z', 'U_x']]  # at the ground (1 meter)
        points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
        values = ap3['U_x']
        points32 = ap32[['x', 'y', 'z']].values  # numpy[cells,3]
        values2 = ap32['U_x']
        xi, yi, zi = np.mgrid[ap3['x'].min():ap3['x'].max():ticks, ap3['y'].min():ap3['y'].max():ticks,
                     ap3['z'].min():ap3['z'].max():ticks]
        print('start griddata cubic', datetime.datetime.now())
        # grid_z1 = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(values), (xi, yi, zi))  # Nearest for keeping zeros that are buildings
        print('mid griddata cubic', datetime.datetime.now())
        grid_z1n = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(values), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        grid_z1n2 = griddata((points32[:, 0], points32[:, 1], points32[:, 2]), np.asarray(values2), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        print('end griddata cubic', datetime.datetime.now())
        cells = 40
        print('cells size',ap3['x'].min(), ap3['x'].max(),ap3['y'].min(), ap3['y'].max())

        z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticks.imag))
        x = np.linspace(ap3['x'].min(), ap3['x'].max(), int(ticks.imag))
        y = np.linspace(ap3['y'].min(), ap3['y'].max(), int(ticks.imag))

        if mainsimulation!='None':
            zunique = np.unique(coordinates[:, 2])
            zvalue = np.zeros_like(zunique)
            for i in range(len(zunique)):
                zvalue[i] = np.average(model[coordinates[:, 2] == zunique[i]])
                
        title1 = shrinktitle(textboxfile1value)       
        title2 = shrinktitle(textboxfile2value)       
                
        map1 = np.zeros([cells, cells])
        map2 = np.zeros([cells, cells])
        for k in range(cells):
            for j in range(cells):
                print('A0',  datetime.datetime.now())
                ci1 = k*int((grid_z1n.shape[0]/float(cells)))
                ci2 = (k+1)*int((grid_z1n.shape[0]/float(cells)))-1
                cj1 = j*int((grid_z1n.shape[1]/float(cells)))
                cj2 = (j+1)*int((grid_z1n.shape[1]/float(cells)))-1
#                print('A:',grid_z1n.shape[0],k,ci1,ci2,grid_z1n.shape[1],j,cj1,cj2)
                grid_z1 = grid_z1n[ci1:ci2,cj1:cj2,:]
                ci1 = k*int((grid_z1n2.shape[0]/float(cells)))
                ci2 = (k+1)*int((grid_z1n2.shape[0]/float(cells)))-1
                cj1 = j*int((grid_z1n2.shape[1]/float(cells)))
                cj2 = (j+1)*int((grid_z1n2.shape[1]/float(cells)))-1
#                print('B:',grid_z1n2.shape[0],k,ci1,ci2,grid_z1n2.shape[1],j,cj1,cj2)
                grid_z12 = grid_z1n2[ci1:ci2,cj1:cj2,:]
        
#                grid_z1[grid_z1==0]=0
#                grid_z12[grid_z12==0]=0
#                grid_z1[grid_z1n==0]=0
#                grid_z12[grid_z1n2==0]=0

                x1 = x[(k)*int((int(ticks.imag)/float(cells)))]
                x2 = x[(k+1)*int((int(ticks.imag)/float(cells)))-1]
                y1 = y[(j)*int((int(ticks.imag)/float(cells)))]
                y2 = y[(j+1)*int((int(ticks.imag)/float(cells)))-1]
                print('A:', k*(cells)+j, datetime.datetime.now(), x1,x2,y1,y2)
                unz = np.zeros_like(z)
                un2z = np.zeros_like(z)
                for i in range(1,int(ticks.imag)):
#                    unz1 = db[(db['z'] > z[i - 1]) & (db['z'] < z[i])]
##                    unz2 = unz1[unz1 < z[i]]
#                    unz1 = unz1[(unz1['x']>=x1) & (unz1['x']<x2) & (unz1['y']>=y1) & (unz1['y']<y2)]
##                    unz1 = unz1[unz1['x']<x2]
##                    unz1 = unz1[(unz1['y']>=y1) & (unz1['y']<y2)]
##                    unz1 = unz1[unz1['y']<y2]
#                    unz3 = unz1.groupby('z')[sel[0].data()].mean()
#                    if len(unz3)>0:
#                        unz[i] = np.average(unz3)
#                    else:
#                        unz[i] = unz[i-1]
#
##                    un2z1 = db2[db2['z'] > z[i - 1]]
#                    un2z1 = db2[(db2['z'] > z[i - 1]) & (db2['z'] < z[i])]
##                    un2z2 = un2z1[un2z1 < z[i]]                   
##                    un2z1 = un2z1[un2z1['x']>=x1]
#                    un2z1 = un2z1[(un2z1['x']>=x1) & (un2z1['x']<x2) & (un2z1['y']>=y1) & (un2z1['y']<y2)]
##                    un2z1 = un2z1[(un2z1['y']>=y1) & (un2z1['y']<y2)]
##                    un2z1 = un2z1[un2z1['y']<y2]
#                    un2z3 = un2z1.groupby('z')[sel[0].data()].mean()
#                    if len(un2z3)>0:
#                        un2z[i] = np.average(un2z3)
#                    else:
#                        un2z[i] = un2z[i-1]
                    g1 = grid_z1[:,:,i]
                    g2 = grid_z12[:,:,i]
                    unz[i] = g1[g1!=0.].mean()
                    un2z[i] = g2[g2!=0.].mean()
        
                print('A2',  datetime.datetime.now())

                if cells > 3:
                    if self.holdon.checkState():
                        fig = plt.figure(100)
                        plt.clf()
                        self.legend=[]
                    else:
                        fig = plt.figure()
                    plt.plot(unz, z, 'r.', label=title1)
                    plt.plot(un2z, z, 'b.', label=title2)
                    # plt.plot(unz, z, 'r*')
                    if mainsimulation!='None':
                        plt.plot(zvalue, zunique, 'k.', label='observations ?')
                    title = title1 + '(' + sel[0].data() + '):' + str(x1) + ',' + str(x2) + ',' + str(y1) + ',' + str(y2)
                    plt.title(title)
                    if self.holdon.checkState():
                        self.legend.append(title)
                        plt.legend(self.legend)
                    plt.legend()
                    fig.canvas.draw()
                    fig.savefig('profile'+str(cells)+'-'+str(k)+'-'+str(j)+'.png')
                    plt.close()
                count1 = 0
                sumd = 0.
                for m in range(1,len(unz)):
                    sumd += math.fabs(un2z[m]-unz[m])
#                    print (sumd)
                    if (math.fabs(unz[m]-un2z[m])>=1.):
                        count1 += 1
#                    print (count1)
                map1[k,j] = count1
                map2[k,j] = sumd
#                print('A5',  datetime.datetime.now())

        print('test')
        plt.figure()
        plt.imshow(map1, origin='lower', interpolation = 'nearest')
        plt.colorbar()
        plt.title('count diff > 1')
        plt.show()

        plt.figure()
        plt.imshow(map2, origin='lower', interpolation = 'nearest')
        plt.colorbar()
        plt.title('sum (diff)')
        plt.show()


        
    def on_click_profile3(self):
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
           for item in listItems:
               files.append(item.text())
        else:
           for item in range(self.listfiles.count()):
              files.append(self.listfiles.item(item).text())

        axisindex = int(self.listaxis.currentRow())
        axispos = float(self.axispos.text())
        learnfrom=self.learnfrom.text().strip()
        textboxtime1value = float(self.itteration1.text())
        startarea = int(self.xcoord.text())
        endarea = int(self.ycoord.text())
        currentitem = self.listfields.currentItem().text()
        xmingui = float(self.xmin.text())
        ymingui = float(self.ymin.text())
        zmingui = float(self.zmin.text())
        xmaxgui = float(self.xmax.text())
        ymaxgui = float(self.ymax.text())
        zmaxgui = float(self.zmax.text())
        
        areaname='yehuda200'
        
        profile3(files,textboxtime1value,axisindex,axispos,learnfrom,areaname,startarea,endarea,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui)

            
    def on_click_update_fields(self):
        textboxfile1value = self.listfiles.item(0).text()
        reader = bse.ReadCase('slice_1', textboxfile1value, CaseType='Decomposed Case')
        print('update:', reader.CellArrays)
        # self.listfields.clear()
        # for field in range(len(reader.CellArrays)):
        #     self.listfields.addItem(reader.CellArrays[field])
        slice_1 = pvsimple.Slice(reader)
        slice_1.SliceType = 'Plane'
        slice_1.SliceOffsetValues = [0.0]
        slice_1.SliceType.Origin = [0.0, 0.0, 2.0]
        slice_1.SliceType.Normal = [0.0, 0.0, 1.0]
        slice_1.UpdatePipeline()
        textboxtime1value = float(self.itteration1.text())

        itr = bse.to_pandas('reader', slice_1, timelist=textboxtime1value)
        dbb = itr.next()
        db = dbb[dbb.keys()[0]]
        print ('AAAA', db.keys())
        print ('AAABB', db)
        print ('AAACC', itr.next)
        self.listfields.clear()
        for field in range(len(db.keys())):
            self.listfields.addItem(db.keys()[field])

    def on_click_diff_files(self):
        if self.use2.checkState():
            itter = [float(self.itteration1.text())]
        else:
            itter = [float(self.itteration1.text()),float(self.itteration2.text())]
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
            for item in listItems:
               files.append(item.text())
        else:
            for item in range(self.listfiles.count()):
                files.append(self.listfiles.item(item).text())
                        
        field = self.listfields.currentItem().text()
        xminc = float(self.xmin.text())
        xmaxc = float(self.xmax.text())
        yminc = float(self.ymin.text())
        ymaxc = float(self.ymax.text())
        zminc = float(self.zmin.text())
        zmaxc = float(self.zmax.text())
        axisindex = int(self.listaxis.currentRow())
        axispos = float(self.axispos.text())
        currentRow = self.listaxis.currentRow()
        heightcontour = self.heightcontour.checkState()
        print('diff_files', files, [itter])
        
        diff_files(files, itter, field, xminc, xmaxc, yminc, ymaxc, zminc, zmaxc, axisindex, axispos, currentRow, parallel=self.parallel.checkState())
        
    def on_click_plot(self):
        items = self.listfields.selectedItems()
        fields = []
        for i in range(len(items)):
            fields.append(str(self.listfields.selectedItems()[i].text()))
            
        itter = float(self.itteration1.text())
        files=[]        
        listItems=self.listfiles.selectedItems()
       
        if listItems:           
            for item in listItems:
               files.append(item.text())
        else:
            for item in range(self.listfiles.count()):
                files.append(self.listfiles.item(item).text())
        axisindex = int(self.listaxis.currentRow())
        axispos = float(self.axispos.text())
        addtotitle=self.addtotitle.text()
        xminc = float(self.xmin.text())
        xmaxc = float(self.xmax.text())
        yminc = float(self.ymin.text())
        ymaxc = float(self.ymax.text())
        zminc = float(self.zmin.text())
        zmaxc = float(self.zmax.text())
        
        if self.parallel.checkState():
            parallel=True
        else:
            parallel=False
        print('parallel plot', parallel)
        
        print(files, fields, itter, xminc, xmaxc, yminc, ymaxc, zminc, zmaxc, axisindex, axispos, addtotitle)
        plot_file(files, fields, itter, xminc, xmaxc, yminc, ymaxc, zminc, zmaxc, axisindex, axispos, addtotitle, parallel=parallel)


    def on_click_diff_times(self):
        files=[]        
        listItems=self.listfiles.selectedItems()
        if listItems:
            for item in listItems:
               files.append(item.text())
        else:
            for item in range(self.listfiles.count()):
                files.append(self.listfiles.item(item).text())
        textboxfile1value = files[0]
        textboxtime1value = float(self.itteration1.text())
        textboxtime2value = float(self.itteration2.text())
        
        print(textboxfile1value, [textboxtime1value, textboxtime2value], int(self.listaxis.currentRow()),
                               float(self.axispos.text()),
                               float(self.xmin.text()), float(self.xmax.text()),
                               float(self.ymin.text()), float(self.ymax.text()),
                               float(self.zmin.text()), float(self.zmax.text()))
    
        
        if self.heightcontour.checkState():
            db = get_slice_height(textboxfile1value, textboxtime1value,
                                  int(self.listaxis.currentRow()), float(self.axispos.text()),
                                  clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                                  clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                                  clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()),
                                  reader=self.reader1)
        else:
            
            db = get_slice(textboxfile1value, [textboxtime1value, textboxtime2value], int(self.listaxis.currentRow()),
                           float(self.axispos.text()),
                           clipxmin=float(self.xmin.text()), clipxmax=float(self.xmax.text()),
                           clipymin=float(self.ymin.text()), clipymax=float(self.ymax.text()),
                           clipzmin=float(self.zmin.text()), clipzmax=float(self.zmax.text()))
        xmin, xmax, ymin, ymax, zi1 = makegrid(db[0], self.listfields.currentItem().text(), self.listaxis.currentRow())
        xmin, xmax, ymin, ymax, zi2 = makegrid(db[1], self.listfields.currentItem().text(), self.listaxis.currentRow())
        zmin = min(zi1.min(), zi2.min())
        zmax = max(zi1.max(), zi2.max())
        # Plot density map.
        plt.figure()
        ax1 = plt.subplot(221)
        plt.imshow(
            zi1, extent=(xmin, xmax, ymin, ymax), origin='lower',
            cmap=plt.get_cmap('GnBu_r'))
        if ~np.isnan(zmin):
            plt.clim(zmin, zmax)
        plt.title('(' + str(textboxtime1value)+')')
        plt.colorbar()

        # Plot density map.
        plt.subplot(222, sharex=ax1, sharey=ax1)
        plt.imshow(
            zi2, extent=(xmin, xmax, ymin, ymax), origin='lower',
            cmap=plt.get_cmap('GnBu_r'))
        if ~np.isnan(zmin):
            plt.clim(zmin, zmax)
        plt.title('(' + str(textboxtime2value) + ')')
        plt.colorbar()

        plt.subplot(223, sharex=ax1, sharey=ax1)
        plt.imshow(
            zi2-zi1, extent=(xmin, xmax, ymin, ymax), origin='lower',
            cmap=plt.get_cmap('GnBu_r'))
        plt.title('dbb-db')
        plt.colorbar()
        plt.subplot(224)
        plt.scatter(zi2.ravel(), zi1.ravel(), s=1)
        line11 = np.linspace(zmin, zmax)
        plt.scatter(line11, line11, s=1, c='r')
        plt.xlabel(textboxtime2value)
        plt.ylabel(textboxtime1value)
        zi1nan = np.isnan(zi1.ravel())
        zi2nan = np.isnan(zi2.ravel())
        corr = np.corrcoef(zi1.ravel()[~(zi2nan | zi1nan)], zi2.ravel()[~(zi2nan | zi1nan)])[1, 0]        
        plt.title('corr'+str(corr))
        plt.show()


def buildtopo():
        
    yehuda200 = choosearea('tlvbig250')
    topostart = [
    "FoamFile\n"
    "{\n"
    "version 2.0;\n"
    "format ascii;\n"
    "class dictionary;\n"
    "object topoSetDict;\n"
    "}\n"
    "\n"
    "\n"
    "actions\n"
    "(\n"
    ]
    topomid=[]
    topoend=[");\n"]
    
    for i in range(len(yehuda200)):
        topotmp=[
      "{\n"
        "name    grid",str(i),";\n"
        "type    cellSet;\n"
        "action  new;\n"
        "source  boxToCell;\n"
        "sourceInfo\n"
        "{\n"
        "name    lower;\n"
          "box (",str(yehuda200[i][0])," ",str(yehuda200[i][2])," 0) (",str(yehuda200[i][1])," ",str(yehuda200[i][3])," ",str(yehuda200[i][6]),");\n"
        "}\n"
      "}\n"
                ]
        topomid.append(''.join(topotmp))
         
    file1 = open('topoSetDict', 'w')
    file1.writelines(topostart)
    file1.writelines(topomid)
    file1.writelines(topoend)
    file1.close()
    
    fvoptionstart= [
    "/*--------------------------------*- C++ -*----------------------------------*\\n"
    "| =========                 |                                                 |\n"
    "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
    "|  \\    /   O peration     | Version:  v1912                                 |\n"
    "|   \\  /    A nd           | Website:  www.openfoam.com                      |\n"
    "|    \\/     M anipulation  |                                                 |\n"
    "\*---------------------------------------------------------------------------*/\n"
    "FoamFile\n"
    "{\n"
        "version     2.0;\n"
        "format      ascii;\n"
    "    class       dictionary;\n"
    "    location    \"constant\";\n"
    "    object      fvOptions;\n"
    "}\n"
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n  "
            ]
    
    fvoptionmid=[]
    
    for i in range(len(yehuda200)):
        if yehuda200[i][6]==0.0:
            lc1=0.0
        else:
            lc1=yehuda200[i][4]/(1.-yehuda200[i][5])/yehuda200[i][6] # 1/lc=lf/(1-lp)/he
#        lc1=yehuda200[i][4]/(1.-yehuda200[i][5])/yehuda200[i][6]/1. # 1/lc=lf/(1-lp)/he
    #    lc1=yehuda200[72][4]/(1.-yehuda200[72][5])/yehuda200[72][6] # 1/lc=lf/(1-lp)/he
        fvoptiontmp=[
    "grid"+str(i)+"\n"
    "{\n"
        "type vectorPorous_fvOptions;\n"
        "active on;\n"
            "selectionMode   cellSet;\n"
            "cellSet         grid"+str(i)+";\n"
        "volumeMode specific;\n"
        "injectionRateSuSp\n"
            "{\n"
                  "U (( 0.00000 0.0000 0 ) "+ str(lc1)+");\n"
            "}\n"
    "}\n"
            ]
    
        fvoptionmid.append(''.join(fvoptiontmp))
         
    file1 = open('fvOptions', 'w')
    file1.writelines(fvoptionstart)
    file1.writelines(fvoptionmid)
    file1.close()
    
    
#     end def buildtopo
    
    
def remapsimulation(fromsim, tosim, itter = 14000):
    fromsim= u'/data4bk/nirb/Simulations/Dans/tlvbig2/'
    tosim = u'/data4bk/nirb/Simulations/Dans/tlvbigmap2/'
    itter = 14000
    
    print ('start saving mapped results')
    print ('we will start with taking the origin file and interpolate it to a regular fine grid with nearest method so it will catch all the zero values inside the buildings')
    print ('!!! we need to decide if we take the cell center or the borders')
    print ('than we will interpolate the fine grid to the coarse grid using mean method')
              
    full = False
    print('start remap',fromsim, tosim, datetime.datetime.now())

    getdatafromtext=True
    
    if getdatafromtext is True:
        db = get_clip(fromsim, [itter],
                  clipxmin=float(0), clipxmax=float(0),
                  clipymin=float(0), clipymax=float(0),
                  clipzmin=float(0), clipzmax=float(0))
        
        turbulence = 'epsilon' #'omega'
    
        if full:
            ap3 = db[0][['x', 'y', 'z', 'U_x', 'U_y', 'U_z','p','k','nut',turbulence]]  # at the ground (1 meter)
        else:
            ap3 = db[0][['x', 'y', 'z', 'U_x', 'U_y', 'U_z']]  # at the ground (1 meter)
        points3 = ap3[['x', 'y', 'z']].values  # numpy[cells,3]
        # values = ap3[sel[0].data()]
        valuesux = ap3['U_x']
        valuesuy = ap3['U_y']
        valuesuz = ap3['U_z']
        if full:
            valuesp = ap3['p']
            valuesk = ap3['k']
            valueso = ap3[turbulence]
            valuesn = ap3['nut']
        xmin=ap3['x'].min()
        xmax=ap3['x'].max()
        ymin=ap3['y'].min()
        ymax=ap3['y'].max()
        zmin=ap3['z'].min()
        zmax=ap3['z'].max()
    
        # ticks = 600j
        grid_points = 27000000
        meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)
    
        resolution = (meters / grid_points)**(1./3)
        if resolution==0:
            print('no clip boundaries in one of the axis')
        # resolution = 0.1
        ticksx = int((xmax-xmin) / resolution) * 1j
        ticksy = int((ymax-ymin) / resolution) * 1j
        ticksz = int((zmax-zmin) / resolution) * 1j
        print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)
    
        xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]
        # grid = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(values), (xi, yi, zi),
        #                    method='nearest')  # Nearest for keeping zeros that are buildings
        gridux = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesux), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        griduy = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuy), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        griduz = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesuz), (xi, yi, zi),
                           method='nearest')  # Nearest for keeping zeros that are buildings
        if full:
            gridrp = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesp), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            gridrk = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesk), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            gridro = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valueso), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
            gridrn = griddata((points3[:, 0], points3[:, 1], points3[:, 2]), np.asarray(valuesn), (xi, yi, zi),
                               method='nearest')  # Nearest for keeping zeros that are buildings
    
        x = np.linspace(ap3['x'].min(), ap3['x'].max(), int(ticksx.imag))
        y = np.linspace(ap3['y'].min(), ap3['y'].max(), int(ticksy.imag))
        z = np.linspace(ap3['z'].min(), ap3['z'].max(), int(ticksz.imag))
    else: # getdatafromtext is False
        print('start to read U file in the time directory')
        slash = textboxfile1value.rfind("/") + 1
        directory = fromsim +str(itter)+"/"
        filex = directory + "Cx"
        filey = directory + "Cy"
        filez = directory + "Cz"
        fileu = directory + "U"
        if full:
            filep = directory + "p"
            filek = directory + "k"
            fileo = directory + "omega"
            filen = directory + "nut"
        fhu = open(fileu, 'r')
        fhx = open(filex, 'r')
        fhy = open(filey, 'r')
        fhz = open(filez, 'r')
        if full:
            fhp = open(filep, 'r')
            fhk = open(filek, 'r')
            fho = open(fileo, 'r')
            fhn = open(filen, 'r')
        txtu = fhu.read()
        txtx = fhx.read()
        txty = fhy.read()
        txtz = fhz.read()
        if full:
            txtp = fhp.read()
            txtk = fhk.read()
            txto = fho.read()
            txtn = fhn.read()
        txtusplit = txtu.split()
        txtxsplit = txtx.split()
        txtysplit = txty.split()
        txtzsplit = txtz.split()
        if full:
            txtpsplit = txtp.split()
            txtksplit = txtk.split()
            txtosplit = txto.split()
            txtnsplit = txtn.split()
        posu = txtusplit.index("internalField")
        posx = txtxsplit.index("internalField")
        posy = txtysplit.index("internalField")
        posz = txtzsplit.index("internalField")
        if full:
            posp = txtpsplit.index("internalField")
            posk = txtksplit.index("internalField")
            poso = txtosplit.index("internalField")
            posn = txtnsplit.index("internalField")
        itemsfrom = int(txtusplit[posu+3])
        # first item is posu+5

        us = []
        xs = []
        ys = []
        zs = []
        if full:
            ps = []
            ks = []
            os = []
            ns = []
        print('start reading U field', datetime.datetime.now())
        for i in range(itemsfrom):
            us.append([float(remove_left(txtusplit[posu+5+i*3])), float(txtusplit[posu+6+i*3]), float(remove_right(txtusplit[posu+7+i*3]))])
            xs.append(float(txtxsplit[posx + 5 + i]))
            ys.append(float(txtysplit[posy + 5 + i]))
            zs.append(float(txtzsplit[posz + 5 + i]))
            if full:
                ps.append(float(txtpsplit[posp + 5 + i]))
                ks.append(float(txtksplit[posk + 5 + i]))
                os.append(float(txtosplit[poso + 5 + i]))
                ns.append(float(txtnsplit[posn + 5 + i]))

        xs = np.asarray(xs)
        ys = np.asarray(ys)
        zs = np.asarray(zs)
        if full:
            ps = np.asarray(ps)
            ks = np.asarray(ks)
            os = np.asarray(os)
            ns = np.asarray(ns)

        xmin = xs.min()
        ymin = ys.min()
        zmin = zs.min()
        xmax = xs.max()
        ymax = ys.max()
        zmax = zs.max()
        grid_points = 27000000
        meters = (xmax-xmin) * (ymax - ymin) * (zmax - zmin)
    
        resolution = (meters / grid_points)**(1./3)
        if resolution==0:
            print('no clip boundaries in one of the axis')
        # resolution = 0.1
        ticksx = int((xmax-xmin) / resolution) * 1j
        ticksy = int((ymax-ymin) / resolution) * 1j
        ticksz = int((zmax-zmin) / resolution) * 1j
        print ('ticks',ticksx,ticksy,ticksz,ticksx*ticksy*ticksz, meters, resolution)

        xi, yi, zi =  np.mgrid[xmin:xmax:ticksx, ymin:ymax:ticksy, zmin:zmax:ticksz]

        print('finish reading U field', datetime.datetime.now())
        gridux = griddata((xs, ys, zs), np.asarray(us)[:,0], (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
        griduy = griddata((xs, ys, zs), np.asarray(us)[:,1], (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
        griduz = griddata((xs, ys, zs), np.asarray(us)[:,2], (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
        if full:
            gridp = griddata((xs, ys, zs), np.asarray(ps), (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
            gridk = griddata((xs, ys, zs), np.asarray(ks), (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
            grido = griddata((xs, ys, zs), np.asarray(os), (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings
            gridn = griddata((xs, ys, zs), np.asarray(ns), (xi, yi, zi), method='nearest')  # Nearest for keeping zeros that are buildings

#        plt.figure()
#        testxnp = gridu.ravel()
#        sc = plt.scatter(gridu.ravel(), gridux.ravel(), c=(gridu.ravel()-gridux.ravel()), s=1) #, cmap=cm)
#        plt.colorbar(sc)
#        plt.xlabel('gridu ordered file')
#        plt.ylabel('grid vtk')
#        plt.title('READING file in order ? ')
#        plt.show()
            
#        print('corr for reading is:', str(np.corrcoef(gridu.ravel(), gridux.ravel())[0, 1]),str(sklearn.metrics.r2_score(gridu.ravel(), gridux.ravel())))
        
    if 5==5:
        # write the new interpolate data to new directory
        #first we need to read Cx,Cy, Cz
        # than we need to interpolate and save
        
        directory = tosim +"0/"
        filex = directory + "Cx"
        filey = directory + "Cy"
        filez = directory + "Cz"
        fhx = open(filex, 'r')
        fhy = open(filey, 'r')
        fhz = open(filez, 'r')
        txtx = fhx.read()
        txty = fhy.read()
        txtz = fhz.read()
        txtxsplit = txtx.split()
        txtysplit = txty.split()
        txtzsplit = txtz.split()
        posx = txtxsplit.index("internalField")
        posy = txtysplit.index("internalField")
        posz = txtzsplit.index("internalField")
        itemsto = int(txtxsplit[posx+3])
        # first item is posu+5
        xs2 = []
        ys2 = []
        zs2 = []
        print('start reading C field', datetime.datetime.now())
        for i in range(itemsto):
            xs2.append(float(txtxsplit[posx + 5 + i]))
            ys2.append(float(txtysplit[posy + 5 + i]))
            zs2.append(float(txtzsplit[posz + 5 + i]))
        xs2 = np.asarray(xs2)
        ys2 = np.asarray(ys2)
        zs2 = np.asarray(zs2)
        
        xs2unique = np.unique(xs2)
        ys2unique = np.unique(ys2)
        zs2unique = np.unique(zs2)
        
        dx = xs2unique[1] - xs2unique[0]
        dy = ys2unique[1] - ys2unique[0]

################################        
        
        casedir = tosim+'0/U'
        fhur = open(casedir, 'r')
#            fhur = open(r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor/2000/U', 'r')
        utxt = fhur.readlines()
        posu = utxt.index('internalField   nonuniform List<vector> \n')
        if full:
            fhpr = open(fromsim+str(int(itter))+'/p', 'r')
            fhkr = open(fromsim+str(int(itter))+'k', 'r')
            fhor = open(fromsim+str(int(itter))+'omega', 'r')
            fhnr = open(fromsim+str(int(itter))+'nut', 'r')
            fhpw = open(tosim+'0/p', 'w')
            fhkw = open(tosim+'0/k', 'w')
            fhow = open(tosim+'0/'+turbulence, 'w')
            fhnw = open(tosim+'0/nut', 'w')
            ptxt = fhpr.readlines()
            ktxt = fhkr.readlines()
            otxt = fhor.readlines()
            ntxt = fhnr.readlines()
            posp = ptxt.index('internalField   nonuniform List<scalar> \n')
            posk = ktxt.index('internalField   nonuniform List<scalar> \n')
            poso = otxt.index('internalField   nonuniform List<scalar> \n')
            posn = ntxt.index('internalField   nonuniform List<scalar> \n')
#        print('before items = ', posu)
#        print('before items2 = ', utxt[posu+1])
#        items = int(utxt[posu+1])
        print('items = ', itemsto)
        
        xiunique = np.unique(xi[:,0,0])
        yiunique = np.unique(yi[0,:,0])
        ziunique = np.unique(zi[0,0,:])
        for i in range(len(xs2)):  # todo items
            if i % 10000 == 0:
                print('remap',i, len(xs2))
            indexi1 = (np.abs(xiunique - xs2[i] + dx/2)).argmin()
            indexi2 = (np.abs(xiunique - xs2[i] - dx/2)).argmin()
            indexj1 = (np.abs(yiunique - ys2[i] + dy/2)).argmin()
            indexj2 = (np.abs(yiunique - ys2[i] - dy/2)).argmin()
            if i>0:
                vz1 = (zs2[i]+zs2[i-1])/2.
            else:
                vz1 = zs2[i]
            if i<len(zs2)-1:
                vz2 = (zs2[i]+zs2[i+1])/2.
            else:
                vz2 = zs2[i]
            indexk0 = (np.abs(ziunique - zs2[i])).argmin()
            if zs2[i]<ziunique[indexk0]:
                indexk1=indexk0
                indexk2=indexk0+1
            else:                
                indexk1=indexk0-1
                indexk2=indexk0
            uu0 =gridux[indexi1:indexi2,indexj1:indexj2,indexk1:indexk2][gridux[indexi1:indexi2,indexj1:indexj2,indexk1:indexk2]!=0.].mean().mean()
            vv0 =griduy[indexi1:indexi2,indexj1:indexj2,indexk1:indexk2][griduy[indexi1:indexi2,indexj1:indexj2,indexk1:indexk2]!=0.].mean().mean()
            ww0 =griduz[indexi1:indexi2,indexj1:indexj2,indexk1:indexk2][griduz[indexi1:indexi2,indexj1:indexj2,indexk1:indexk2]!=0.].mean().mean()
            if np.isnan(uu0):
                uu0=0.0
            if np.isnan(vv0):
                vv0=0.0
            if np.isnan(ww0):
                ww0=0.0
            u0 = str(uu0)
            u1 = str(vv0)
            u2 = str(ww0)
            utxt[posu+3+i]='('+u0+' '+u1+' '+u2+')\n'
#            if full:
#                p = str(labelspredictp[i])
#                k = str(labelspredictk[i])
#                o = str(labelspredicto[i])
#                n = str(labelspredictn[i])
#                ptxt[posp + 3 + i] = p + '\n'
#                ktxt[posk + 3 + i] = k + '\n'
#                otxt[poso + 3 + i] = o + '\n'
#                ntxt[posn + 3 + i] = n + '\n'
        fhuw = open(tosim + '/0/U', 'w')            
        fhuw.writelines(utxt)
        fhuw.close()
#        if full:
#            fhpw.writelines(ptxt)
#            fhkw.writelines(ktxt)
#            fhow.writelines(otxt)
#            fhnw.writelines(ntxt)
#            fhpw.close()
#            fhkw.close()
#            fhow.close()
#            fhnw.close()
        print('fin new U, p, k, omega, nut', datetime.datetime.now())

        
if __name__ == "__main__":


#    textboxfile1value = r'/ibdata2/nirb/openFOAM/vortex/tlvs2/tlvs2.foam'
#    textboxtime1value = 999999999
#    db = get_clip(textboxfile1value, [textboxtime1value],
#              clipxmin=float(0), clipxmax=float(0),
#              clipymin=float(0), clipymax=float(0),
#              clipzmin=float(0), clipzmax=float(0))



    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())

#    from shapely.geometry import Polygon, LineString
#    import geopandas as gpd
#    line1 = LineString([(0.0, 0.0), (1.0, 0.0), (1., 1.), (0., 1.0),(0.0, 0.0)])
#    line11 = LineString([(-1., -1.), (-2, -2)])
#    line2 = LineString([(-0.5, 0.6), (1.7, 0.6)])
#    print(line1.intersectssssssssssssssssssssssssssion(line2))
#    res = line1.intersection(line2)
#    res2 = line1.intersection(line11)
#    print(res.is_empty)
#    print(res2.is_empty)

plt.close('all')   
filename=r'/ibdata2/nirb/openFOAM/vortex/tlv/tlv.foam'
corrlast, corrnear, ctimes = corrall(filename, 'U_x', 1 , 664000, 0, 0, 0, 0, 0, 0, False, end=74)
filename=r'/ibdata2/nirb/openFOAM/vortex/tlv21/tlv2.foam'
corrlast, corrnear, ctimes = corrall(filename, 'U_x', 1 , 664000, 0, 0, 0, 0, 0, 0, False)
filename=r'/ibdata2/nirb/openFOAM/vortex/caseB2/caseB1.foam'
corrlast, corrnear, ctimes = corrall(filename, 'U_x', 1 , 0, 0, 0, 0, 0, 0, 0, False, end=5)

################# start eddie
itter=1000000000
files=[
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomega/windAroundCube.foam',
      r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz50/windAroundCube.foam',
      r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz35/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegafloor/windAroundCube.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegafloorz20x5grid2/windAroundCube.foam',
# Check vs 40 loco size 
#      r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz40/windAroundCube.foam'
#      r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz20x9/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz20x5grid2/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundtest/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegarule/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegagrid2/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegagrid3/michaelstadt3.foam'
       ]        
axisindex = 1
axispos = 100000
xmin=300000
xmax=600000
ymin=0
ymax=0
zmin=0
zmax=35000
currentRow=1
    
eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow)
plot_file([files[1]], ['U_x'], itter, xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, '')
corrlast, corrnear, ctimes = corrall(files[1], 'U_x', 1 , axispos, xmin, xmax, ymin, ymax, zmin, zmax, False, end=744)
diff_files(files, [itter], 'U_x', xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, currentRow)


itter=1000000000
files=[
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomega/windAroundCube.foam',
       r'/data3/nirb/Simulations/Vortex/windAroundurbanMichelstadt2zomegafloorprofile/windAroundCube.foam',
#       r'/data3/nirb/Simulations/Vortex/windAroundurbanMichelstadt2zomegafloorprofilez20/windAroundCube.foam'
       r'/data3/nirb/Simulations/Vortex/windAroundurbanMichelstadt2zomegafloorprofilez20x20grid2/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegafloorz20x5grid2/windAroundCube.foam',
# Check vs 40 loco size 
#      r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz20/windAroundCube.foam'
#      r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz20x9/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz20x5grid2/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundtest/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegarule/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegagrid2/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegagrid3/michaelstadt3.foam'
       ]        
axisindex = 1
axispos = 100
xmin=00000
xmax=00000
ymin=0
ymax=0
zmin=0
zmax=70
currentRow=1
    
#eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow)
#plot_file([files[0]], ['U_x'], itter, xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, '')
corrlast, corrnear, ctimes = corrall(files[1], 'U_x', 1 , axispos, xmin, xmax, ymin, ymax, zmin, zmax, False, end=35)
diff_files(files, [itter], 'U_x', xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, currentRow)


itter=1000000000
files=[
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomega/windAroundCube.foam',
       r'/data3/nirb/Simulations/Vortex/windAroundurbanMichelstadt2zomegafloorprofile/windAroundCube.foam',
       r'/data3/nirb/Simulations/Vortex/windAroundurbanMichelstadt2zomegafloorprofilez20x20grid2/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegafloorz20x5grid2/windAroundCube.foam',
# Check vs 40 loco size 
#      r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz20/windAroundCube.foam'
#      r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz20x9/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegaz20x5grid2/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundtest/windAroundCube.foam'
#       '/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegarule/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegagrid2/windAroundCube.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/windAroundurbanMichelstadt2zomegagrid3/michaelstadt3.foam'
       ]        
axisindex = 1
axispos = 100
xmin=-750
xmax=750
xmin=375 #-250
xmax=520 # 500 520
ymin=0
ymax=0
zmin=0
zmax=35
currentRow=1
plt.close('all')   
 
eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow)
plot_file([files[0]], ['U_x'], itter, xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, '')
corrlast, corrnear, ctimes = corrall(files[1], 'U_x', 1 , axispos, xmin, xmax, ymin, ymax, zmin, zmax, False, end=744)
diff_files(files, [itter], 'U_x', xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, currentRow)


itter=1000000000
files=[r'/ibdata2/nirb/openFOAM/ml/windAroundurbanMichelstadt2zomegablock1/windAroundCube.foam']
axisindex = 1
axispos = -250000#100000
xmin=0#300000
xmax=0#300000
ymin=0
ymax=0
zmin=0
zmax=0
currentRow=1
    
eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow)



itter=40000
files=[
       r'/data4bk/nirb/Simulations/Dans/tlv0/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/tlv3/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/tlv2/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/tlvdans3/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/tlvbig0/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/tlvbig3/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/tlvbig2/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/tlvdans3big/tlvs2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlv/tlv.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlv3/tlv3.foam'
       r'/ibdata2/nirb/openFOAM/vortex/tlv1/tlv1.foam'

#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2z20x80/tlvs2.foam'
       ]        
axisindex = 1
axispos = 664200
xmin=178750
xmax=179000
ymin=0
ymax=0
zmin=0
zmax=40
currentRow=1
    
eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow)

diff_files(files, [itter], 'U_x', xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, currentRow)
diff_files(files, [itter], 'U_z', xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, currentRow)

corrlast, corrnear, ctimes = corrall(files[4], 'U_x', 1 , axispos, xmin, xmax, ymin, ymax, zmin, zmax, False, end=2774)
plot_file([files[1]], ['U_x'], itter, xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, '')


itter=40000
files=[
#       r'/ibdata2/nirb/openFOAM/vortex/tlv/tlv.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs/tlvs.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs1/tlvs1.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2/tlvs2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2g/tlvs2.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2shm2/tlvs2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2shm/tlvs2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2dx/tlvs2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2dx2/tlvs2.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2z20/tlvs2.foam',
       r'/ibdata2/nirb/openFOAM/vortex/tlvs2z20x80b/tlvs2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2z20x20b/tlvs2.foam'
       r'/ibdata2/nirb/openFOAM/vortex/tlvrule/tlvs2.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2z40/tlvs2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlvs2z50/tlvs2.foam'
#       r'/ibdata2/nirb/openFOAM/vortex/tlv3/tlv3.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlv2/tlv2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/tlv1/tlv1.foam'
       ]        
xmin=177750
xmax=178000
bounds=[178750,17900,664000,664250]
vprofile(files, bounds, itter)


itter=1000000000
files=[
       r'/ibdata2/nirb/openFOAM/vortex/caseB1/caseB1.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/caseB2/caseB2.foam',
#       r'/ibdata2/nirb/openFOAM/vortex/caseB3/caseB3.foam'
       r'/ibdata2/nirb/openFOAM/vortex/caseB4/caseB4.foam'
       ]        
axisindex = 1
axispos = 0
xmin=0
xmax=0
ymin=0
ymax=0
zmin=0
zmax=0
currentRow=1
    
eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow)

################ end eddie


import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(1,100,101)
y = np.linspace(1,100,101)
xx,yy = np.meshgrid(x,y)
size = 396
w=-np.sin(xx/size)  
u=np.sin(yy/size) #eddies
#u=xx/xx #waves
scl=100
plt.rcParams["figure.figsize"] = (20,10)

from scipy import integrate
intx=integrate.cumtrapz(w,xx,axis=1,initial=0)[0]
inty=integrate.cumtrapz(u,yy,axis=0,initial=0)
psi1=intx-inty               
plt.figure()
plt.subplot(1,2,1)
plt.imshow(psi1,origin='lower')
plt.colorbar()
plt.quiver(xx,yy,u,w, scale=scl)
plt.title('u')
plt.subplot(1,2,2)
plt.imshow(w,origin='lower')
plt.colorbar()
plt.quiver(xx,yy,u,w, scale=scl)
plt.title('w')
plt.show()      

   # integrate to make an intial guess
intx=integrate.cumtrapz(w,xx,axis=1,initial=0)[0]
inty=integrate.cumtrapz(u,yy,axis=0,initial=0)
psi1=intx-inty
plt.figure()
plt.imshow(psi1,origin='lower')
plt.title('psi1')
print('shae intx',intx.shape)
print('shae inty',inty.shape)

intx=integrate.cumtrapz(w,xx,axis=1,initial=0)
inty=integrate.cumtrapz(u,yy,axis=0,initial=0)[:,0][:,None]
psi2=intx-inty
plt.figure()
plt.imshow(psi2,origin='lower')
plt.title('psi2')

psi=0.5*(psi1+psi2)

intx=integrate.cumtrapz(u,xx,axis=1,initial=0)[0]
inty=integrate.cumtrapz(w,yy,axis=0,initial=0)
chi1=intx+inty

intx=integrate.cumtrapz(u,xx,axis=1,initial=0)
inty=integrate.cumtrapz(w,yy,axis=0,initial=0)[:,0][:,None]
chi2=intx+inty

plt.figure()
plt.imshow(psi1+psi2,origin='lower')
plt.title('psi1+2')
   
plt.figure()
plt.imshow(chi2+chi1,origin='lower')
plt.title('chi2+1')

itter=40000
files=[
       r'/data4bk/nirb/Simulations/Dans/natanya0/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/natanya1/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/natanya2/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/natanya3/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/natanyadans/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/natanyadans1/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/natanyadans0a/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/natanyadans0b/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/natanyadans0e/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/natanyadans0b1/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/natanyadans0e1/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/natanyadans111111111111111111111111111110f/tlvs2.foam'
       ]     
areaname="natanya"
axisindex = 1
axispos = 693000
xmin=184000
xmax=192000
xmin=186000 #184000
xmax=188000 #192000
ymin=0
ymax=0
zmin=0
zmax=100
currentRow=1
startarea=0
endarea=999999
currentitem =  u'U_x'
profile4(files,itter,axisindex,axispos,'',areaname,startarea,endarea,currentitem,0,0,0,0,0,0)
corrlast, corrnear, ctimes = corrall(files[9], 'U_x', 1 , axispos, xmin, xmax, ymin, ymax, zmin, zmax, False, end=24)
diff_files([r'/data4bk/nirb/Simulations/Dans/natanyao2/tlvs2.foam', r'/data4bk/nirb/Simulations/Dans/natanya2/tlvs2.foam'], [itter], 'U_x', xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, currentRow)


itter=40000
files=[
       r'/data4bk/nirb/Simulations/Dans/ashkelon0/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/ashkelon1/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/ashkelon2/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/ashkelondans/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/ashkelondans1/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/ashkelondans0a/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/ashkelondans0b/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/ashkelondans0e/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/ashkelondans0b1/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/ashkelondans0e1/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/ashkelondans0f/tlvs2.foam'
       ]      
areaname='ashkelon'
axisindex = 1
axispos = 620000
xmin=156000
xmax=164000
ymin=0
ymax=0
zmin=0
zmax=100
currentRow=1
startarea=0
endarea=9999
currentitem =  u'U_x'
corrlast, corrnear, ctimes = corrall(files[9], 'U_x', 1 , axispos, xmin, xmax, ymin, ymax, zmin, zmax, False, end=20)
profile4(files,itter,axisindex,axispos,'',areaname,startarea,endarea,currentitem,0,0,0,0,0,0)
plot_file([files[0]], ['U_x'], itter, xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, '')

itter=40000
files=[
#       r'/data4bk/nirb/Simulations/Dans/bs0/tlvs2.foam',      
       r'/data4bk/nirb/Simulations/Dans/bs0a/tlvs2.foam',      
#       r'/data4bk/nirb/Simulations/Dans/bs1/tlvs2.foam',      
#       r'/data4bk/nirb/Simulations/Dans/bs1a/tlvs2.foam',      
#       r'/data4bk/nirb/Simulations/Dans/bs2/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/bs2a/tlvs2.foam',
       
#       r'/data4bk/nirb/Simulations/Dans/bsdans/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/bsdans1/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/bsdans0a/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/bsdans0b/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/bsdans0e/tlvs2.foam',
#       r'/data4bk/nirb/Simulations/Dans/bsdans0f/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/bsdans0a1/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/bsdans0b1/tlvs2.foam',
       r'/data4bk/nirb/Simulations/Dans/bsdans0e1/tlvs2.foam'

       ]       
areaname='bs'
axisindex = 1
axispos = 574000
xmin=181000
xmax=182000
xmin=0
xmax=0
ymin=0
ymax=0
zmin=0
zmax=1000
currentRow=1
startarea=0 #586
endarea=41000000 #587
currentitem =  u'U_x'
xmingui=xmin
ymingui=ymin
zmingui=zmin
xmaxgui=xmax
ymaxgui=ymax
zmaxgui=zmax
verbose=True
learnfrom=''
corrlast, corrnear, ctimes = corrall(files[2], 'U_x', 1 , axispos, xmin, xmax, ymin, ymax, zmin, zmax, False, end=27)
profile4(files,itter,axisindex,axispos,'',areaname,startarea,endarea,currentitem,0,0,0,0,0,0)
plot_file([files[0],files[1]], ['U_x'], itter, xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, '')
diff_files(files, [itter], 'U_x', xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, currentRow)
#profile3(files,itter,axisindex,axispos,'',areaname,40,41,currentitem,0,0,0,0,0,0, verbose=True)
diff_files([r'/data4bk/nirb/Simulations/Dans/bso2/tlvs2.foam', r'/data4bk/nirb/Simulations/Dans/bs2/tlvs2.foam'], [itter], 'U_x', xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, currentRow)
profile4(files,itter,axisindex,axispos,'',areaname,183,184,currentitem,0,0,0,0,0,0,verbose=True)



itter=40000
files = [
#        r'/ibdata2/nirb/openFOAM/vortex/tlv/tlv.foam',
#        r'/ibdata2/nirb/openFOAM/vortex/tlv0/tlv.foam',
                r'/data4bk/nirb/Simulations/Dans/tlva/tlv.foam',
                r'/data4bk/nirb/Simulations/Dans/tlvsmall2/tlv.foam',
#        r'/data4bk/nirb/Simulations/Dans/tlvsmall3/tlv.foam',
                 r'/data4bk/nirb/Simulations/Dans/tlvsmall0b2/tlv.foam',
                 r'/data4bk/nirb/Simulations/Dans/tlvsmall0f2/tlv.foam'

#         u'/data4bk/nirb/Simulations/Dans/tlvbig0/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbig1/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbig2/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigmap/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigmap2/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans1/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans1a/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans1b/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans0a/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans0d1/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans0e1/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans0b1/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans0c1/tlvs2.foam',
#         u'/data4bk/nirb/Simulations/Dans/tlvbigdans0f/tlvs2.foam'

]
areaname='tlvbig250'
textboxtime1value = 300000.0
axisindex = 1
axispos = 664200.0
learnfrom = ''
startarea =  0
endarea = 106000
currentitem =  u'U_x'
xmingui = 0.0
ymingui = 0.0
zmingui = 0.0
xmaxgui = 0.0
ymaxgui = 0.0
zmaxgui = 0.0
currentRow=1

#profile3(files,textboxtime1value,axisindex,axispos,learnfrom,areaname,startarea,endarea,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui)
corrlast, corrnear, ctimes = corrall(files[0], 'U_x', 1 , axispos, xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, False, end=270, last=7)
plot_file([files[0]], ['U_x'], itter, xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, axisindex, axispos, '')
diff_files([files[0], files[1]], [itter], 'U_x', xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, axisindex, axispos, currentRow)
plot_file([files[0]], ['U_x'], itter, xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, 2, 100, '')

#600 601
profile4(files,textboxtime1value,axisindex,axispos,learnfrom,areaname,0,600001,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui, verbose=False)

fromsim= u'/data4bk/nirb/Simulations/Dans/tlvbig2/'
tosim = u'/data4bk/nirb/Simulations/Dans/tlvbigmap2/'
remapsimulation(fromsim, tosim)


plot_file(r'/ibdata2/nirb/openFOAM/windAroundurbanMichelstadt2zomega1/windAroundCube.foam', ['U_x'], 9999999, 0, 0, 0, 0, 0, 0, 2, 10000, '', parallel=False)

files=[r'/ibdata2/nirb/openFOAM/porous/sinthetic0a/log',
       r'/ibdata2/nirb/openFOAM/porous/sinthetic0s/log']

plot_file(files, ['U_x'], 9999999, 0, 0, 0, 0, 0, 0, 2, 50, '', parallel=False)



itter=40000
files = [ r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor1/windAroundcaseE.foam',
         r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor1ml/windAroundcaseE.foam'
]
areaname='ml'
textboxtime1value = 300000.0
axisindex = 1
axispos = -200.0
learnfrom = ''
currentitem =  u'U_x'
xmingui = 0.0
ymingui = 0.0
zmingui = 0.0
xmaxgui = 0.0
ymaxgui = 0.0
zmaxgui = 0.0
currentRow=1
xmin=625
xmax=750
ymin=0
ymax=0
zmin=0
zmax=20

#profile3(files,textboxtime1value,axisindex,axispos,learnfrom,areaname,startarea,endarea,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui)
corrlast, corrnear, ctimes = corrall(files[0], 'U_x', 1 , axispos, xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, False, end=300, last=400)
eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow,parallel=False)
diff_files([files[0], files[1]], [textboxtime1value], 'U_x', xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, axisindex, axispos, currentRow)
diff_files([files[0], files[1]], [textboxtime1value], 'U_x', xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, 2, 10, 2, parallel=False)
plot_file([files[0]], ['U_z'], itter, xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, '', parallel=False)

ml(files)
diff_files([files[0], files[1]], [textboxtime1value], 'U_z', 0, 0, -0, 150, 0, zmaxgui, 2, 10, 2, parallel=False)
diff_files([files[0], files[1]], [textboxtime1value], 'U_x', -400, 400, -400, 400, zmingui, zmaxgui, 2, 10, 2, parallel=False)
diff_files([files[0], files[1]], [textboxtime1value], 'U_z', -500, 500, 0, 0, 0, 30, 1, 10, 1, parallel=False)


def quad(x):
    # typical drone speed - up to 20m/s there are few that can go 80m/s
    #typical wind speed at 1000m a.g.l. - 10 m/s
    
    # 2add - drone speed direction
    # 2add - choose best wind level (wind goes logaritmic till certain height)
    # 2add - wind direction can change with height - take from GFS
    # check if there is a patent for boats, drones
    # use slingshot with towers (like with spaceship) 
    # agricalture and delivery (transportation)
    # can calculate maximum distance according to wind
    # unmanned aerial vehicle (UAV) 
    # 
    import numpy as np
    import math
    mul=10.  # multiplication of the calculation for smoothing
    windangle = 90 # 0 is in the drone direction
    logprofile = 1.*np.log((np.linspace(50,1000,20)-20)/1.)
    wind = 2/mul
#    wind = logprofile[1]   
    wind = logprofile[1]    
    windx = wind/mul*math.cos(windangle/180.*math.pi)
    windy = wind/mul*math.sin(windangle/180.*math.pi)
#    speed = 15/mul*math.cos(windangle)
    speed = 10/mul
    upspeed = 5/mul
    startx=0
    starty=0
    endx=40000
    endy=0
    far = (endx-startx)**2.+(endy-starty)**2.

    option = 1
    if option==1:  
        posx=startx
        posy=starty
        fuel=0
        far1 = 0
        print ('option edvance start @ ',startx,starty)
        angle = math.asin(-(windy)/(speed))
        # angle = math.asin(((endy-starty)/(endx-startx)*(sx+windx)-windy)/(speed))
        while posx<endx:
            posx+= windx + speed * math.cos(angle)
            posy+= windy + speed * math.sin(angle)
            far1+=1
        print ('pos ', far1, ' @ ',posx,posy,angle/math.pi*180.)
            
    
    if option==2:
        posx=startx
        posy=starty
        fuel=0
        far2 = 0
    #    while far<100:
        print ('option corrections start @ ',startx,starty)
        while posx<endx:
            angle = math.asin((posy-endy)/(posx-endx))
            posx+=windx + speed * math.cos(angle) 
            posy+=windy + speed * math.sin(angle)
            far2+=1
        print ('pos ', far2, ' @ ',posx,posy,angle/math.pi*180.)

    print ('ratio ', 1.-float(far1)/far2, 1800.*(1.-float(far1)/far2))
            
    
    
    return

files=[u'/ibdata2/nirb/openFOAM/porous/windAroundBuildings-davidsona/windAroundBuildings-davidsona.OpenFOAM',
       u'/ibdata2/nirb/openFOAM/porous/windAroundBuildings-davidson0/windAroundCube.foam']
axisindex = 2
axispos = 1.15
xminc=0
xmaxc=0
yminc=0
ymaxc=0
zminc=0
zmaxc=0
currentRow=1  
itter=444444
fields=['U_x']
filesindex=1
parallel=False
f=0

textboxfile1value = files[filesindex] #self.filename1.text()
db = get_slice(textboxfile1value, itter, axisindex, axispos,
               clipxmin=float(xminc), clipxmax=float(xmaxc),
               clipymin=float(yminc), clipymax=float(ymaxc),
               clipzmin=float(zminc), clipzmax=float(zmaxc), parallel=parallel)

plt.figure()
for f in range(len(fields)):
    xmin, xmax, ymin, ymax, zi1 = makegrid(db, fields[f], axisindex, ticks=5000, method='linear')
    zi1[zi1 == 0] = np.NaN
    plt.imshow(zi1, origin='lower')  # jet, Paired
    plt.title(shrinktitle(textboxfile1value)+'-'+fields[f])
    plt.colorbar()
plt.show()

plt.figure()
plt.plot(np.mean(zi1,axis=0)/zi1[2500,0])
plt.show()

diff_files(files, [itter], 'U_x', xminc, xmaxc, yminc, ymaxc, zminc, zmaxc, axisindex, axispos, currentRow)


plot_file([files[0]], ['U_x'], 9999999, 0, 0, 0, 0, 0, 0, 2, 1.150000, '', parallel=True)




itter=40000
files = [ r'/data4bk/nirb/Simulations/Mala/mala2b/mala2b.foam',
]
areaname='mala'
textboxtime1value = 300000.0
axisindex = 1
axispos = 563735.0
learnfrom = ''
currentitem =  u'U_x'
xmingui = 0.0
ymingui = 0.0
zmingui = 0.0
xmaxgui = 0.0
ymaxgui = 0.0
zmaxgui = 0.0
currentRow=1
xmin=155200
xmax=155400
ymin=0
ymax=0
zmin=0
zmax=0

#profile3(files,textboxtime1value,axisindex,axispos,learnfrom,areaname,startarea,endarea,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui)
corrlast, corrnear, ctimes = corrall(files[0], 'U_x', 1 , axispos, xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, False, end=300, last=400, parallel=True)
eddies(files, itter, axisindex, axispos, xmin, xmax, ymin, ymax, zmin, zmax, currentRow,parallel=False)
diff_files([files[0], files[1]], [textboxtime1value], 'U_x', xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, axisindex, axispos, currentRow)
plot_file([files[0]], ['U_z'], itter, xmin, xmax, ymin, ymax, zmin, zmax, axisindex, axispos, '', parallel=False)
profile4(files,textboxtime1value,axisindex,axispos,learnfrom,areaname,0,1,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui, verbose=False, parallel=True)

    

itter=40000
files = [ r'/data4bk/nirb/Simulations/michaelstadtfloor1ml/windAroundcaseE.foam',
          r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor1/windAroundcaseE.foam'
]
files = [ r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor1/windAroundcaseE.foam',
          r'/ibdata2/nirb/openFOAM/ml/michaelstadtfloor1ml/windAroundcaseE.foam'
]
areaname='michaelstadt'
textboxtime1value = 300000.0
axisindex = 0# 1
axispos = 0 #-200.0
learnfrom = ''
currentitem =  u'U_x'
xmingui = -400.0
ymingui = -300.0
zmingui = 0.0
xmaxgui = 400.0
ymaxgui = 300.0
zmaxgui = 60.0
currentRow=0 #1
xmin=625
xmax=750
xmin=0
xmax=0
ymin=0
ymax=0
zmin=0
zmax=20
zmax=0
profile4(files,textboxtime1value,axisindex,axispos,learnfrom,areaname,0,1,currentitem,xmingui,xmaxgui,ymingui,ymaxgui,zmingui,zmaxgui, verbose=True, parallel=True)
corrlast, corrnear, ctimes = corrall(files[0], 'U_x', 1 , axispos, xmingui, xmaxgui, ymingui, ymaxgui, zmingui, zmaxgui, False, end=3, last=4, parallel=True)
