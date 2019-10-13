# -*- coding: utf-8 -*-
"""
Specify what data you want to upload from the Inner Shelf Project 

written by Jack McSweeney (help from Zech Thurman) 
October 12, 2019 
"""

#def DataPull():
from scipy.io import loadmat
import numpy as np 

dir= '/Volumes/InnerShelf1/JackAnalysis/McSweeney_data/InnerShelf_AlongShoreVariability/'

fnames= np.array(['MS100','OC50','OC40N'])

#loadedmats={}
fext='.mat'

for moor in fnames:
    filepath=(dir+'Moorings/'+moor+fext)
    print(filepath)
#    loadedmats["{0}".format(moor)]=loadmat(filepath)
    globals()[moor] = loadmat(filepath)
    

    
#    print()
    
#    return 
