# -*- coding: utf-8 -*-
"""
Specify what data you want to upload from the Inner Shelf Project 

Written by Jack McSweeney (Tech Support: Zech Thurman)
October 12, 2019 
"""


import numpy as np
from scipy.io import loadmat


matdir = '/Volumes/InnerShelf1/JackAnalysis/McSweeney_data/InnerShelf_AlongShoreVariability/Moorings/'
moorings = np.array(['MS100', 'OS50', 'OC40N'])

for mooring in moorings:
    filepath = matdir + mooring
    mooring = loadmat(filepath, appendmat=True)


# def DataPull(matdir: str, entities: np.array):
#     for entity in entities:
#         filepath = matdir + entity
#         entity = loadmat(filepath, appendmat=True)
#     return
