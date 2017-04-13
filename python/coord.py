# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 12:00:16 2016

@author: Super
"""
import numpy as np

def nu2uv(lat,lon,hr,dec):
#   UNTITLED Summary of this function goes here
#   Detailed explanation goes here
    nu2ecf = [ [-np.sin(lat) * np.cos(lon), -np.sin(lon), np.cos(lat) * np.cos(lon)],\
                       [-np.sin(lat) * np.sin(lon), np.cos(lon), np.cos(lat) * np.sin(lon)],\
                   [np.cos(lat), 0, np.sin(lat)] ]
    ecf2uv = [ [np.sin(hr), np.cos(hr), 0],\
                   [-np.sin(dec) * np.cos(hr), np.sin(dec) * np.sin(hr), np.cos(dec)],\
                   [np.cos(dec) * np.cos(hr), -np.cos(dec) * np.sin(hr), np.sin(dec)] ]
    ma = np.dot(ecf2uv, nu2ecf)
    return ma