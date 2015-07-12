#! /usr/bin/env python
import os
from utility import *
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------functions--------------------------------------------

def grepInt(log_file,keyname,default=""):
    data = grepNumberWithKey(log_file,keyname)
    if len(data) <= 0:
        return default
    return str(int(data[0]))


def grepFloat(log_file,keyname,default=""):
    data = grepNumberWithKey(log_file,keyname)
    if len(data) <= 0:
        return default
    return str('%.5g'%(data[0]))


def grepStr(log_file,keyname,default=""):
    data = grepStrWithKey(log_file,keyname)
    if len(data) <= 0:
        return default
    return str(data[0])

def grepFloatList(log_file,keyname,default=""):
    data = grepNumberWithKey(log_file,keyname)
    if len(data) <= 0:
        return default
    d = ""
    for var in data:
        d = d +str('%.3g'%(var))+", "
    return d

def count(log_file,key):
    data_file = open(log_file)
    data_lines = data_file.readlines()
    num = 0
    for line in data_lines:
        if line.find(key) >= 0:
            num = num+1
    return num

def sum(log_file,key):
    s = 0
    data = grepNumberWithKey(log_file,key)
    for var in data:
        s = s + var
    return s

def sumInt(log_file,key):
    s = 0
    data = grepNumberWithKey(log_file,key)
    for var in data:
        s = s + int(var)
    return s

def averageInt(log_file,key):
    s = 0
    data = grepNumberWithKey(log_file,key)
    for var in data:
        s = s + int(var)
    n = 1;
    if len(data) > 0:
        n = len(data)
    return str(s/n)

def averageFloat(log_file,key):
    s = 0
    data = grepNumberWithKey(log_file,key)
    for var in data:
        s = s + float(var)
    n = 1;
    if len(data) > 0:
        n = len(data)
    return str( '%.3g'%(s/n) )

#-----------------------------main------------------------------------------------
log_f = "/home/simba/bunny_two.txt"
os.system("/home/simba/Workspace/CollisionHandle/bin/release/collision_handle ./two.ini  > "+log_f)

cg = sumInt(log_f,"cg steps: ")
prop = sumInt(log_f,"prop steps: ")
exp = sumInt(log_f,"exp steps: ")
cg_it = float(cg)/float(cg+prop+exp)
total_time = averageFloat(log_f,"total solving:")

print "t: "+str(total_time)
print "cg/it: "+str(cg_it)
