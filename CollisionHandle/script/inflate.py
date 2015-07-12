#! /usr/bin/env python
import os

inflate_app = "./script/inflator"
simp_mesh = "./data/bunny/model"
inflated_mesh = simp_mesh+"_infalted.obj"
os.system( inflate_app+" "+ simp_mesh +" "+inflated_mesh+" "+" 3 100 ")
