#! /usr/bin/env python

import os
import glob

def addSubFold(file_path):
    file_name = os.path.basename(file_path)
    file_dir = os.path.dirname(file_path)
    sub_fold = file_dir + "/" + file_name[0:5]
    if not os.path.exists(sub_fold):
        os.makedirs(sub_fold)
    return sub_fold + "/" + file_name

def changeName(old_name):
    new_name = old_name
    if old_name[-6] == '_':
        new_name = old_name[0:-5] + "00" + old_name[-5:len(old_name)]
    elif old_name[-7] == '_':
        new_name = old_name[0:-6] + "0" + old_name[-6:len(old_name)]
    return new_name

def reName(dir_name,ext = "*.obj"):
    files = glob.glob(dir_name+ext)
    for old_name in files:
        new_name = addSubFold(changeName(old_name))
        print new_name
        if new_name != old_name:
            os.system("mv " + old_name + " " + new_name)

reName("./data/dino/tempt_cubes_test/obj/")
# reName("./data/dragon/tempt_stair/obj/")
# reName("./data/bunny/tempt_two/obj/")
# reName("./data/longcube/tempt_ica/obj/")
# reName("./data/longcube/tempt_ica2/obj/")
