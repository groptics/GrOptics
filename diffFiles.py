#!/usr/bin/env python

import os
import sys
import filecmp
import shutil

makecopy = "true"
#makecopy = "false"

diffDevCmdBase = "../../GrOptics"
diffDevConfig = diffDevCmdBase + "/Config"
diffDevInc = diffDevCmdBase + "/include"
diffDevSrc = diffDevCmdBase + "/src"

fConfig = filecmp.dircmp("./Config",diffDevConfig)
fBase = filecmp.dircmp("./",diffDevCmdBase)
fInc =  filecmp.dircmp("./include",diffDevInc)
fSrc =  filecmp.dircmp("./src",diffDevSrc)

for name in fBase.diff_files:
    namedir = diffDevConfig + "/" + name
    print namedir
for name in fConfig.diff_files:
    print name
for name in fInc.diff_files:
    print name
for name in fSrc.diff_files:
    print name

testDir = "./Copy/"

for name in fBase.diff_files:
    namedir = diffDevCmdBase + "/" + name
    namefi  = "./" + name
    print namedir, "  ", namefi
    if makecopy == "true":
        shutil.copyfile(namedir,namefi)

for name in fConfig.diff_files:
    namedir = diffDevConfig + "/" + name
    namefi  = "./Config/" + name
    print namedir,"  ",namefi 
    if makecopy == "true":
        shutil.copyfile(namedir,namefi)
    
for name in fInc.diff_files:
    namedir = diffDevInc + "/" + name
    namefi  = "./include/" + name
    print namedir,"  ",namefi 
    if makecopy == "true":
        shutil.copyfile(namedir,namefi)
    
for name in fSrc.diff_files:
    namedir = diffDevSrc + "/" + name
    namefi  = "./src/" + name
    print namedir,"  ",namefi 
    if makecopy == "true":
        shutil.copyfile(namedir,namefi)
    




#print fcmp.diff_files
configFi = ["arrayConfig.cfg", \
"geoConfig.cfg", \
"headerForVeritasConfigFile", \
"opticsSimulation.pilot", \
"photon.cph", \
'stdSCTelescopes.cfg', \
'stdSegSCTelescopes.cfg', \
'veritas.cfg']

#print configFi[0]
