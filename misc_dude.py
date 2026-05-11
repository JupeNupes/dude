import os
import fabio
from readMDA import *
import numpy as np


MotorMNE = {\
"26idcnpi:m10.VAL": "fomx",\
"26idcnpi:m11.VAL": "fomy",\
"26idcnpi:m12.VAL": "fomz",\
"26idcnpi:m17.VAL": "samy",\
"26idcnpi:Y_HYBRID_SP.VAL": "hybridy",\
"26idcnpi:X_HYBRID_SP.VAL": "hybridx",\
"26idbATTO:m4.VAL": "attox",\
"26idbATTO:m3.VAL": "attoz",\
"26idbATTO:m2.VAL": "phi",\
"26idbATTO:m1.VAL": "chi",\
'26idbATTO:PIC867:1:m1.VAL': "samth",\
'26idcDET:base:Theta.VAL': "TwoTheta",\
'26idb:DAC1_1.VAL' : "DAC1",\
}


def Display2Data(Axe, x, y):
    """display coordinate to data coordinate"""

    return Axe.transData.inverted().transform(np.array([(x,y)]))[0]


def MotorMNE_Parser(longname):
    """ return the motor nmenitic mnemotic from the description string"""
    
    try:
        return MotorMNE[longname]
    except:
        return longname.split(":")[-1].split(".")[-2]


def image_loader(fname):

    return fabio.open(fname).data


def mda_loader(fname):

    # Use basename to avoid dots in directory names (e.g. smith.scott) breaking the split
    scannum = os.path.splitext(os.path.basename(fname))[0].split("_")[-1]
    try:
        data = readMDA(fname, verbose=0)
    except Exception as e:
        print(fname, e)
        return [str(int(scannum))]+["-"]*9+["","","",fname,False]
    ctime = -1
    ndim = len(data)-1 if data[0]["dimensions"][-1] != 2048 else len(data)-2
    try:
        if ndim == 1: # 1D Scan 
            for d in data[1].d:
                if d.name == "26idcXMAP:PresetReal":
                    ctime = d.data[0]
            return [str(int(scannum)),\
                    MotorMNE_Parser(data[1].p[0].name),\
                    "{0:.3f}".format(data[1].p[0].data[0]),\
                    "{0:.3f}".format(data[1].p[0].data[-1]),\
                    str(data[0]["dimensions"][0]),\
                    "", "", ""," ",\
                    "{0:.1f}".format(ctime), "", "", "", fname, False]
        elif ndim == 2: # 2D Scan
            for d in data[2].d:
                if d.name == "26idcXMAP:PresetReal":
                    ctime = d.data[0][0]
            return [str(int(scannum)),\
             MotorMNE_Parser(data[1].p[0].name),\
             "{0:.3f}".format(data[1].p[0].data[0]),\
             "{0:.3f}".format(data[1].p[0].data[-1]),\
             str(data[0]["dimensions"][0]),\
             MotorMNE_Parser(data[2].p[0].name),\
             "{0:.3f}".format(data[2].p[0].data[0][0]),\
             "{0:.3f}".format(data[2].p[0].data[-1][-1]),\
             str(data[0]["dimensions"][1]),\
            "{0:.1f}".format(ctime), "", "", "", fname, False]
    except Exception as e:
        print (e)
        return [str(int(scannum))]+["-"]*9+["","","",fname,False]

