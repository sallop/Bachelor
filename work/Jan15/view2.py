#!/usr/bin/python

import os

def FileData(fin, fout,
             title, xlabel, ylabel, option):
    return {'fin'   : fin   ,
            'fout'  : fout  ,
            'title' : title ,
            'xlabel': xlabel,
            'ylabel': ylabel,
            'option': option }

files = [ FileData("r3zv.dat", "r3zv_zz.gif"   , "zz"     , "[sec]", "[m]"      , "u 1:2 w l"),
          FileData("r3zv.dat", "r3zv_zv.gif"   , "zv"     , "[sec]", "[m/sec]"  , "u 1:3 w l"),
          FileData("r3zv.dat", "r3zv_za.gif"   , "za"     , "[sec]", "[m/sec^2]", "u 1:4 w l"),
          FileData("r3zv.dat", "r3zv_v1x.gif"  , "v1x"    , "[sec]", "[m]"      , "u 1:5 w l"),
          FileData("r3zv.dat", "r3zv_v1xv.gif" , "v1xv"   , "[sec]", "[m/sec]"  , "u 1:6 w l"),
          FileData("r3zv.dat", "r3zv_v1xa.gif" , "v1xa"   , "[sec]", "[m/sec^2]", "u 1:7 w l"),
          FileData("r3zv.dat", "r3zv_r2_x2.gif", "r^2-x^2", "[sec]", "[m/sec^2]", "u 1:8 w l"),]

config = [ "set grid"         ,
           "set size 0.6, 0.6",
#           "set terminal gif" ,]
           ]

gp = os.popen("gnuplot", "w")
#gp = os.popen("gnuplot","w+")

for cmd in config:
    print >> gp, cmd

#print >> gp, "plot sin(x)\n"
#print >> gp, "pause 2.0\n"
#

for f in files:
#    print >> gp, "set output '%s'" % ("dir-zcurve/pderiv/"+f['fout'])
    print >> gp, "set output '%s'" % ("Foo/"+f['fout'])
    print >> gp, "set xlabel '%s'" % (f['xlabel'])
    print >> gp, "set ylabel '%s'" % (f['ylabel'])
    print >> gp, "set title  '%s'" % (f['title'])
    print >> gp, "plot '%s' %s title '%s'" % ("dir-etc/"+f['fin'],
                                              f['option']        ,
                                              f['title'] )
    print >> gp, "pause 2.0"
    gp.flush()


### def gifname(fname):
###     return fname+".gif"
### 
### def datname(fname):
###     return fname+".dat"
### 
### def plot(config, file, option):
###     gp = os.popen("gnuplot", "w")
### #gp = os.popen("gnuplot","w+")
###     for cmd in config:
###         print >> gp, cmd
### 
###     for f in files:
### #    print >> gp, "set output '%s'" % ("dir-zcurve/pderiv/"+f['fout'])
###         print >> gp, "set output '%s'" % ("Foo/"+f['fout'])
###         print >> gp, "set xlabel '%s'" % (f['xlabel'])
###         print >> gp, "set ylabel '%s'" % (f['ylabel'])
###         print >> gp, "set title  '%s'" % (f['title'])
###         print >> gp, "plot '%s' %s title '%s'" % ("dir-etc/"+f['fin'],
###                                                   f['option']        ,
###                                                   f['title'] )
###         print >> gp, "pause 2.0"
###         gp.flush()
