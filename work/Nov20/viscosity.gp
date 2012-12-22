#
set size 0.6, 0.6
set grid

set key right top

#set terminal gif
#set key off
#
set macros
style1 = "with lines lw 1 lt"
style2 = "with lines lw 2 lt"

range1 = "using 1:1" # x:y:w
range2 = "using 1:2" # x:y:w
range3 = "using 1:3" # x:y:w
range4 = "using 1:4" # x:y:w
range5 = "using 1:5" # x:y:w
range6 = "using 1:6" # x:y:w
range7 = "using 1:7" # x:y:w
range8 = "using 1:8" # x:y:w
range9 = "using 1:9" # x:y:w

RED    = "1"
GREEN  = "2"
BLUE   = "3"
PURPLE = "4"
AQUA   = "5"

set xtics 0.1

fin(file) = sprintf("dir-dat/%s.dat", file)
fout(file) = sprintf("dir-gif/%s.gif", file)

file1  = "sinc"
file2  = "saiteki_00060"
#file2  = "saiteki_00101"
title1 = "sin()"
title2 = "IDP"
############################################
set xlabel '[sec]'
set ylabel '[m]'
prefix="z"
set title prefix
plot fin(file1) @range2 @style1 @RED  title title1,\
     fin(file2) @range2 @style1 @BLUE title title2,\
     0 lt AQUA
     
set output fout(prefix)
set terminal gif
replot
set terminal wxt
###########################################
set xlabel '[sec]'
set ylabel '[m/s]'
prefix = "zv"
set title prefix
plot fin(file1) @range3 @style1 @RED  title "sin()",\
     fin(file2) @range3 @style1 @BLUE title "IDP",\
     0 lt AQUA
set output fout(prefix)
set terminal gif
replot
set terminal wxt
###########################################
set xlabel '[sec]'
set ylabel '[m/s^2]'
prefix = "za"
set title prefix
plot fin(file1) @range4 @style1 @RED  title title1,\
     fin(file2) @range4 @style1 @BLUE title title2,\
     0 lt AQUA
set output fout(prefix)
set terminal gif
replot
set terminal wxt
############################################
set xlabel '[sec]'
set ylabel '[volt]'
prefix = "ev"
set title prefix
plot fin(file1) @range5 @style1 @RED  title title1,\
     fin(file2) @range5 @style1 @BLUE title title2,\
     0 lt AQUA     

set output fout(prefix)
set terminal gif
replot
set terminal wxt
#############################################
set xlabel '[sec]'
set ylabel '[amp]'
prefix = "ia"
set title prefix
plot fin(file1) @range6 @style1 @RED  title title1,\
     fin(file2) @range6 @style1 @BLUE title title2,\
     0 lt AQUA     

set output fout(prefix)
set terminal gif
replot
set terminal wxt
#############################################
set xlabel '[sec]'
set ylabel '[Nm]'
prefix = "tau"
set title prefix
plot fin(file1) @range7 @style1 @RED  title title1,\
     fin(file2) @range7 @style1 @BLUE title title2,\
     0 lt AQUA

set output fout(prefix)
set terminal gif
replot
set terminal wxt
#############################################
set xlabel '[sec]'
set ylabel '[J]'
prefix = "energy"
set title prefix
plot fin(file1) @range8 @style1 @RED title title1,\
     fin(file2) @range8 @style1 @BLUE title title2,\
     0 lt AQUA
set output fout(prefix)
set terminal gif
replot
set terminal wxt

#############################################
set xlabel '[sec]'
set ylabel '[J/s]'
prefix = "ev_ia"
set title prefix
plot fin(file1) @range9 @style1 @RED  title title1,\
     fin(file2) @range9 @style1 @BLUE title title2,\
     0 lt AQUA
set output fout(prefix)
set terminal gif
replot
set terminal wxt
################################################
# dir-dat/enemin2.dat
################################################
set xtics 5
set key right bottom
set xlabel '[Time]'
set ylabel '[J]'
prefix = "repeat"
set title prefix
plot fin("enemin2") @range2 @style1 @RED title "energy",\
     0 lt AQUA
set output fout(prefix)
set terminal gif
replot
set terminal wxt
################################################
set key right top
set xlabel '[Time]'
set ylabel '[m]'
prefix = "kiza"
set title prefix
plot fin("enemin2") @range3 @style1 @RED title "kiza",\
     0 lt AQUA
set output fout(prefix)
set terminal gif
replot
set terminal wxt
################################################
