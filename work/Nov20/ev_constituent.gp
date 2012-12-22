#
set size 0.6, 0.6
set grid
#set terminal gif
#set key off
#
set macros
style1 = "with lines lw 1 lt"
style2 = "with lines lw 2 lt"

range1 = "using 1:1" # x:y:w t
range2 = "using 1:2" # x:y:w ev
range3 = "using 1:3" # x:y:w b1*thv
range4 = "using 1:4" # x:y:w b2*tha
range5 = "using 1:5" # x:y:w b3*tau
range6 = "using 1:6" # x:y:w
range7 = "using 1:7" # x:y:w
range8 = "using 1:8" # x:y:w
range9 = "using 1:9" # x:y:w

RED    = "1"
GREEN  = "2"
BLUE   = "3"
PURPLE = "4"
AQUA   = "5"

#set xtics 0.1

fin(file) = sprintf("dir-dat/%s.dat", file)
fout(file) = sprintf("dir-gif/%s.gif", file)

file1  = "sinc_ev"
file2  = "saiteki_ev_00060"
#file2  = "saiteki_00101"
title1 = "sin()"
title2 = "IDP"

######################################################
set title "line graph sin()"
set xlabel "[sec]"
set ylabel "[V]"
plot fin(file1) u 1:3 w l lt @RED   title "b1*thv",\
     fin(file1) u 1:4 w l lt @GREEN title "b2*tha",\
     fin(file1) u 1:5 w l lt @BLUE  title "b3*tau",\
     fin(file1) u 1:2 w l lt @PURPLE title "volt"
set output fout(file1."line")
set terminal gif
replot
set terminal wxt
######################################################
set title "line graph IDP"
set xlabel "[sec]"
set ylabel "[V]"
plot fin(file2) u 1:3 w l lt @RED   title "b1*thv",\
     fin(file2) u 1:4 w l lt @GREEN title "b2*tha",\
     fin(file2) u 1:5 w l lt @BLUE  title "b3*tau",\
     fin(file2) u 1:2 w l lt @PURPLE title "volt"

set output fout(file2."line")
set terminal gif
replot
set terminal wxt
##################################################


set style data histograms
set style histogram rowstacked
set style fill solid 1.00 noborder
########################################################
set title "histograms sin()"
set xlabel ""
set ylabel "[V]"
plot fin(file1) u 3 lt @RED   title "b1*thv",\
     fin(file1) u 4 lt @GREEN title "b2*tha",\
     fin(file1) u 5 lt @BLUE  title "b3*tau",\
     fin(file1) u 2 with lines lt @PURPLE title "volt"

set output fout(file1)
set terminal gif
replot
set terminal wxt
########################################################
set title "histograms IDP"
set xlabel ""
set ylabel "[V]"
plot fin(file2) u 3 lt @RED   title "b1*thv",\
     fin(file2) u 4 lt @GREEN title "b2*tha",\
     fin(file2) u 5 lt @BLUE  title "b3*tau",\
     fin(file2) u 2 with lines lt @PURPLE title "volt"

set output fout(file2)
set terminal gif
replot
set terminal wxt