set terminal gif
set size 0.6, 0.6

set xlabel "[x]"
set ylabel "[y]"
set title "endeffector_position"
set output "dir-graph/v1_endeffector.gif"
plot "dir-etc/nagato_sinc_endeffector.dat" u 2:3 title "postion" pt 1

set xlabel "[sec]"
set ylabel "[m/sec]"
set title "endeffector_xyv"
set output "dir-graph/v1_endeffector_xyv.gif"
plot "dir-etc/nagato_sinc_endeffector.dat" u 1:5 title "xv" w l,\
     "dir-etc/nagato_sinc_endeffector.dat" u 1:6 title "yv" w l ls 3

set xlabel "[sec]"
set ylabel "[m/sec^2]"
set title "endeffector_xya"
set output "dir-graph/v1_endeffector_xya.gif"
plot "dir-etc/nagato_sinc_endeffector.dat" u 1:8 title "xa" w l,\
     "dir-etc/nagato_sinc_endeffector.dat" u 1:9 title "ya" w l ls 3
