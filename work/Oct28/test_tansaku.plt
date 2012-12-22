set terminal gif
set grid
set size 0.6, 0.6

set output "t_th.gif"
set xlabel "[sec]"
set ylabel "[rad]"
plot "t_th.dat" u 1:2 w l title "t_th"

set output "t_thv.gif"
set xlabel "[sec]"
set ylabel "[rad/sec]"
plot "t_th.dat" u 1:3 w l title "t_thv"

set output "t_tha.gif"
set xlabel "[sec]"
set ylabel "[rad/sec^2]"
plot "t_th.dat" u 1:4 w l title "t_tha"

set output "t_Ene.gif"
set xlabel "[sec]"
set ylabel "[J]"
plot "t_th.dat" u 1:5 w l title "t_Ene"

set output "t_ev.gif"
set xlabel "[sec]"
set ylabel "[volt]"
plot "t_e3.dat" u 1:2 w l

set output "t_ia.gif"
set xlabel "[sec]"
set ylabel "[ampere]"
plot "t_e3.dat" u 1:3 w l

set output "t_tau.gif"
set xlabel "[sec]"
set ylabel "[Nm]"
plot "t_e3.dat" u 1:4 w l

set output "t_energy.gif"
set xlabel "[sec]"
set ylabel "[J]"
plot "t_e3.dat" u 1:5 w l
