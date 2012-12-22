set terminal gif
set size 0.6, 0.6

set xlabel "[sec]"
set ylabel "[rad]"
set title "r1r2_th"
set output "dir-graph/r1r2_th.gif"
plot "dir-etc/nagato_r1_sinc_st.dat" u 1:2 w l title "r1.th",\
     "dir-etc/nagato_r2_sinc_st.dat" u 1:2 w l title "r2.th" ls 3

set xlabel "[sec]"
set ylabel "[rad/sec]"
set title "r1r2_thv"
set output "dir-graph/r1r2_thv.gif"
plot "dir-etc/nagato_r1_sinc_st.dat" u 1:3 w l title "r1.thv",\
     "dir-etc/nagato_r2_sinc_st.dat" u 1:3 w l title "r2.thv" ls 3

set xlabel "[sec]"
set ylabel "[rad/sec^2]"
set title "r1r2_tha"
set output "dir-graph/r1r2_tha.gif"
plot "dir-etc/nagato_r1_sinc_st.dat" u 1:4 w l title "r1.tha",\
     "dir-etc/nagato_r2_sinc_st.dat" u 1:4 w l title "r2.tha" ls 3

set xlabel "[sec]"
set ylabel "[volt]"
set title "r1r2_ev"
set output "dir-graph/r1r2_ev.gif"
plot "dir-etc/nagato_r1_sinc_st.dat" u 1:5 w l title "r1.ev",\
     "dir-etc/nagato_r2_sinc_st.dat" u 1:5 w l title "r2.ev" ls 3

set xlabel "[sec]"
set ylabel "[amp]"
set title "r1r2_ia"
set output "dir-graph/r1r2_ia.gif"
plot "dir-etc/nagato_r1_sinc_st.dat" u 1:6 w l title "r1.ia",\
     "dir-etc/nagato_r2_sinc_st.dat" u 1:6 w l title "r2.ia" ls 3

set xlabel "[sec]"
set ylabel "[Nm]"
set title "r1r2_tau"
set output "dir-graph/r1r2_tau.gif"
plot "dir-etc/nagato_r1_sinc_st.dat" u 1:7 w l title "r1.tau",\
     "dir-etc/nagato_r2_sinc_st.dat" u 1:7 w l title "r2.tau" ls 3

set key left top
set xlabel "[sec]"
set ylabel "[J]"
set title "r1r2_energy"
set output "dir-graph/r1r2_energy.gif"
plot "dir-etc/nagato_r1_sinc_st.dat" u 1:8 w l title "r1.energy",\
     "dir-etc/nagato_r2_sinc_st.dat" u 1:8 w l title "r2.energy" ls 3

set key left top
set xlabel "[sec]"
set ylabel "[watt]"
set title "r1r2_evia"
set output "dir-graph/r1r2_evia.gif"
plot "dir-etc/nagato_r1_sinc_st.dat" u 1:9 w l title "r1.ev_ia",\
     "dir-etc/nagato_r2_sinc_st.dat" u 1:9 w l title "r2.ev_ia" ls 3
