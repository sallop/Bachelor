set size 0.6, 0.6
set grid
set terminal gif

set output "dir-gif/saiteki_energy.gif"
ncol = 2
set xlabel "(sec)"
set ylabel "(J)"
set title "saiteki_energy"
plot "dir-dat/saiteki_00001.dat" u 1:ncol w l title "nnn=1" ,\
     "dir-dat/saiteki_00010.dat" u 1:ncol w l title "nnn=10",\
     "dir-dat/saiteki_00020.dat" u 1:ncol w l title "nnn=20",\
     "dir-dat/saiteki_00030.dat" u 1:ncol w l title "nnn=30",\
     "dir-dat/saiteki_00040.dat" u 1:ncol w l title "nnn=40",\
     "dir-dat/saiteki_00050.dat" u 1:ncol w l title "nnn=50",\
     "dir-dat/saiteki_00060.dat" u 1:ncol w l title "nnn=60"

set output "dir-gif/saiteki_zpos.gif"
ncol = 3     
set xlabel "(sec)"
set ylabel "(m)"
set title "position z"
plot "dir-dat/saiteki_00001.dat" u 1:ncol w l title "nnn=1" ,\
     "dir-dat/saiteki_00010.dat" u 1:ncol w l title "nnn=10",\
     "dir-dat/saiteki_00020.dat" u 1:ncol w l title "nnn=20",\
     "dir-dat/saiteki_00030.dat" u 1:ncol w l title "nnn=30",\
     "dir-dat/saiteki_00040.dat" u 1:ncol w l title "nnn=40",\
     "dir-dat/saiteki_00050.dat" u 1:ncol w l title "nnn=50",\
     "dir-dat/saiteki_00060.dat" u 1:ncol w l title "nnn=60"
     
set output "dir-gif/saiteki_zv.gif"
ncol = 4
set xlabel "(sec)"
set ylabel "(m/sec)"
set title "zv"
plot "dir-dat/saiteki_00001.dat" u 1:ncol w l title "nnn=1" ,\
     "dir-dat/saiteki_00010.dat" u 1:ncol w l title "nnn=10",\
     "dir-dat/saiteki_00020.dat" u 1:ncol w l title "nnn=20",\
     "dir-dat/saiteki_00040.dat" u 1:ncol w l title "nnn=40",\
     "dir-dat/saiteki_00050.dat" u 1:ncol w l title "nnn=50",\
     "dir-dat/saiteki_00060.dat" u 1:ncol w l title "nnn=60"

set output "dir-gif/saiteki_za.gif"
ncol = 5
set xlabel "(sec)"
set ylabel "(m/sec^2)"
set title "za"
plot "dir-dat/saiteki_00001.dat" u 1:ncol w l title "nnn=1" ,\
     "dir-dat/saiteki_00010.dat" u 1:ncol w l title "nnn=10",\
     "dir-dat/saiteki_00020.dat" u 1:ncol w l title "nnn=20",\
     "dir-dat/saiteki_00030.dat" u 1:ncol w l title "nnn=30",\
     "dir-dat/saiteki_00040.dat" u 1:ncol w l title "nnn=40",\
     "dir-dat/saiteki_00050.dat" u 1:ncol w l title "nnn=50",\
     "dir-dat/saiteki_00060.dat" u 1:ncol w l title "nnn=60"

set output "dir-gif/enemin2.gif"
ncol = 2
set xlabel "(nnn)"
set ylabel "(J)"
set title "za"
plot "dir-dat/enemin2.dat" u 1:ncol w l title "enemin2"

set output "dir-gif/kiza1.gif"
ncol = 3
set xlabel "(nnn)"
set ylabel "(m)"
set title "za"
plot "dir-dat/enemin2.dat" u 1:ncol w l title "kiza1"
