# using gnuplot file

set size 0.6, 0.6
set grid
set terminal gif

### thv,tha,tau
set output "plot/t_thv1_tha1_tau1.gif"
set xlabel "(sec)"
plot	"motor.dat" u 1:3 w l title "thv1" 1,\
	"motor.dat" u 1:4 w l title "tha1" 3,\
	"motor.dat" u 1:5 w l title "tau1" 4

set output "plot/t_thv1.gif"
set xlabel "(sec)"
set ylabel "(rad/sec)
plot	"motor.dat" u 1:3 w l title "thv1" 1

set output "plot/t_tha1.gif"
set xlabel "(sec)"
set ylabel "(rad/sec**2)
plot	"motor.dat" u 1:4 w l title "tha1" 3
	
set output "plot/t_thv2_tha2_tau2.gif"
set xlabel "(sec)"
plot	"motor.dat" u 1:7 w l title "thv2" 1,\
	"motor.dat" u 1:8 w l title "tha2" 3,\
	"motor.dat" u 1:9 w l title "tau2" 4

set output "plot/t_thv2.gif"
set xlabel "(sec)"
set ylabel "(rad/sec)
plot	"motor.dat" u 1:7 w l title "thv2" 1

set output "plot/t_tha2.gif"
set xlabel "(sec)"
set ylabel "(rad/sec**2)
plot	"motor.dat" u 1:8 w l title "tha2" 3

	
### energy ###     
## Ene_[12]
set output "plot/t_energy1.gif"
set xlabel "(sec)"
set ylabel "(J)"
plot "energy1.dat" u 1:2 w l title "energy1"

set output "plot/t_energy2.gif"
set xlabel "(sec)"
set ylabel "(J)"
plot "energy2.dat" u 1:2 w l title "energy2"

set output "plot/t_energy12.gif"
set xlabel "(sec)"
set ylabel "(J)"
plot	"energy1.dat" u 1:2 w l title "energy1" 1,\
	"energy2.dat" u 1:2 w l title "energy2" 3
## ev[12]
set output "plot/t_ev1.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot	"energy1.dat" u 1:3 w l title "ev1"

set output "plot/t_ev2.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot	"energy2.dat" u 1:3 w l title "ev2"

set output "plot/t_ev12.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot	"energy1.dat" u 1:3 w l title "ev1" 1,\
	"energy2.dat" u 1:3 w l title "ev2" 3

## ia[12]
set output "plot/t_ia1.gif"
set xlabel "(sec)"
set ylabel "(Ampere)"
plot "energy1.dat" u 1:4 w l title "ia1"

set output "plot/t_ia2.gif"
set xlabel "(sec)"
set ylabel "(Ampere)"
plot "energy2.dat" u 1:4 w l title "ia2"

set output "plot/t_ia12.gif"
set xlabel "(sec)"
set ylabel "(Ampere)"
plot	"energy1.dat" u 1:4 w l title "ia1" 1,\
	"energy2.dat" u 1:4 w l title "ia2" 3

## tau[12]
set output "plot/t_tau1.gif"
set xlabel "(sec)"
set ylabel "(N*m)"
plot "energy1.dat" u 1:4 w l title "tau1" 4

set output "plot/t_tau2.gif"
set xlabel "(sec)"
set ylabel "(N*m)"
plot "energy2.dat" u 1:4 w l title "tau2" 4

set output "plot/t_tau12.gif"
set xlabel "(sec)"
set ylabel "(N*m)"
plot	"energy1.dat" u 1:4 w l title "tau1" 1,\
	"energy2.dat" u 1:4 w l title "tau2" 3

### factor.dat ###
## fact1
set output "plot/t_b1*thv1.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot "factor1.dat" u 1:2 w l title "b1*thv1" 1

set output "plot/t_b1*thv2.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot "factor2.dat" u 1:2 w l title "b1*thv2" 1

set output "plot/t_b1*thv1_b1*thv2.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot	"factor1.dat" u 1:2 w l title "b1*thv1" 1,\
	"factor2.dat" u 1:2 w l title "b1*thv2" 3

## fact2
set output "plot/t_b2*tha1.gif" 
set xlabel "(sec)"
set ylabel "(Volt)"
plot "factor1.dat" u 1:3 w l title "b2*tha1" 3

set output "plot/t_b2*tha2.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot "factor2.dat" u 1:3 w l title "b2*tha2" 3

set output "plot/t_b2*tha1_b2*tha2.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot	"factor1.dat" u 1:3 w l title "b2*tha1" 1,\
	"factor2.dat" u 1:3 w l title "b2*tha2" 3

## fact3
set output "plot/t_b3*tau1.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot "factor1.dat" u 1:4 w l title "b3*tau1" 4

set output "plot/t_b3*tau2.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot "factor2.dat" u 1:4 w l title "b3*tau2" 4

set output "plot/t_b3*tau1_b3*tau2.gif"
set xlabel "(sec)"
set ylabel "(Volt)"
plot	"factor1.dat" u 1:4 w l title "b3*tau1" 1,\
	"factor2.dat" u 1:4 w l title "b3*tau2" 3

## fact ia 1
set output "plot/t_ev1.gif"
set xlabel "(sec)"
set ylabel "(volt)"
plot "factor1.dat" u 1:5 w l title "ev1"

set output "plot/t_ev2.gif"
set xlabel "(sec)"
set ylabel "(volt)"
plot "factor2.dat" u 1:5 w l title "ev2"

set output "plot/t_ev1_ev2.gif"
set xlabel "(sec)"
set ylabel "(volt)"
plot	"factor1.dat" u 1:5 w l title "ev1" 1,\
	"factor2.dat" u 1:5 w l title "ev2" 3

## fact ia 2
set output "plot/t_Kv*thv1.gif"
set xlabel "(sec)"
set ylabel "(volt)"
plot "factor1.dat" u 1:6 w l title "Kv*thv1"

set output "plot/t_Kv*thv2.gif"
set xlabel "(sec)"
set ylabel "(volt)"
plot "factor2.dat" u 1:6 w l title "Kv*thv2"

set output "plot/t_Kv*thv1_Kv*thv2.gif"
set xlabel "(sec)"
set ylabel "(volt)"
plot	"factor1.dat" u 1:6 w l title "Kv*thv1" 1,\
	"factor2.dat" u 1:6 w l title "Kv*thv2" 3

## fact ia 1_2
set output "plot/t_factor1_ia.gif"
set xlabel "(sec)"
set ylabel "(volt)"
set output "plot/t_factor1_ia.gif"
plot	"factor1.dat" u 1:5 w l title "ev1" 1,\
	"factor1.dat" u 1:6 w l title "Kv*thv1" 3

set output "plot/t_factor2_ia.gif"
set xlabel "(sec)"
set ylabel "(volt)"
set output "plot/t_factor2_ia.gif"
plot	"factor2.dat" u 1:5 w l title "ev2" 1,\
	"factor2.dat" u 1:6 w l title "Kv*thv2" 3

#### compare
set output "plot/t_ev_factor1.gif"
set xlabel "(sec)"
set ylabel "(volt)"
plot	"factor1.dat" u 1:2 w l title "b1*thv1" 1,\
	"factor1.dat" u 1:3 w l title "b2*tha1" 3,\
	"factor1.dat" u 1:4 w l title "b3*tau1" 4

set output "plot/t_ev_factor2.gif"
set xlabel "(sec)"
set ylabel "(volt)"
plot	"factor2.dat" u 1:2 w l title "b1*thv2" 1,\
	"factor2.dat" u 1:3 w l title "b2*tha2" 3,\
	"factor2.dat" u 1:4 w l title "b3*tau2" 4
