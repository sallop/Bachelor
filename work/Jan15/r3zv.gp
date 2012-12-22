# config
set terminal gif
set size 0.6, 0.6
set grid
set macros
#input director##################################################
##IDir = "dir-zcurve"
IDir = "dir-etc"
#Onput director##################################################
# pderiv
ODir = "dir-zcurve/pderiv"	# parametric derivation
set title "pderiv"
# deriv
#ODir = "dir-zcurve/deriv"    # total differential
#set title "deriv"
# chain
#ODir = "dir-zcurve/chain"    # total differential
#set title "chain"
##################################################
base1 = "r3zv"

# making file name
fname(_dir, _base, _ext) = sprintf("%s/%s.%s", _dir, _base, _ext)
##################################################
set output fname(ODir, "r3zv_zz", "gif")
op2 = "u 1:2 w l title 'zz'"
set xlabel '[sec]'
set ylabel '[m]'
#set terminal gif
plot fname(IDir, base1, "dat") @op2
pause 0.8

##################################################
set output fname(ODir, "r3zv_zv", "gif")
op3 = "u 1:3 w l title 'zv'"
set xlabel '[sec]'
set ylabel '[m/sec]'
#set terminal gif
plot fname(IDir, base1, "dat") @op3
pause 0.8

##################################################
set output fname(ODir, "r3zv_za", "gif")
op4 = "u 1:4 w l title 'za'"
set xlabel '[sec]'
set ylabel '[m/sec^2]'
#set terminal gif
plot fname(IDir, base1, "dat") @op4
pause 0.8

##################################################
set output fname(ODir, "r3zv_v1x", "gif")
op5 = "u 1:5 w l title 'v1x '"
set xlabel '[sec]'
set ylabel '[m]'
#set terminal gif
plot fname(IDir, base1, "dat") @op5
pause 0.8

##################################################
set output fname(ODir, "r3zv_v1xv", "gif")
op6 = "u 1:6 w l title 'v1xv'"
set xlabel '[sec]'
set ylabel '[m/sec]'
#set terminal gif
plot fname(IDir, base1, "dat") @op6
pause 0.8

##################################################
set output fname(ODir, "r3zv_v1xa", "gif")
op7 = "u 1:7 w l title 'v1xa'"
set xlabel '[sec]'
set ylabel '[m/sec^2]'
#set terminal gif
plot fname(IDir, base1, "dat") @op7
pause 0.8

##################################################
set output fname(ODir, "r3zv_r2x2", "gif")
op8 = "u 1:8 w l title 'r2_x2'"
set xlabel '[sec]'
set ylabel 'r2_x2[m]'
#set terminal gif
plot fname(IDir, base1, "dat") @op8
pause 0.8

