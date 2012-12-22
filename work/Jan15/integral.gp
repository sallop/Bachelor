set terminal gif
set size 0.6, 0.6
set grid

set macros

##IDir = "dir-zcurve"
IDir = "dir-etc"

ODir = "dir-zcurve/pderiv"	# parametric derivation
set title "pderiv"
#ODir = "dir-zcurve/deriv"    # total differential
#set title "deriv"
#ODir = "dir-zcurve/chain"    # total differential
#set title "chain"

base1 = "zcurve_r3_sinc_st"
base2 = "compare"
# base1
op1 = "u 1:2 w l title 'zz'"
op2 = "u 1:3 w l title 'zv'"
op3 = "u 1:4 w l title 'za'"
# base2
op4 = "u 1:2 w l title 'zz'    ls 1"
op5 = "u 1:3 w l title 'izvdt' ls 3"

# making file name
fname(_dir, _base, _ext) = sprintf("%s/%s.%s", _dir, _base, _ext)

# zcurve_zz
set xlabel '[sec]'
set ylabel '[m]'
set terminal gif
set output fname(ODir, "zcurve_zz", "gif")
plot fname(IDir, base1, "dat") @op1
pause 0.3

# zcurve_zv
set xlabel '[sec]'
set ylabel '[m/sec]'
set output fname(ODir, "zcurve_zv", "gif")
plot fname(IDir, base1, "dat") @op2

pause 0.3

# zcurve_za
set xlabel '[sec]'
set ylabel '[m/sec^2]'
set output fname(ODir, "zcurve_za", "gif")
plot fname(IDir, base1, "dat") @op3

pause 0.3

# zcurve_zz_izvdt
set xlabel '[sec]'
set ylabel '[m]'
set output fname(ODir, "zz_izvdt", "gif")
plot fname(IDir, base2, "dat") @op4,\
     fname(IDir, base2, "dat") @op5
