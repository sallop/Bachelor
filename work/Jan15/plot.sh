#!/bin/bash

WHICH="pderiv"
#WHICH="deriv"
#WHICH="chain"

[ -d dir-dat ] || mkdir dir-dat
[ -d dir-gif ] || mkdir dir-gif
[ -d dir-etc ] || mkdir dir-etc
[ -d dir-zcurve/${WHICH} ] || mkdir -p dir-zcurve/${WHICH}
# L261 pderiv(parametric) or deriv(total differential)
echo "L261 now=${WHICH} |pderiv or deriv or chain|"
read
vim +261 IDPZCurve.cpp
make pulley && ./pulley || echo "./pulley failed $?" 

#if [ -d dir-zcurve ]; then
#	view +${LINENO} plot.sh
#	mv dir-etc dir-zcurve || echo "mv dir-dat dir-zcurve"
#else
#	view +${LINENO} plot.sh
#fi
#if [ -d dir-zcurve/${WHICH} ]; then
#	view +15 plot.sh
#	echo "dir-zcurve/${WHICH} already exist."
#else
#	view +18 plot.sh
#	mkdir dir-zcurve/${WHICH}
#	echo "mkdir dir-zcurve/${WHICH}"
#fi

echo "${WHICH} integral_zv.py"
read
vim integral_zv.py
python integral_zv.py

#echo "${WHICH} dir-zcurve/compare.dat"
echo "${WHICH} dir-etc/compare.dat"
read
#vim dir-zcurve/compare.dat
vim dir-etc/compare.dat

echo "${WHICH} integral.gp"
read
vim integral.gp
gnuplot integral.gp

echo "${WHICH} r3zv3gp"
read
vim r3zv.gp
gnuplot r3zv.gp

if [ $WHICH = "pderiv" ]; then
	echo $WHICH
	eog dir-zcurve/pderiv &
elif [ $WHICH = "deriv" ]; then
	echo $WHICH
	eog dir-zcurve/deriv &
elif [ $WHICH = "chain" ]; then
    echo $WHICH
    eog dir-zcurve/chain &
fi
