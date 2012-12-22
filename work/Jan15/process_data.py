#!/usr/bin/python
import sys, os

for i in range(len(sys.argv)):
	print i, sys.argv[i]

for i in range(1, len(sys.argv)):
	print i, sys.argv[i]
	if not os.access(sys.argv[1], os.R_OK):
		print "Can't open %s" % (sys.argv[1])
		os.exit()
	
	f = open(sys.argv[i], "r")
	mat = [ map(float, li.strip().split() ) for li in f.readlines()]
	fmt = " ".join( ["%8.4lf" for i in range(len(mat[0]))] )
	# t, z, zv, za, ev, ia, tau, energy, ev_ia
	
	for row in mat:
		print fmt%tuple(row)
	
	f.close()
	# end for i range(1, len(sys, argv)

