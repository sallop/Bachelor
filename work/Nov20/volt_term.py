#!/usr/bin/python
import sys, os
import subprocess

for i in range(len(sys.argv)):
	print i, sys.argv[i]

for current_file in sys.argv[1:]:
	if not os.access(current_file, os.R_OK):
		print "Can't open %s" % (current_file)
		os.exit()
	
	fr = open(current_file, "r")
	mat = [ map(float, li.strip().split() ) for li in fr.readlines() ]
	fr.close()

	fmt = " ".join( ["%8.4lf" for spam in range(len(mat[0]))] )
	# t, volt, b1, b2, b3
	basename = current_file[ current_file.index('/')+1: ]
	fname = "dir-dat/hist_%s" % (basename)
	fw = open(fname, "w")
	for row in mat:
		(t, volt, b1, b2, b3) = row
		print >> fw, fmt%(t, b1, b1+b2, b1+b2+b3, volt)
	fw.close()
	# end for i range(1, len(sys, argv)

