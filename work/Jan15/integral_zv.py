#!/usr/bin/python

#fi = open('dir-zcurve/zcurve_r3_sinc_st.dat')
#fi = open('dir-etc/zcurve_r3_sinc_st.dat')
fi = open('dir-etc/zcurve_r3_saiteki_st_00009.dat')
#fi = open('dir-etc/r3zv.dat','r')
lines = fi.readlines()
fi.close()

for line in lines:
    line = line.rstrip('\n')
#    nums = line.split(' ')
    nums = line.split()
    print nums

#strnums = [line.rstrip('\n').split(' ')  for line in lines]
strnums = [line.rstrip('\n').split()  for line in lines]

#nums = [( float(strnum[0]),
#          float(strnum[1]),
#          float(strnum[2]),
#          float(strnum[3]) ) for strnum in strnums ]

nums = [( float(strnum[0]),
          float(strnum[1]),
          float(strnum[2]),
          float(strnum[3]) ) for strnum in strnums ]

(tt, zz, zv, za, zvdt, zadt, izvdt, izadt) = [], [], [], [], [], [], [], []
for strnum in strnums:
    tt.append( float(strnum[0]) )
    zz.append( float(strnum[1]) )
    zv.append( float(strnum[2]) )
    za.append( float(strnum[3]) )

print len(tt)
print len(zz)
print len(zv)
print len(za)

dt = 0.001

# maybe this program right.

for t, z, v in zip(tt, zz, zv):
#    print "%lf %lf %lf"%(t, z, v)
    zvdt.append( v*dt )
for i in range( len(zvdt) ):
    sum = 0
    for j in range(i):
        sum += zvdt[j]
    izvdt.append(sum)

for t, v, a in zip(tt, zv, za):
    zadt.append( a*dt )
for i in range( len(zadt) ):
    sum = 0
    for j in range(i):
        sum += zadt[j]
    izadt.append(sum)


#fo = open("dir-zcurve/compare.dat",'w')
fo = open("dir-etc/compare.dat",'w')

for t, z, iz, v, iv in zip(tt, zz, izvdt, zv, izadt):
#    print >> fo,"%lf %lf %lf"%(t, x, y)
    print >> fo,"%lf %lf %lf %lf %lf"%(t, z, iz, v, iv)
    
fo.close()
