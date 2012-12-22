#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;

open( my $fh, "<", "dir-zcurve/zcurve_r3_sinc_st.dat")
  or die "Cannot open file: $!";

my $dt = 0.001;
my (@tt, @zz, @zv, @za);
#my @zv = ();
#my @za = ();
# @data = <FH>;
while( my $line = readline $fh ){
    chomp $line;
    my @list = split(/ /, $line);
    push(@tt, $list[0]);
    push(@zz, $list[1]);
    push(@zv, $list[2]);
    push(@zz, $list[3]);
}

close($fh);


open(my $fo, ">", "dir-zcurve/integral_zv.dat");
my (@zvdt, @zadt);

foreach my $zv( @zv ){
    push(@zvdt, $zv*$dt);
}

my @foo = (\@zz,\@zvdt);


for(my $i=0; $i <= $#zz; $i++){
    print $zz[$i]," ", $zvdt[$i], "\n";
}

print "@tt=",$#tt,"\n";
print "@zz=",$#zz,"\n";
print "@zv=",$#zv,"\n";
print "@za=",$#za,"\n";


#print @data;
