#!/usr/bin/ruby

# Keep it simple stupid.

require 'rubygems'
require 'gnuplot'

$odir   = "dir-zcurve/pderiv"
$idir   = "dir-etc"
$prefix = "zcurve"

def fstat(fname, xlabel, ylabel, title, option)
  { "fname"   => fname,
    "xlabel"  => xlabel,
    "ylabel"  => ylabel,
    "title"   => title ,
    "option"  => option }
end

# input file data
a_flag = [ 1, 2, 3, 4, 5, 6, 7, 8, 9 ]
file_st  = [#["%s/%s_sinc_st.dat"          , []     ],
            ["%s/%s_r1_sinc_st.dat"        , []     ],
            ["%s/%s_r2_sinc_st.dat"        , []     ],
            ["%s/%s_r3_sinc_st.dat"        , []     ],
            #["%s/%s_saiteki_st_%05d.dat"  , a_flag ],
            ["%s/%s_r1_saiteki_st_%05d.dat", a_flag ],
            ["%s/%s_r2_saiteki_st_%05d.dat", a_flag ],
            ["%s/%s_r3_saiteki_st_%05d.dat", a_flag ],   ]
file_ev  = [#["%s/%s_sinc_ev.dat"          , []     ],
            ["%s/%s_r1_sinc_ev.dat"        , []     ],
            ["%s/%s_r2_sinc_ev.dat"        , []     ],
            ["%s/%s_r3_sinc_ev.dat"        , []     ],
            #["%s/%s_saiteki_ev_%05d.dat"  , a_flag ],
            ["%s/%s_r1_saiteki_ev_%05d.dat", a_flag ],
            ["%s/%s_r2_saiteki_ev_%05d.dat", a_flag ],
            ["%s/%s_r3_saiteki_ev_%05d.dat", a_flag ],   ]
file_ee  = [["%s/%s_sinc_endeffector.dat"        , []    ],
            ["%s/%s_saiteki_endeffector_%05d.dat", a_flag],]
file_etc = [["%s/%s_enemin2.dat"                 , []    ],
            ["%s/%s_energy_sinc.dat"             , []    ],
            ["%s/%s_energy_saiteki_%05d.dat"     , a_flag],]

pdata = []
file_ee.each do |ftemp, nums|
  if nums.empty?
    pdata += [ sprintf(ftemp, $idir, $prefix) ]
  elsif
    nums.collect{|num|
      pdata += [ sprintf(ftemp, $idir, $prefix, num) ]
    }
  end
end

# pdata = file_ee.collect{|ftemp, nums|
#   if nums.empty?
#     sprintf(ftemp, $idir, $prefix)
#   elsif
#     nums.collect{|num| sprintf(ftemp, $idir, $prefix, num)}.join(',')
#   end
# }.join(',')
# p pdata
#exit

$options = {
  'trjc' => {'opt'=>'u 2:4 w l', 'ylabel'=>'[m]'      , 'title'=>''},
  #'trjc' => {'opt'=>'u 2:3:4 w l', 'ylabel'=>'[m]'      , 'title'=>''},
  'xv'   => {'opt'=>'u 1:5 w l', 'ylabel'=>'xv[m/sec]', 'title'=>''},
  'yv'   => {'opt'=>'u 1:6 w l', 'ylabel'=>'yv[m/sec]', 'title'=>''},
  'zv'   => {'opt'=>'u 1:7 w l', 'ylabel'=>'zv[m/sec]', 'title'=>''},
}
  
# draw graph
def draw_trajectory(gp, fname)
  opt = $options['trjc']['opt']  #"u 2:3:4 w l"

  base = fname.gsub(/#{$idir}/,"").gsub(/\.dat$/,"")
  if base =~ /(\d+)/
    num = $1
  else
    num = "sinc"
  end
  ofname = "#{$odir}/#{base}_trjc.gif"
  [
   "set terminal gif"         ,
   "set output '#{ofname}'"   ,
   # "set xtics 0.02"         ,
   # "set ytics 0.01"         ,
   # "set xrange [-0.06:0.06]",
   # "set yrange [ 0.00:0.14]",
   # "set zlabel '[m]'"
   # "set view 67.0, 340"     ,
   "set xlabel 'x_axis[m]'"   ,
   "set ylabel 'z_axis[m]'"   , 
   "plot '#{fname}' u 2:4 w l title '#{num}'"  ,
   #"splot '#{fname}' u 2:4 w l title '#{num}'",
  ].each{|cmd| gp.puts(cmd) }
end

def draw_velocity(gp, fname, xyz)
  opt = $options[xyz]['opt']
  ylabel = $options[xyz]['ylabel']
  base = fname.gsub(/#{$idir}/,"").gsub(/\.dat$/,"")
  if base =~ /(\d+)/
    num = $1
  else
    num = "sinc"
  end
  
  ofname = "#{$odir}/#{base}_#{xyz}.gif"
  [
   "set output '#{ofname}'",
   "set xlabel '[sec]'"    ,
   "set ylabel '#{ylabel}'",
   #"set xrange [0.0:0.7]" ,
   "plot '#{fname}' #{opt} title '#{num}'",
  ].each{|cmd| gp.puts(cmd) }
end

def draw_velo_cmp(gp, fnames, xyz, nums)
  opt = $options[xyz]['opt']
  ylabel = $options[xyz]['ylabel']
  
  base = "cmpvelo_#{nums}"
  if base =~ /(\d+)/
    num = $1
  else
    num = "sinc"
  end
  
  ofname = "#{$odir}/#{base}_#{xyz}.gif"
  sendcmd = nums.collect{|num|
    "'#{fnames[num]}' #{opt} title '#{xyz}#{num}'"
  }.join(',')
  
  [
   "set output '#{ofname}'"    ,
   "set xlabel '[sec]'"        ,
   "set ylabel '#{xyz}[m/sec]'",
   "plot #{sendcmd}"           ,
  ].each{|cmd| gp.puts(cmd) }
end

def draw_trjc_cmp(gp, fnames, nums)
  opt = $options['trjc']['opt']
  ylabel = $options['trjc']['ylabel']
  
  base = "trjc_cmp#{nums}"
  fnames[0].gsub(/#{$idir}/,"").gsub(/\.dat$/,"")

  ofname = "#{$odir}/#{base}.gif"
  sendcmd = nums.each{|num|
    "'#{fnames[num]}' #{opt} title '#{num}'"}.join(',')
  [
   "set terminal gif"      ,
   "set output '#{ofname}'",
   "set xlabel 'x[m]'",
   "set ylabel 'z[m]'",
   #"set zlabel '[m]'"     ,
   "plot #{sendcmd}",
  ].each{|cmd| gp.puts(cmd) }
end


Gnuplot.open{|gp|
  gp.puts("set grid")
  gp.puts("set size 0.6, 0.6\n")
  # one line
  pdata.each do |fin|
    draw_trajectory(gp, fin)
    ['xv','yv','zv'].each{|xyz| draw_velocity(gp, fin, xyz) }
  end

  # compare
  cmpnum = [0, 3, 6, 9]
  draw_trjc_cmp(gp, pdata, cmpnum)
  ['xv','yv','zv'].each{|xyz| draw_velo_cmp(gp, pdata, xyz, cmpnum) }
  
  gp.puts "pause 1.8"
}
