#!/usr/bin/ruby

# Keep it simple stupid.

require 'rubygems'
require 'gnuplot'

$odir   = "dir-zcurve/pderiv"
$idir   = "dir-etc"
$prefix = "zcurve"

# input file data
a_flag = [ 1, 2, 3, 4, 5, 6, 7, 8, 9 ] # Animation number
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
file_st.each do |ftemp, nums|
  if nums.empty?
    #pdata += [ sprintf(ftemp, $idir, $prefix) ]
    pdata += nums.collect{|num| sprintf(ftemp, $idir, $prefix) }
  elsif
    #nums.each{|num| pdata += [ sprintf(ftemp, $idir, $prefix, num) ] }
    pdata += nums.collect{|num| sprintf(ftemp, $idir, $prefix, num)}
  end
end

# constituent  t, ev, term1, term2, term3
# standrad     t,  z, zv   , za   , ev   , ia, tau , energy, ev_ia
$options = {
  # constituent
  "t_ev"     => {"opt"=>"u 1:2 w l", "ylabel" => "[volt]"  },
  "t_b1-tha" => {"opt"=>"u 1:3 w l", "ylabel" => "b1[volt]"},
  "t_b2-thv" => {"opt"=>"u 1:4 w l", "ylabel" => "b2[volt]"},
  "t_b3-tau" => {"opt"=>"u 1:5 w l", "ylabel" => "b3[volt]"},
  # standard
  "t_z"      => {"opt"=>"u 1:2 w l", "ylabel" => "[m]"      },
  "t_zv"     => {"opt"=>"u 1:3 w l", "ylabel" => "[m/sec]"  },
  "t_za"     => {"opt"=>"u 1:4 w l", "ylabel" => "[m/sec^2]"},
  "t_th"     => {"opt"=>"u 1:2 w l", "ylabel" => "[rad]"      },
  "t_thv"    => {"opt"=>"u 1:3 w l", "ylabel" => "[rad/sec]"  },
  "t_tha"    => {"opt"=>"u 1:4 w l", "ylabel" => "[rad/sec^2]"},
  "t_ev"     => {"opt"=>"u 1:5 w l", "ylabel" => "ev[volt]"   },
  "t_ia"     => {"opt"=>"u 1:6 w l", "ylabel" => "ia[amp]"    },
  "t_tau"    => {"opt"=>"u 1:7 w l", "ylabel" => "[Nm]"       },
  "t_energy" => {"opt"=>"u 1:8 w l", "ylabel" => "[J]"        },
  "t_evia"   => {"opt"=>"u 1:9 w l", "ylabel" => "[watt]"     },  
}

# draw graph
def draw_graph(gp, fname, which)
  opt    = $options[which]['opt']
  ylabel = $options[which]['ylabel']
  
  base = fname.gsub(/#{$idir}/,"").gsub(/\.dat$/,"")
  if base =~ /(\d+)/
    num = $1
  else
    num = "sinc"
  end
  ofname = "#{$odir}/#{base}_#{which}.gif"
  # ["set xtics 0.02"         ,
  #  "set ytics 0.01"         ,
  #  "set xrange [-0.06:0.06]",
  #  "set yrange [ 0.00:0.14]",
  #  "set view 67.0, 340"     ,
  [
   "set terminal gif"      ,
   "set output '#{ofname}'",
   "set xlabel '[sec]'"    ,
   "set ylabel '#{ylabel}'",
   "plot '#{fname}' #{opt} title '#{num}'"
   #"splot '#{fname}' #{opt} title '#{num}'"
  ].each{|cmd| gp.puts(cmd) }
end

def draw_graph_cmp(gp, fnames, nums, which)
  opt    = $options[which]['opt']
  ylabel = $options[which]['ylabel']
  
  base = fnames[0].gsub(/#{$idir}/,"").gsub(/\.dat$/,"")
  if base =~ /(\d+)/
    num = $1
  else
    num = "sinc"
  end
  
  ofname  = "#{$odir}/#{base}_#{which}.gif"
  sendcmd = nums.collect{|num|
    "'#{fnames[num]}' #{opt} title '#{num}'"
  }.join(',')
  
  ["set output '#{ofname}'",
   "set xlabel '[sec]'"    ,
   "set ylabel '#{ylabel}'",
   "plot #{sendcmd}"       ,
  ].each{|cmd| gp.puts(cmd) }
end

############################################################
Gnuplot.open{|gp|
  [
   "set grid"         ,
   "set size 0.6, 0.6",
  ].each{|cmd| gp.puts(cmd) }

  # oneline
  #pdata.each {|file|
  #  $options.keys.each{|data| draw_graph(gp, file, data) }
  #}
  # compare
  cmpnum = [0, 3, 6, 9]
  #$options.keys.each{|data| draw_graph_cmp(gp, pdata, cmpnum, data) }
  
  cmpfiles = [
              sprintf("%s/%s_r3_sinc_st.dat", $idir, $prefix),
              cmpnum.map{|num|
                sprintf("%s/%s_r3_saiteki_st_%05d.dat", $idir, $prefix, num)
              }
             ].flatten()

  sincfile = sprintf("%s/%s_r3_sinc_st.dat", $idir, $prefix)
  #  idpfiles = cmpnum.map{|num|
  idpfiles = a_flag.map{|num|
    sprintf("%s/%s_r3_saiteki_st_%05d.dat", $idir, $prefix, num)
  }
  result   = idpfiles.pop
  
  l_opt = ["t_z" , "t_zv" , "t_za" ,
           "t_th", "t_thv", "t_tha",
           "t_ev", "t_ia" , "t_tau",
           "t_energy", "t_evia"]

  
  sendcmd  = idpfiles.map{|f| "'#{f}' u 1:3 w l ls 2 title ''"}.join(',')
  sendcmd += ",'#{sincfile}' u 1:3 w l ls 1 title 'sin'"
  sendcmd += ", '#{result}' u 1:3 w l ls 3 title 'last'"
#  sendcmd += ",'#{sincfile}' u 1:3 w l ls 1 title 'sin'"
  
  print sendcmd
  [
#   "set terminal gif"                          ,
#   "set output '#{$odir}/#{$prefix}_cmpzv.gif'",
   "set xlabel '[sec]'"                        ,
   "set ylabel '[m/sec]'"                      ,
   "plot #{sendcmd}"                           ,
   "pause 0.9"                                 ,
  ].each{|cmd| gp.puts(cmd) }

  
  
  gp.puts "pause 1.8"
}
