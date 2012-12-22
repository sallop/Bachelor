#!/usr/bin/ruby

# Keep it simple stupid.

require 'rubygems'
require 'gnuplot'

$odir   = "dir-zcurve/pderiv"
$idir   = "dir-etc"
$prefix = "zcurve"

# input file data
a_flag = [ 1, 2, 3, 4, 5, 6, 7, 8, 9 ] # IDP repeat number. Using make animation.
file_st  = [["%s/%s_r1_sinc_st.dat"        , []     ],
            ["%s/%s_r2_sinc_st.dat"        , []     ],
            ["%s/%s_r3_sinc_st.dat"        , []     ],
            ["%s/%s_r1_saiteki_st_%05d.dat", a_flag ],
            ["%s/%s_r2_saiteki_st_%05d.dat", a_flag ],
            ["%s/%s_r3_saiteki_st_%05d.dat", a_flag ],   ]

file_ev  = [["%s/%s_r1_sinc_ev.dat"        , []     ],
            ["%s/%s_r2_sinc_ev.dat"        , []     ],
            ["%s/%s_r3_sinc_ev.dat"        , []     ],
            ["%s/%s_r1_saiteki_ev_%05d.dat", a_flag ],
            ["%s/%s_r2_saiteki_ev_%05d.dat", a_flag ],
            ["%s/%s_r3_saiteki_ev_%05d.dat", a_flag ],   ]

file_ee  = [["%s/%s_sinc_endeffector.dat"        , []    ],
            ["%s/%s_saiteki_endeffector_%05d.dat", a_flag],]

file_etc = [["%s/%s_enemin2.dat"                 , []    ],
            ["%s/%s_energy_sinc.dat"             , []    ],
            ["%s/%s_energy_saiteki_%05d.dat"     , a_flag],]

def pArg(opt, xlabel, ylabel, *list)
  if list.size % 2 != 0
    puts "list is not even"
    return {"opt"=>opt, "xlabel"=>xlabel, "ylabel"=>ylabel}
  end
  { "opt"=>opt, "xlabel"=>xlabel, "ylabel"=>ylabel }.merge( Hash[*list] )
end

$options_ev = {
  "t_ev"     => pArg("u 1:2 w l", "[sec]", "[volt]"  ),
  "t_b1-tha" => pArg("u 1:3 w l", "[sec]", "b1[volt]"),
  "t_b2-thv" => pArg("u 1:4 w l", "[sec]", "b2[volt]"),
  "t_b3-tau" => pArg("u 1:5 w l", "[sec]", "b3[volt]"),
}

$options_st = {
  # standard
  "t_z"      => pArg("u 1:2 w l", "[sec]", "[m]"        ),
  "t_zv"     => pArg("u 1:3 w l", "[sec]", "[m/sec]"    ),
  "t_za"     => pArg("u 1:4 w l", "[sec]", "[m/sec^2]"  ),
  "t_th"     => pArg("u 1:2 w l", "[sec]", "[rad]"      ),
  "t_thv"    => pArg("u 1:3 w l", "[sec]", "[rad/sec]"  ),
  "t_tha"    => pArg("u 1:4 w l", "[sec]", "[rad/sec^2]"),
  "t_ev"     => pArg("u 1:5 w l", "[sec]", "ev[volt]"   ),
  "t_ia"     => pArg("u 1:6 w l", "[sec]", "ia[amp]"    ),
  "t_tau"    => pArg("u 1:7 w l", "[sec]", "[Nm]"       ),
  "t_energy" => pArg("u 1:8 w l", "[sec]", "[J]"        ),
  "t_evia"   => pArg("u 1:9 w l", "[sec]", "[watt]"     ),  
}

$options_ee = {
  # endeffector
  "t_x"  => pArg("u 1:2  w l", "[sec]", "[m]"      ),
 #"t_y"  => pArg("u 1:3  w l", "[sec]", "[m]"      ),
  "t_z"  => pArg("u 1:4  w l", "[sec]", "[m]"      ),
  "t_xv" => pArg("u 1:5  w l", "[sec]", "[m/sec]"  ),
  "t_yv" => pArg("u 1:6  w l", "[sec]", "[m/sec]"  ),
  "t_zv" => pArg("u 1:7  w l", "[sec]", "[m/sec]"  ),
  "t_xa" => pArg("u 1:8  w l", "[sec]", "[m/sec^2]"),
  "t_ya" => pArg("u 1:9  w l", "[sec]", "[m/sec^2]"),
  "t_za" => pArg("u 1:10 w l", "[sec]", "[m/sec^2]"),
  "x_z"  => pArg("u 2:4  w l", "[m]"  , "[m]"      ),
}

$options_etc = {
  # enemin2
  "n_e1e2e3" => pArg("u 1:2 w l title 'n_enemin'", "[times]", "[J]"),
  "n_e1"     => pArg("u 1:3 w l", "[times]", "[J]"),
  "n_e2"     => pArg("u 1:4 w l", "[times]", "[J]"),
  "n_e3"     => pArg("u 1:5 w l", "[times]", "[J]"),
  "n_kiza"   => pArg("u 1:6 w l", "[times]", "[m]"),
  # energy
  "t_e1e2e3" => pArg("u 1:2 w l", "[sec]"  , "[J]"),
  "t_e1"     => pArg("u 1:3 w l", "[sec]"  , "[J]"),
  "t_e2"     => pArg("u 1:4 w l", "[sec]"  , "[J]"),
  "t_e3"     => pArg("u 1:5 w l", "[sec]"  , "[J]"),
}

$options = {
  "ev"  => $options_ev ,
  "st"  => $options_st ,
  "ee"  => $options_ee ,
  "etc" => $options_etc,
}

def getInputFiles(ftemps)
  files=[]
  ftemps.each {|ftemp|
    if ftemp[1].empty?
      files += [ sprintf(ftemp[0], $idir, $prefix) ]
    elsif
      files += ftemp[1].collect{|num| sprintf(ftemp[0], $idir, $prefix, num) }
    end
  }
  # ftemps.map {|ftemp|
  #   if ftemp[1].empty?
  #     sprintf(ftemp[0], $idir, $prefix)
  #   elsif
  #     ftemp[1].map{|no| sprintf(ftemp[0], $idir, $prefix, no) }
  #   end
  # }
  return files
end

infile_st    = getInputFiles(file_st)
infile_r1_st = infile_st.select{|f| /r1/ =~ f }
infile_r2_st = infile_st.select{|f| /r2/ =~ f }
infile_r3_st = infile_st.select{|f| /r3/ =~ f }
# infile_r123_st = 
infile_ev    = getInputFiles(file_ev)
infile_r1_ev = infile_ev.select{|f| /r1/ =~ f }
infile_r2_ev = infile_ev.select{|f| /r2/ =~ f }
infile_r3_ev = infile_ev.select{|f| /r3/ =~ f }
# infile_r123_ev = 
infile_ee  = getInputFiles(file_ee)
infile_etc = getInputFiles(file_etc)
infile_enemin2= infile_etc.select{|f| /enemin2/ =~ f}
infile_energy = infile_etc.select{|f| /energy/ =~ f}

# constituent  t, ev, term1, term2, term3
# standrad     t,  z, zv   , za   , ev   , ia, tau , energy, ev_ia

def getOutputFileName(base)
  "#{$odir}/#{base}.gif"
end
def iF2base(infname)
  infname.gsub(/#{$idir}/,"").gsub(/\.dat/,"")
end


#p key_ev  = $options_ev.keys
#p $options_hash["ev"].keys

#def plotWxt(gp, ifname, ofname, xlabel, ylabel, opt, nline)
def plotWxt(gp, sendcmd, ofname, xlabel, ylabel)
  [
   "set output '#{ofname}'",
   "set xlabel '#{xlabel}'",
   "set ylabel '#{ylabel}'",
   #"splot '#{fname}' #{opt} title '#{nline}'"   
   #"plot '#{fname}' #{opt}  title '#{nline}'"
   "plot #{sendcmd}",
   "pause 1.2",
  ].each{|cmd| gp.puts(cmd) }
end

#def plotGif(gp, ifname, ofname, xlabel, ylabel, opt, nline)
def plotGif(gp, sendcmd, ofname, xlabel, ylabel)
  [
   "set terminal gif"      ,
   "set output '#{ofname}'",
   "set xlabel '#{xlabel}'",
   "set ylabel '#{ylabel}'",
   "plot #{sendcmd}"
  ].each{|cmd| gp.puts(cmd) }
  p ofname
end

def sendCmd_mFaG(files, opt)
  # multi files make a graph.
  max = 0
  files.each{|f|
    /saiteki.*_(\d+)\.dat/ =~ f
    max = $1.to_i if max < $1.to_i
  }
  
  #  sendcmd = files.map{|f|
  files.map{|f|
    if /sinc/ =~ f
      "'#{f}' #{opt} ls 1 title 'sin()'"
    elsif /saiteki.*_(\d+)\.dat/ =~ f
      if max != $1.to_i
        "'#{f}' #{opt} ls 2 title ''"
      else
        "'#{f}' #{opt} ls 3 title 'repeat=#{max}'"
      end
    else
      "'#{f}' #{opt}"
    end
  }.join(',')
end

def sendCmds_aFmG(file, which, options_keys)
  # a file make multi graph.
  # options = u #{c1}:#{c2} w l title #{no}

  sendcmds = options_keys.map{|key|
    opt = $options_hash[which][key]['opt']
    "'#{file}' #{opt} title #{key}"}
end

def getRepeatNum(idp_fname)
  if idp_fname =~ /(\d)/
    return "R#{$1}"
  else
    return "sin()"
  end
end

def drawGraph(gp, fname, which)
  #  opt    = $options[which]['opt']
  #  xlabel  =
  #  ylabel = $options[which]['ylabel']
  
  base   = iF2base(fname)
  num    = getRepeatNum(base)
  ofname = "#{$odir}/#{base}_#{which.upcase}.gif"
  
  plotWxt(gp, fname, ofname, "[sec]", ylabel, opt, num)
end


def draw_graph_cmp(gp, fnames, nums, which)
  opt    = $options[which]['opt']
  ylabel = $options[which]['ylabel']
  
  base = iF2base(fnames[0])
  num  = getRepeatNum(base)
  ofname  = "#{$odir}/#{base}_#{which}.gif"
  sendcmd = nums.collect{|num|
    "'#{fnames[num]}' #{opt} title '#{num}'"
  }.join(',')
  
  plotGif(gp, ifname, ofname, xlabel, ylabel, nline)
end

############################################################
Gnuplot.open{|gp|
  [
   "set grid"         ,
   "set size 0.6, 0.6",
  ].each{|cmd| gp.puts(cmd) }

  sendcmd = sendCmd_mFaG( infile_r1_st, $options["st"]["t_z"]['opt'] )
  plotWxt(gp, sendcmd, "hoge",
          $options["st"]["t_z"]["xlabel"],
          $options["st"]["t_z"]["ylabel"] )

  sendcmd = sendCmd_mFaG( infile_r1_st, $options["st"]["t_energy"]['opt'] )
  plotWxt(gp, sendcmd, "hoge",
          $options["st"]["t_energy"]["xlabel"],
          $options["st"]["t_energy"]["ylabel"] )

  
  $options['ev'].keys.each{|k|
    sendcmd = sendCmd_mFaG( infile_r1_ev, $options['ev'][k]['opt'] )
    
    plotGif(gp, sendcmd, $odir+"/"+"r1_ev_"+k+".gif",
            $options['ev'][k]['xlabel'],
            $options['ev'][k]['ylabel'] )
  }

  $options['ev'].keys.each{|k|
    sendcmd = sendCmd_mFaG( infile_r2_ev, $options['ev'][k]['opt'] )
    plotGif(gp, sendcmd, $odir+"/"+"r2_ev_"+k+".gif",
            $options['ev'][k]['xlabel'],
            $options['ev'][k]['ylabel'] )
  }

  $options['ev'].keys.each{|k|
    sendcmd = sendCmd_mFaG( infile_r3_ev, $options['ev'][k]['opt'] )
    
    plotGif(gp, sendcmd, $odir+"/"+"r3_ev_"+k+".gif",
            $options['ev'][k]['xlabel'],
            $options['ev'][k]['ylabel'] )
  }

  exit




  
  $options['st'].keys.each{|k|
    sendcmd = sendCmd_mFaG( infile_r1_st, $options['st'][k]['opt'] )
    
    plotGif(gp, sendcmd, $odir+"/"+"r1_st_"+k+".gif",
            $options['st'][k]['xlabel'],
            $options['st'][k]['ylabel'] )
  }

  $options['st'].keys.each{|k|
    sendcmd = sendCmd_mFaG( infile_r2_st, $options['st'][k]['opt'] )
    plotGif(gp, sendcmd, $odir+"/"+"r2_st_"+k+".gif",
            $options['st'][k]['xlabel'],
            $options['st'][k]['ylabel'] )
  }
  
  $options['st'].keys.each{|k|
    sendcmd = sendCmd_mFaG( infile_r3_st, $options['st'][k]['opt'] )
    plotGif(gp, sendcmd, $odir+"/"+"r3_st_"+k+".gif",
            $options['st'][k]['xlabel'],
            $options['st'][k]['ylabel'] )
  }











  
    
  exit
  
  sendcmd = sendCmd_mFaG( infile_r2_st, $options["st"]["t_z"]['opt'] )
  plotWxt(gp, sendcmd, "hoge",
          $options["st"]["t_z"]["xlabel"],
          $options["st"]["t_z"]["ylabel"] )
  
  sendcmd = sendCmd_mFaG( infile_r3_st, $options["st"]["t_energy"]['opt'] )
  plotWxt(gp, sendcmd, "hoge",
          $options["st"]["t_energy"]["xlabel"],
          $options["st"]["t_energy"]["ylabel"] )


  sendcmd = sendCmd_mFaG(infile_energy, $options["etc"]["t_e1e2e3"]["opt"])
  plotWxt(gp, sendcmd, "hoge",
          $options["etc"]["t_e1e2e3"]["xlabel"],
          $options["etc"]["t_e1e2e3"]["ylabel"])


  sendcmd = sendCmd_mFaG(infile_enemin2, $options['etc']['n_e1e2e3']['opt'])
  plotWxt(gp, sendcmd, "hoge",
          $options['etc']['n_e1e2e3']['xlabel'],
          $options['etc']['n_e1e2e3']['ylabel'])
  

  exit
}
