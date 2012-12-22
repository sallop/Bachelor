#!/usr/bin/ruby

require 'rubygems'
require 'gnuplot'

#class FileData
class FD
  @IDIR = "dir-etc"
  @ODIR = "dir-zcurve/pderiv"
  def initialize(fname, xlabel, ylabel, title, option)
    @fname  = fname
    @xlabel = xlabel
    @ylabel = ylabel
    @title  = title
    @option = option
  end

  def setOutputData()
    #
  end
  def outFname()
    #
  end
  def inFname()
    #
  end
end

def fstat(fname, xlabel, ylabel, title, option)
  { "fname"   => fname,
    "xlabel" => xlabel,
    "ylabel" => ylabel,
    "title"  => title ,
    "option" => option }
end

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

def f2mat(fname)
  mat = []
  if File.exist?(fname)
    File::open(fname, "r"){|f|
      f.each{|line|
        nums = line.split().map{|str| str.to_f}
        mat += [nums]
      }
    }
  elsif
    puts "file not exist."
  end
  
  return mat
  #volt_consist [t, ev, term1, term2, term3]
  #standard     [t, z, zv, za, ev, ia, tau, energy, ev_ia]
  #endeffector  [t, r3x, r3y, r3z, r3xv, r3yv, r3zv, r3xa, r3ya, r3za]
end

mat_ee = Hash.new
file_ee.each do |ftemp, nums|
  if nums.empty?
    fname = sprintf(ftemp, "dir-etc", "zcurve")
    mat_ee[fname] = [ f2mat(fname) ]
  elsif
    nums.each{|num|
      fname = sprintf(ftemp, "dir-etc", "zcurve", num)
      mat_ee[fname] = [ f2mat(fname) ]
      #mat_ee += [ f2mat(fname) ]
    }
  end
end

pdata = [ FD.new("trajectory", "[m]"  , "[m]"    , "trjc", ""),
          FD.new("r3xv"      , "[sec]", "[m/sec]", "r3xv", ""),
          FD.new("r3yv"      , "[sec]", "[m/sec]", "r3yv", ""),
          FD.new("r3zv"      , "[sec]", "[m/sec]", "r3zv", ""),]

mat_ee.each do |fname, mat| # |key, value|
  ts, r3xs, r3ys, r3zs, r3xvs, r3yvs, r3zvs = [],[],[],[],[],[],[]
  mat.each {|row|
    row.each{|t, r3x, r3y, r3z, r3xv, r3yv, r3zv, r3xa, r3ya, r3za|
      ts    += [t]
      r3xs  += [r3x]
      r3ys  += [r3y]
      r3zs  += [r3z]
      r3xvs += [r3xv]
      r3yvs += [r3yv]
      r3zvs += [r3zv]
    }
  }
  
  Gnuplot.open{|gp|
    gp.puts("set grid")
    gp.puts("set size 0.6, 0.6\n")
    gifname = fname.gsub(/\.dat$/, ".gif")
    puts gifname

    #   gp.puts("set output '#{gifname}'\n", fname.gsub("\.dat$", ".gif"))
    #    if fname =~ /trajectory/
    if true
      gp.puts "set xlabel '[m]'"
      gp.puts "set ylabel '[m]'"
      gp.puts "set zlabel '[m]'"
      gp.puts "splot '#{fname}' u 2:3:4 w l"
    else
      gp.puts "set xlabel '[m]'"
      gp.puts "set ylabel '[m]'"
      gp.puts "plot '#{fname}' u 1:2 w l"
    end
    gp.puts "pause 1.8"
  }
end
# Gnuplot

# plot data

p pdata
exit

exist_lst   = []
nothing_lst = []

file_st.each do |files|
  if files[1].empty?
    fname = sprintf(files[0], "dir-etc", "zcurve")
    if File::exist?(fname)
      printf("exist #{fname}\n")
    elsif
      nothing_lst += [fname]
      printf("Nothing #{fname}\n")
    end
  elsif
    files[1].each {|num|
      fname = sprintf(files[0], "dir-etc", "zcurve", num)
      if File::exist?(fname)
        printf("exist #{fname}\n")
      elsif
        nothing_lst += [fname]        
        printf("Nothing #{fname}\n")
      end
    }
  end
end

file_ev.each do |files|
  if files[1].empty?
    fname = sprintf(files[0], "dir-etc", "zcurve")
    if File::exist?(fname)
      printf("exist #{fname}\n")
    elsif
      nothing_lst += [fname]      
      printf("Nothing #{fname}\n")      
    end
  elsif
    files[1].each {|num|
      fname = sprintf(files[0], "dir-etc", "zcurve", num)
      if File::exist?(fname)
        printf("exist #{fname}\n")
      elsif
        nothing_lst += [fname]        
        printf("Nothing #{fname}\n")
      end
    }
  end
end

file_ee.each do |files|
  if files[1].empty?
    fname = sprintf(files[0], "dir-etc", "zcurve")
    if File::exist?(fname)
      printf("exist #{fname}\n")
    elsif
      nothing_lst += [fname]
      printf("Nothing #{fname}\n")      
    end
  elsif
    files[1].each {|num|
      fname = sprintf(files[0], "dir-etc", "zcurve", num)
      if File::exist?(fname)
        printf("exist #{fname}\n")
      elsif
        nothing_lst += [fname]        
        printf("Nothing #{fname}\n")
      end
    }
  end
end


file_etc.each do |files|
  if files[1].empty?
    fname = sprintf(files[0], "dir-etc", "zcurve")
    if File::exist?(fname)
      printf("exist #{fname}\n")
    elsif
      nothing_lst += [fname]      
      printf("nothing #{fname}\n")      
    end
  elsif
    files[1].each {|num|
      fname = sprintf(files[0], "dir-etc", "zcurve", num)
      if File::exist?(fname)
        printf("exist #{fname}\n")
      elsif
        nothing_lst += [fname]        
        printf("nothing #{fname}\n")
      end
    }
  end
end

nothing_lst.each{|no| puts "#{no} is not found."}


#### file_ev.each do |files|
####   if files[1].empty?
####     fname = sprintf(files[0], "dir-etc", "zcurve")
####     mat = f2_volt_consist(fname)
####     next if mat.empty?
####     
####     Gnuplot.open do |gp|
####       gp.printf("set grid\n")
####       gp.printf("set size 0.6, 0.6\n")
####       Gnuplot::Plot.new(gp) do |plot|
####         plot.title "#{fname}"
####         plot.ylabel "[sec]"
####         plot.xlabel "[]"
####         xs, ys = [], []
####         mat.each {|x, y| xs += [x]; ys += [y] }
####         plot.data << Gnuplot::DataSet.new( [xs, ys] ) do |ds|
####           #ds.with = "linespoints"
####           ds.with = "lines"
####           ds.notitle #ds.title = ""
####         end
####         gp.printf("pause 1.8\n")
####       end
####     end
####     
####   elsif
####     files[1].each{|num|
####       fname = sprintf(files[0], "dir-etc", "zcurve", num)
####       mat = f2_volt_consist(fname)
####       next if mat.empty?
####       
####       Gnuplot.open do |gp|
####         gp.printf("set grid\n")
####         gp.printf("set size 0.6, 0.6\n")
####         Gnuplot::Plot.new(gp) do |plot|
####           plot.title "#{fname}"
####           plot.ylabel "[sec]"
####           plot.xlabel "[]"
####           xs, ys = [], []
####           mat.each {|x, y| xs += [x]; ys += [y] }
####           plot.data << Gnuplot::DataSet.new( [xs, ys] ) do |ds|
####             #ds.with = "linespoints"
####             ds.with = "lines"
####             ds.notitle #ds.title = ""
####           end
####           gp.printf("pause 1.8\n")
####         end
####       end
####     }    
####   end
#### end

#fi.close

io = IO.popen("gnuplot", "w")
## files1.each do |f|
##   puts f
##   str = sprintf(f, 'dir-etc', 'zcurve')
##   io.printf("set grid\n")
##   io.printf("plot '%s' u 1:2 w l\n", str)
##   io.puts "pause 1.8\n"
## end
## 
## files2.each do |f|
##   puts f
##   animation_flag.each {|num|
##     str = sprintf(f, 'dir-etc', 'zcurve', num)
##     io.printf("set grid\n")
##     io.printf("plot '%s' u 1:2 w l\n", str)
##     io.puts "pause 1.8\n"
##   }
## end
io.close

### Gnuplot.open do |gp|
###   Gnuplot::Plot.new(gp) do |plot|
###     plot.title "Array Plot Example"
###     plot.ylabel "x"
###     plot.xlabel "x^2"
### 
###     x = (0..50).collect {|v| v.to_f }
###     y = x.collect {|v| v**2 }
### 
###     plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
###       ds.with = "linespoints"
###       ds.notitle
###     end
###   end
### end

#File.open("gnuplot.dat","w") do |gp|
#  Gnuplot::Plot.new( gp ) do |plot|
#    plot.xrange "[-10:10]"
#    plot.title "Sin Wave Example"
#    plot.ylabel "x"
#    plot.xlabel "sin(x)"
#
#    x = (0..50).collect {|v| v.to_f }
#    y = x.collect {|v| v**2 }
#
#    plot.data = [
#                 Gnuplot::DataSet.new("sin(x)"){|ds|
#                   ds.with = "lines"
#                   ds.title = "String function"
#                   ds.linewidth = 4
#                 },
#                 Gnuplot::DataSet.new([x, y]){|ds|
#                   ds.with = "linespoints"
#                   ds.title = "Array data"
#                 }
#                ]
#  end
#end

#### Gnuplot::Plot.new(gp){|plot|
      ####   plot.title "#{fname}"
      ####   plot.xlabel "[sec]"
      ####   plot.ylabel "[m/sec]"
      ####   plot.data << Gnuplot::DataSet.new( [ts, r3xvs] ) do |ds|
      ####     #ds.with = "linespoints"
      ####     ds.with = "lines"
      ####     ds.notitle #ds.title = ""
      ####   end
      ####   gp.printf("pause 1.8\n")
      #### }
