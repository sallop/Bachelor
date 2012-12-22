#!/usr/bin/ruby

require 'rubygems'
require 'gnuplot'

class DataFile
  @@IDir="dir-etc"
  @@ODir="dir-zcurve/pderiv"
  def initialize(fname, title, xlabel, ylabel, option)
    @fname   = fname
    @title   = title
    @xlabel  = xlabel
    @ylabel  = ylabel
    @option  = option
  end
  def giffile
    @fname+".gif"
  end
  def datfile
    @fname+".dat"
  end  
end

h = {
  "title"=>"zz",
  "xlabel"=>"[sec]",
  "ylabel"=>"[m/sec]",
}

files = [ DataFile.new("zz", "[sec]", "[m]", "", ""),
          DataFile.new("zv", "[sec]", "[m/sec]", "",""),
          DataFile.new("za", "[sec]", "[m/sec^2]", "", "")]

files.each do |f|
  print f.giffile()+"\n"
  print f.datfile()+"\n"  
end


Gnuplot.open do |gp|
  Gnuplot::Plot.new(gp) do |plot|
    plot.title "Array Plot Example"
    plot.ylabel "x"
    plot.xlabel "x^2"

    x = (0..50).collect {|v| v.to_f }
    y = x.collect {|v| v**2 }

    plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
      ds.with = "linespoints"
      ds.notitle
    end
  end
end


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
