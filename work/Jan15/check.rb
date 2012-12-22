#!/usr/bin/ruby
require 'rubygems'
require 'gnuplot'

a_flag = [ 1, 2, 3, 4, 5, 6, 7, 8, 9 ]

file_st  = [["%s/%s_sinc_st.dat"           , []     ],
            ["%s/%s_r1_sinc_st.dat"        , []     ],
            ["%s/%s_r2_sinc_st.dat"        , []     ],
            ["%s/%s_r3_sinc_st.dat"        , []     ],
            ["%s/%s_saiteki_st_%05d.dat"   , a_flag ],
            ["%s/%s_r1_saiteki_st_%05d.dat", a_flag ],
            ["%s/%s_r2_saiteki_st_%05d.dat", a_flag ],
            ["%s/%s_r3_saiteki_st_%05d.dat", a_flag ],   ]
file_ev  = [["%s/%s_sinc_ev.dat"           , []     ],
            ["%s/%s_r1_sinc_ev.dat"        , []     ],
            ["%s/%s_r2_sinc_ev.dat"        , []     ],
            ["%s/%s_r3_sinc_ev.dat"        , []     ],
            ["%s/%s_saiteki_ev_%05d.dat"   , a_flag ],
            ["%s/%s_r1_saiteki_ev_%05d.dat", a_flag ],
            ["%s/%s_r2_saiteki_ev_%05d.dat", a_flag ],
            ["%s/%s_r3_saiteki_ev_%05d.dat", a_flag ],   ]
file_ee  = [["%s/%s_sinc_endeffector.dat"        , []    ],
            ["%s/%s_saiteki_endeffector_%05d.dat", a_flag],]
file_etc = [["%s/%s_enemin2.dat"                 , []    ],
            ["%s/%s_energy_sinc.dat"             , []    ],
            ["%s/%s_energy_saiteki_%05d.dat"     , a_flag],]

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
