#!/usr/bin/gnuplot 

set term post enhanced color eps 20 ''
set output 'h2.eps'

set key font ",20"
set key box top right

#set format y "%1.1f"
set format x "%1.1f"

set ylabel "ERPA-FCI difference (mHa)"
set xlabel "R_{O-H}  (a_0)"

set style line 1 lc rgb 'red' lt 1 lw 4 pt 5 ps 1.3
set style line 2 lc rgb 'dark-green' lt 1 lw 2 pt 6 ps 1.3
set style line 3 lc rgb 'blue' lt 1 lw 4 pt 4 ps 1.3
set style line 4 lc rgb 'black' lt 1 lw 4 pt 3 ps 1.3

plot "h2.data3" u 1 : (1000*$2)  with linespoints ls 1 title "ERPA-FCI difference"
#     "h2.data" u (1.889725989 * $1) : 5  with linespoints ls 2 title "AC corrected-CASSCF(NO basics)",\
 #    "h2.data" u (1.889725989 * $1) : 3  with linespoints ls 3 title "AC corrected-CASSCF(MO basics)",\
  #   "h2.data" u (1.889725989 * $1) : 4  with linespoints ls 4 title "CASPT2"
