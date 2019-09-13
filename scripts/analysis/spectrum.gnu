# Gnuplot script for entropy analysis

reset
set terminal postscript eps enhanced color 14 size 16cm, 9cm
set output "spectrum.eps"
set multiplot layout 1,2 title 'halo-based power spectrum analysis - m000.full.mpicosmo.499 - 1073726359 particles (2% halos)'

# ---------------
set title ''
set size ratio 1
set xlabel "k (h/Mpc)"
#set ylabel "P_0(k) (Mpc/h)^3"
set ylabel "power spectrum P(k)"
set key Left reverse below maxcols 2
set grid
set logscale x 10
set logscale y 10
set xrange [0.01:10]

plot 'pk-full-raw.dat'            using 1:2:3 title 'full-raw'      w errorbars lc rgb '#FF00FF',\
     'pk-full-sz.dat'             using 1:2:3 title 'full-sz'       w errorbars lc rgb '#800080',\
     'pk-full-fpzip.dat'          using 1:2 title   'full-fpzip'    w lines     lc rgb '#000080',\
     'pk-combined-zip-20bits.dat' using 1:2:3 title 'merged-20bits' w errorbars lc rgb '#000000',\
     'pk-combined-zip-22bits.dat' using 1:2:3 title 'merged-22bits' w errorbars lc rgb '#006400',\
     'pk-combined-zip-24bits.dat' using 1:2:3 title 'merged-24bits' w errorbars lc rgb '#0000FF',\
     'pk-combined-zip.dat'        using 1:2:3 title 'merged-26bits' w errorbars lc rgb '#CB0707'
#      'pk-combined-zip-25bits.dat' using 1:2:3 title 'merged-25bits' w errorbars lc rgb '#A0522D',\

# --------------
#unset logscale
#unset key
#set title 'compression level'
#set size ratio 0.7
#set xlabel 'scalars'
#set ylabel 'gain'
#set yrange [0:6]
#set key right maxcols 1
#set boxwidth 1
#set style fill solid 1.00
#
#set grid
#set style data histograms
#plot 'metrics.dat' using 2:xtic(1) lt rgb '#800080' title 'halo [sz], abs:0.003',\
#     'metrics.dat' using 3:xtic(1) lt rgb '#CB0707' title 'non-halo [fpzip], bits: 26'
# ---------------

set title ''
set xlabel "k (h/Mpc)"
#set ylabel "P_{zip}/P_{raw}"
set ylabel "discrepancy"
set key Left reverse below maxcols 2

set grid
unset logscale

set xrange [0:10]
set yrange [0.98:1.02]

plot 'ratio.dat' using 1:($10/$2) title 'merged-20bits' w lines lc rgb '#000000',\
     'ratio.dat' using 1:($11/$2) title 'merged-22bits'w lines lc rgb '#006400',\
     'ratio.dat' using 1:($12/$2) title 'merged-24bits'w lines lc rgb '#0000FF',\
     'ratio.dat' using 1:($9/$2) title 'merged-26bits'w lines lc rgb '#CB0707',\
     'ratio.dat' using 1:( $4/$2) title 'full-sz'      w lines lc rgb '#800080',\
     'ratio.dat' using 1:( $5/$2) title 'full-fpzip'   w lines lc rgb '#000080'

#
#     'ratio.dat' using 1:( $9/$2) title 'merged-raw' w lines lc rgb '#FF00FF',\
#     'ratio.dat' using 1:($10/$2) title 'merged-25bits'w lines lc rgb '#A0522D',\
unset multiplot
