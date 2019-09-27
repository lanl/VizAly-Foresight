# Gnuplot script for entropy analysis

reset
set terminal postscript eps enhanced color 14 size 24cm, 9cm
set output "spectrum.eps"
set multiplot layout 1,3 title 'halo-based power spectrum analysis - m000.full.mpicosmo.499 - 1073726359 particles (2% halos)'

# next colors '#FF00FF' and '#A0522D'

# ---------------
set title ''
set size ratio 1
set xlabel "k (h/Mpc)"
set ylabel "primordial power spectrum P_0(k) (Mpc/h)^3"
set key Left reverse below maxcols 2
set grid
set logscale x 10
set logscale y 10
set xrange [0.01:10]

plot 'data/pk-full-raw.dat'      using 1:2:3 title 'full-raw'     w errorbars lc rgb '#FF00FF',\
     'data/pk-merged-20bits.dat' using 1:2:3 title 'lossy-20bits' w errorbars lc rgb '#000000',\
     'data/pk-merged-24bits.dat' using 1:2:3 title 'lossy-24bits' w errorbars lc rgb '#0000FF',\
     'data/pk-merged-26bits.dat' using 1:2:3 title 'lossy-26bits' w errorbars lc rgb '#CB0707',\
     'data/pk-sampled-0.80.dat'  using 1:2:3 title 'sampled-0.80' w errorbars lc rgb '#800080',\
     'data/pk-sampled-0.90.dat'  using 1:2:3 title 'sampled-0.90' w errorbars lc rgb '#000080'
#     'data/pk-sampled-0.95.dat'  using 1:2:3 title 'sampled-0.95' w errorbars lc rgb '#006400'

# --------------
reset
unset logscale
unset key
set title 'compression level'
set size ratio 0.8
set xlabel 'scalars'
set ylabel 'gain'
set yrange [0:8]
set key right maxcols 2
set boxwidth 1
set style fill solid 2.00
#
set grid
set style data histograms
plot 'stats.dat' using 2:xtic(1) lt rgb '#000000' title 'lossy-20bits',\
     'stats.dat' using 3:xtic(1) lt rgb '#0000FF' title 'lossy-25bits',\
     'stats.dat' using 4:xtic(1) lt rgb '#CB0707' title 'lossy-26bits',\
     'stats.dat' using 5:xtic(1) lt rgb '#800080' title 'sampled-0.80',\
     'stats.dat' using 6:xtic(1) lt rgb '#800080' title 'sampled-0.90',\
     'stats.dat' using 7:xtic(1) lt rgb '#006400' title 'sampled-0.95'
# ---------------
reset
set title ''
set xlabel "wavenumber k (h/Mpc)"
set size ratio 1
#set ylabel "P_{zip}/P_{raw}"
set ylabel "discrepancy"
set key Left reverse below maxcols 2

set grid
unset logscale

set xrange [0:10]
#set yrange [0.98:1.02]

plot 'data/ratio.dat' using 1:($ 9/$2) title 'lossy-20bits' w lines lc rgb '#000000',\
     'data/ratio.dat' using 1:($10/$2) title 'lossy-24bits' w lines lc rgb '#0000FF',\
     'data/ratio.dat' using 1:($11/$2) title 'lossy-26bits' w lines lc rgb '#CB0707',\
     'data/ratio.dat' using 1:($12/$2) title 'sampled-0.80'  w lines lc rgb '#800080',\
     'data/ratio.dat' using 1:($13/$2) title 'sampled-0.90'  w lines lc rgb '#000080'
#     'data/ratio.dat' using 1:($14/$2) title 'sampled-0.95'  w lines lc rgb '#006400'

unset multiplot
