set terminal png
set pointsize 2
set size 0.7,0.7

set datafile separator ","

set output "TipPosition.png"
set xlabel "Time"
set ylabel "Tip position"
set xtics 500
set ytics 100

set xrange [0:1500]
set key off

plot 'tip_position.csv' u 1:2 w lines lw 2



