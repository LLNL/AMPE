set terminal png
set size 0.7,0.7

set datafile separator ","

set output "Energy.png"
set xlabel "Time"
set ylabel "Free Energy"
set xtics 500
set ytics 0.5e6

set xrange [0:1500]
set key off

plot 'free_energy.csv' u 1:2 w lines



