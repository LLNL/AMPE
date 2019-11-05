set terminal png
set size 0.7,0.7

set datafile separator ","

set output "SolidFraction.png"
set xlabel "Time"
set ylabel "Solid Fraction"
set xtics 500
set ytics 0.005

set xrange [0:1500]
set key off

plot 'solid_fraction.csv' u 1:2 w lines



