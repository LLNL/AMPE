set terminal png
set size 0.7,0.7

#set datafile separator ","

set output "phase_field_1500.png"
set xlabel "x"
set ylabel "y"
set xtics 100
set ytics 100

set key off

plot 'phase_field_1500.csv' u 1:2 w lines



