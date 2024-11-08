set terminal svg standalone fname 'Arial bold' fsize 20
set pointsize 3
set datafile separator ","
set key autotitle columnhead

set output "CvsT.svg"
set xlabel "Temperature (K)"
set ylabel "Composition (\%)"

plot 'CvsTliquid.csv' u 1:2 t 'liquid' w lines lw 2, \
     'CvsTsolid.csv'  u 1:2 t 'solid' w lines lw 2
