set terminal png
set pointsize 1
set output "dtvst.png"

set logscale y

plot "dtvst.dat" lt 1 pt 5 with linespoints

