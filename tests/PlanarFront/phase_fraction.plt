set terminal png
set pointsize 1
set output "volfraction.png"

plot "volfraction.dat" lt 1 pt 5 with linespoints

