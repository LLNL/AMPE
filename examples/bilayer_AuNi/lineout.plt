set terminal png
set pointsize 1
set output "lineout.png"

plot "lineout_extract.dat" i 0 w lines lt 1, \
     "lineout_extract.dat" i 1 w lines lt 2, \
     "lineout_extract.dat" i 2 w lines lt 3, \
     "lineout_extract.dat" i 3 w lines lt 4, \
     "lineout_extract.dat" i 4 w lines lt 5, \
     "lineout_extract.dat" i 5 w lines lt 6, \
     "lineout_extract.dat" i 6 w lines lt 7, \
     "lineout_extract.dat" i 7 w lines lt 8

