set terminal png
set pointsize 1
set output "lineout.png"

plot "lineout_extract.dat" i 0 lt 1, \
     "lineout_extract.dat" i 1 lt 2, \
     "lineout_extract.dat" i 2 lt 3, \
     "lineout_extract.dat" i 3 lt 4, \
     "lineout_extract.dat" i 4 lt 5, \
     "temperature1.dat" lt 6, \
     "temperature5.dat" lt 7
