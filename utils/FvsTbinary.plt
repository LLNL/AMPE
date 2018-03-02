set terminal png
#set pointsize 2
set size 0.7,0.7

set output "FvsT.png"
set xlabel "Temperature (K)"
set ylabel "Energy (J/mol)"

plot "FsC0vsT.dat" u ($1):($2) title 'FsC0' w lines lt 1, \
     "FsC1vsT.dat" u ($1):($2) title 'FsC1' w lines lt 2, \
     "FlC0vsT.dat" u ($1):($2) title 'FlC0' w lines lt 4, \
     "FlC1vsT.dat" u ($1):($2) title 'FlC1' w lines lt 5


