set terminal png
set pointsize 1
set output "phase_lineout.png"

plot "phase_lineout_extract.dat" i 0 lt 1, \
     "phase_lineout_extract.dat" i 1 lt 2, \
     "phase_lineout_extract.dat" i 2 lt 3, \
     "phase_lineout_extract.dat" i 3 lt 4, \
     "phase1.dat" u ($1*1.e6):($2) lt 6, \
     "phase2.dat" u ($1*1.e6):($2) lt 7, \
     "phase3.dat" u ($1*1.e6):($2) lt 8, \
     "phase4.dat" u ($1*1.e6):($2) lt 9
