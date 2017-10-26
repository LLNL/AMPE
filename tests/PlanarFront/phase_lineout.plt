set terminal png
set pointsize 1
set output "phase_lineout.png"

plot "phase_lineout_extract.dat" i 0 lt 1, \
     "phase_lineout_extract.dat" i 1 lt 2, \
     "phase_lineout_extract.dat" i 2 lt 3, \
     "phase_lineout_extract.dat" i 3 lt 4, \
     "phase_lineout_extract.dat" i 4 lt 5, \
     "phase_lineout_extract.dat" i 5 lt 6, \
     "phase1.dat" lt 8, \
     "phase3.dat" lt 8, \
     "phase6.dat" lt 8
