set terminal pngcairo enhanced
set pointsize 1.7
set output "lineout.png"

set yrange [0.35:0.65]
set xrange [1.:1.12]

ceql=0.604113
ceqs=0.457175
c0=0.457175
v=24686.
DL=1.e3
x0=1.051772928

lambda=0.0105

c(x)  = (x>0.) ? c0+(ceql-c0)*exp(-v*x/DL) : ceqs

set xlabel "position in {/Symbol m}m"
set ylabel "mole fraction of Cu"
set ytics 0.05 nomirror
set y2tics 0.1 nomirror 
set y2label "phase variable" 

plot "lineout_extract.dat" u ($1):($3) w lines lt 2 lw 2 title "c", \
     "lineout_extract.dat" u ($1):($4) w lines lt 3 lw 2 title "cl", \
     "lineout_extract.dat" u ($1):($5) w lines lt 4 lw 2 title "cs", \
     ceql w lines lt 5, \
     ceqs w lines lt 6, \
     c(x-x0) w lines lt 7 lw 2 title "c sharp", \
     "lineout_extract.dat" u ($1):($2) w lines lt 1 lw 2 title "phi" axes x1y2

