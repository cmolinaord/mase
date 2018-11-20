# Datafile columns
# "Load (N)" "Displacement (mm)"

l0 = 80
a = 4.45
b = 10.4


f(x) = m*x + n

set title "PP + 30\% weight Glass fiber"
set xlabel "Strain (mm/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0:0.025]
fit f(x) "PP1.dat" u ($2/l0):($1/a/b) via m,n
set xrange [0:0.08]
set yrange [0:100]

plot \
"PP1.dat" u ($2/l0):($1/a/b) w l lw 2 title "Data", \
f(x) w l dt 2 title "Linear fit"


set term pngcairo enhanced size 1000,600
set output "PP1.png"

replot
