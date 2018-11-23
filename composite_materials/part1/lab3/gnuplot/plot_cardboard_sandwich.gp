# Datafile columns
# "Load (N)" "Displacement (mm)"

span = 114
w = 57
b = 30.7

f(x) = m*x + n

set title "Cardboard sandwich"
set xlabel "Strain (1/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0.000005:0.00004]
fit f(x) "Sandwich_cardboard.dat" u ($2/span**2):(3*$1*span/2/w/b**2) via m,n
set xrange [-0.0001:0.0013]
set yrange [0:1.2]
set key top left

plot \
"Sandwich_cardboard.dat" u ($2/span**2):(3*$1*span/2/w/b**2) w l lw 2 title "Data", \
f(x) w l dt 2 title sprintf("y = %3.2fx + %3.2f", m, n)
