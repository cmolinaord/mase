# Datafile columns
# "Load (N)" "Displacement (mm)"

span = 70.4
w = 10.5
b = 4.4

f(x) = m*x + n

set title "Reinforced PP with 30% short GF (3)"
set xlabel "Strain (1/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0:0.0002]
fit f(x) "reinf_PP3.dat" u ($2/span**2):(3*$1*span/2/w/b**2) via m,n
set xrange [0:0.0016]
set yrange [0:150]
set key top left

plot \
"reinf_PP3.dat" u ($2/span**2):(3*$1*span/2/w/b**2) w l lw 2 title "Data", \
f(x) w l dt 2 title sprintf("y = %3.2fx + %3.2f", m, n)
