# Datafile columns
# "Load (N)" "Displacement (mm)"

span = 80
w = 40.0
b = 19.6

f(x) = m*x + n

set title "Foam sandwich panel"
set xlabel "Strain (1/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0.00004:0.00014]
fit f(x) "Sandwich_foamed.dat" u ($2/span**2):(3*$1*span/2/w/b**2) via m,n
set xrange [-0.0001:0.0006]
set yrange [0:40]
set key top left

plot \
"Sandwich_foamed.dat" u ($2/span**2):(3*$1*span/2/w/b**2) w l lw 2 title "Data", \
f(x) w l dt 2 title sprintf("y = %3.2fx + %3.2f", m, n)
