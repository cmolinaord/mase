set term pngcairo enhanced size 1000,600
set key top left
##################################
set output "neat_PP1.png"

span = 68.8
w = 10.5
b = 4.3

f(x) = m*x + n

set title "Unreinforced PP(1)"
set xlabel "Strain (1/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0:0.0005]
fit f(x) "neat_PP1.dat" u ($2/span**2):(3*$1*span/2/w/b**2) via m,n
set xrange [0:0.003]
set yrange [0:50]
set key top left

plot \
"neat_PP1.dat" u ($2/span**2):(3*$1*span/2/w/b**2) w l lw 2 title "Data", \
f(x) w l dt 2 title sprintf("y = %3.2fx + %3.2f", m, n)

#################
set output "neat_PP2.png"

span = 70.4
w = 10.5
b = 4.4

f(x) = m*x + n

set title "Unreinforced PP(2)"
set xlabel "Strain (1/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0:0.0005]
fit f(x) "neat_PP2.dat" u ($2/span**2):(3*$1*span/2/w/b**2) via m,n
set xrange [0:0.0032]
set yrange [0:50]
set key top left

plot \
"neat_PP2.dat" u ($2/span**2):(3*$1*span/2/w/b**2) w l lw 2 title "Data", \
f(x) w l dt 2 title sprintf("y = %3.2fx + %3.2f", m, n)

#################
set output "neat_PP3.png"

span = 70.4
w = 10.5
b = 4.4

f(x) = m*x + n

set title "Unreinforced PP(3)"
set xlabel "Strain (1/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0:0.0005]
fit f(x) "neat_PP3.dat" u ($2/span**2):(3*$1*span/2/w/b**2) via m,n
set xrange [0:0.0034]
set yrange [0:50]
set key top left

plot \
"neat_PP3.dat" u ($2/span**2):(3*$1*span/2/w/b**2) w l lw 2 title "Data", \
f(x) w l dt 2 title sprintf("y = %3.2fx + %3.2f", m, n)

#################
set output "reinf_PP1.png"

span = 70.4
w = 10.5
b = 4.4

f(x) = m*x + n

set title "Reinforced PP with 30% short GF (1)"
set xlabel "Strain (1/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0:0.0002]
fit f(x) "reinf_PP1.dat" u ($2/span**2):(3*$1*span/2/w/b**2) via m,n
set xrange [0:0.0016]
set yrange [0:150]
set key top left

plot \
"reinf_PP1.dat" u ($2/span**2):(3*$1*span/2/w/b**2) w l lw 2 title "Data", \
f(x) w l dt 2 title sprintf("y = %3.2fx + %3.2f", m, n)

#################
set output "reinf_PP2.png"

span = 70.4
w = 10.5
b = 4.4

f(x) = m*x + n

set title "Reinforced PP with 30% short GF (2)"
set xlabel "Strain (1/mm)"
set ylabel "Stress (GPa)"

set grid
set xrange [0:0.0002]
fit f(x) "reinf_PP2.dat" u ($2/span**2):(3*$1*span/2/w/b**2) via m,n
set xrange [0:0.0016]
set yrange [0:150]
set key top left

plot \
"reinf_PP2.dat" u ($2/span**2):(3*$1*span/2/w/b**2) w l lw 2 title "Data", \
f(x) w l dt 2 title sprintf("y = %3.2fx + %3.2f", m, n)

#################
set output "reinf_PP3.png"

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

#################
set output "foam_sandwich.png"

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

#################
set output "cardboard_sandwich.png"

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
