# plot_curve.gnuplot
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'curve_plot.png'

set title "Implicitly Defined Curve"
set xlabel "x"
set ylabel "y"
set grid

plot 'curve_points.txt' using 1:2 with linespoints title 'Curve'
