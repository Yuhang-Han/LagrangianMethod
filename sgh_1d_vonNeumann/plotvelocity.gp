set datafile separator ","

set xlabel "(Eulerian) x"
set ylabel "Velocity"

set border linewidth 1.5
set tics scale 0.7
set grid ytics

set key top right
set xrange [*:*]
set yrange [*:*]

# ------------ PDF ----------------

set terminal pdfcairo enhanced font "Times New Roman,16" size 5in,3.5in
set output "velocity.pdf"

plot "build/output/staggeredGrid_M1000_N1000_U.csv" using 1:2 \
	with lines lw 2 \
	title "Velocity"

set output

# ------------ PNG ----------------

set terminal pngcairo size 1200,800 enhanced font "Times New Roman,16"
set output "velocity.png"

plot "build/output/staggeredGrid_M1000_N1000_U.csv" using 1:2 \
	with lines lw 2 \
	title "Velocity"

set output
