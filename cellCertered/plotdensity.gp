set datafile separator ","

set xlabel "(Eulerian) x"
set ylabel "Density"

set border linewidth 1.5
set tics scale 0.7
set grid ytics

set key top right
set xrange [*:*]
set yrange [0:*]

# ------------ PDF ----------------

set terminal pdfcairo enhanced font "Times New Roman,16" size 5in,3.5in
set output "density.pdf"

plot "build/output/cellCentered_N1000_Rho.csv" using 1:2 \
	with lines lw 2 \
	title "Density"

set output

# ------------ PNG ----------------

set terminal pngcairo size 1200,800 enhanced font "Times New Roman,16"
set output "density.png"

plot "build/output/cellCentered_N1000_Rho.csv" using 1:2 \
	with lines lw 2 \
	title "Density"

set output
