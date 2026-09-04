set datafile separator ","

set xlabel "(Eulerian) x"
set ylabel "Pressure"

set border linewidth 1.5
set tics scale 0.7
set grid ytics

set key top right
set xrange [*:*]
set yrange [*:*]

# ----------------- PDF ------------------------
set terminal pdfcairo enhanced font "Times New Roman,16" size 5in,3.5in
set output "pressure.pdf"

plot "build/output/cellCentered_N1000_P.csv" using 1:2 \
	with lines lw 2 \
	title "Pressure"

set output

# ----------------- PNG ------------------------
set terminal pngcairo size 1200,800 enhanced font "Times New Roman,16"
set output "pressure.png"

plot "build/output/cellCentered_N1000_P.csv" using 1:2 \
	with lines lw 2 \
	title "Pressure"

set output
