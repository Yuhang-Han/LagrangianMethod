# Requires gnuplot 5.4 or newer and the pngcairo terminal.
#
# Usage:
#   gnuplot -c plot_quad_mesh.gp input.csv [output.png]

if (ARGC < 1) {
    print "Usage: gnuplot -c plot_quad_mesh.gp input.csv [output.png]"
    exit
}

csv_file = ARG1
png_file = (ARGC >= 2) ? ARG2 : csv_file . ".png"

set datafile separator ","
set datafile columnheaders

# Column 1 contains the same time on every node record.
stats csv_file using 1 name "Time" nooutput

set terminal pngcairo size 1000,800 enhanced font "Sans,12"
set output png_file

set title sprintf("Quad mesh at t = %.10g", Time_min)
set xlabel "x"
set ylabel "y"
unset key
unset colorbox
unset ztics

# View the fictitious z=0 surface from above.
set view map
set size ratio -1
set zrange [-1:1]

# Columns 5 and 6 are x and y.
# Each blank-line-separated mesh row is one scan.
# Equal-length scans connect corresponding nodes in both directions.
splot csv_file using 5:6:(0.0) with lines \
    linecolor rgb "#244a73" linewidth 1.2

unset output
