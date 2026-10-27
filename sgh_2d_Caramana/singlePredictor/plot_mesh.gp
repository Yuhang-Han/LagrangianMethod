set terminal pngcairo size 1200,1200
set output "build/output/mesh_624.png"

set xlabel "x"
set ylabel "y"
set size ratio -1
set key off

file = "build/output/nodes_000624.csv"

Ni = 102
Nj = 102

plot \
    for [i=1:Ni] sprintf("< awk -F, 'NR>1 && $2==%d {print $4,$5}' %s", i, file) \
        using 1:2 with lines lw 1, \
    for [j=1:Nj] sprintf("< awk -F, 'NR>1 && $3==%d {print $4,$5}' %s", j, file) \
        using 1:2 with lines lw 1
