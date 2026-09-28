# Animate mesh files:
#   build/output/nodes_000001.csv
#   ...
#   build/output/nodes_000500.csv
#
# Run with:
#   gnuplot animate_mesh.gp

set terminal gif animate delay 5 loop 0 optimize size 1200,1200
set output "build/output/mesh_animation.gif"

set xlabel "x"
set ylabel "y"
set size ratio -1
set key off

Ni = 102
Nj = 102

# If the mesh extent changes between frames, fixed ranges prevent the
# axes from rescaling and the animation from visually jumping.
# Uncomment and adjust these if you know the global domain limits:
# set xrange [xmin:xmax]
# set yrange [ymin:ymax]

start_frame = 1
end_frame   = 500

# One `plot` command produces one frame in an animated GIF.
do for [frame=start_frame:end_frame] {
    file = sprintf("build/output/nodes_%06d.csv", frame)

    set title sprintf("Mesh - frame %d / %d", frame, end_frame)

    plot \
        for [i=1:Ni] sprintf("< awk -F, 'NR>1 && $2==%d {print $4,$5}' %s", i, file) \
            using 1:2 with lines lw 1, \
        for [j=1:Nj] sprintf("< awk -F, 'NR>1 && $3==%d {print $4,$5}' %s", j, file) \
            using 1:2 with lines lw 1
}

unset output
