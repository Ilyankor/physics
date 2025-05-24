set terminal pngcairo size 3840,2160 font "Helvetica,48"

set xrange [-0.2:1.2]
set yrange [-0.2:0.2]
unset ytics

n = 5
loc = "0.00 0.25 0.50 0.75 1.00"
sizes = "12 36 24 12 24"

N = 300
do for [i=0:N] {
    plot \
        'output.txt' \
        every ::i::i \
        using ($2 + word(loc, 1)):(0):($6 + word(loc, 5)):(0) \
        with vectors \
        nohead \
        linewidth 8 \
        linecolor "black" \
        notitle, \
        for [j=1:n] 'output.txt' \
            every ::i::i \
            using (column(j+1) + word(loc, j)):(0) \
            with points \
            pointtype 7 \
            pointsize word(sizes, j) \
            notitle
}
