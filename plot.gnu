set term pdf enhanced color font "Ubuntu"

set style line 1 linecolor rgb '#e6e6e6' lw 3  ps 0.5
set style line 2 linecolor rgb '#3abbc9' lw 3  ps 0.5
set style line 3 linecolor rgb '#9bca3e' lw 3  ps 0.5
set style line 4 linecolor rgb '#feeb51' lw 3  ps 0.5
set style line 5 linecolor rgb '#ffb92a' lw 3  ps 0.5
set style line 6 linecolor rgb '#ed5314' lw 3  ps 0.5

set size square
set key bottom right
set xtics 0.0,0.04,0.2
set xrange [0:0.2]


set output 'voter_01.pdf'

set xlabel 'Candidate 2 information'
plot 'voter_m.txt' u 1:2:3 w errorbars ls 3 t 'Cand. 2 votes',\
5000 w l ls 6 notitle

set output

set output 'voter_1.pdf'

set xlabel 'Candidate 2 information'
plot 'voter_m2.txt' u 1:2:3 w errorbars ls 3 t 'Cand. 2 votes',\
5000 w l ls 6 notitle

set output