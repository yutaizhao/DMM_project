set terminal pngcairo

set output 'iter_H_pcg.png'

set title "Primal CG: iter vs H with h/H constant"


set xlabel "H"
set ylabel "iterations"


plot 'data1.txt' using 4:7 with points title 'Primal CG: iter vs H with h/H constant'


set output 'iter_hH_pcg.png'

set title "rimal CG: iter vs h/H with H constant"


set xlabel "h/H"
set ylabel "iterations"
set logscale x


plot 'data2.txt' using 6:7 with points title 'Primal CG: iter vs h/H with H constant'