set terminal pngcairo

set output 'time_H.png'

set title "time vs H for different methods"


set xlabel "H"
set ylabel "time"
set logscale x


plot 'data3.txt' using 4:8 with points title 'Dual direct', \
     'data3.txt' using 4:9 with points title 'Primal CG', \
     'data3.txt' using 4:10 with points title 'Primal BDD', \
     'data3.txt' using 4:11 with points title 'Dual TEFI'

