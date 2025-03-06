set terminal pngcairo

set output 'error.png'

set title "relativer errors for different methods"


set xlabel "methods"
set ylabel "relative errors"
set logscale y
unset xtics


plot 'data4.txt' using (0):7 with points title 'Primal direct', \
     'data4.txt' using (1):8 with points title 'Dual direct', \
     'data4.txt' using (2):9 with points title 'Primal CG', \
     'data4.txt' using (3):10 with points title 'Primal BDD', \
     'data4.txt' using (4):11 with points title 'Dual TEFI'
