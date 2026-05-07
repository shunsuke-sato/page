set term postscript eps enhanced color font "Arial, 24"


fs = 0.024189

set xlabel "Time (fs)"
set ylabel "Dipole moment (arb. units)"

unset key
set xrange [0:5]
set output "dipole_moment.eps"
p "td.general/multipoles" u ($2*fs):6 w l lt 1 lc rgb "red" lw 4

unset output