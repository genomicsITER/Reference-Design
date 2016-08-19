#################################################################
# Copyright (c) Intel Corporation
# PLEASE DO NOT SHARE WITHOUT NDA
#################################################################


clear
reset
print "averagie cpu utilization"
set terminal pngcairo transparent enhanced font "arial,25" fontscale 1.0 size 1920, 1080
set key outside bottom center box title "Workflow Phase(s)" enhanced
set key maxrows 4
set key font ",25" spacing 1 samplen 2.9 width 2 height 1
set xlabel "Time (hours)" font ",25"
set ylabel "Utilization (%)" font ",25"

set output "/home/jpinsonault/workflow_profiler_testing/workflow_stats_parser/_test_multi_stage/output_average_cpu_utilization_plot.png"
set title "Average CPU Utilization (%) per Phase\n{/*0.5 sample\\_multistage\\_input}" font ",35"
set datafile separator ","
#set xdata time
set timefmt "%Y-%m-%d %H:%M:%S"
#set xtics format "%d:%H:%M" font ",25"
set ytics font ",25"

set style line 1 lt 1 lc rgb "red" lw 4
set style line 2 lt 1 lc rgb "orange" lw 4
set style line 3 lt 1 lc rgb "brown" lw 4
set style line 4 lt 1 lc rgb "green" lw 4
set style line 5 lt 1 lc rgb "cyan" lw 4
set style line 6 lt 1 lc rgb "blue" lw 4
set style line 7 lt 1 lc rgb "violet" lw 4
set style line 8 lt 1 lc rgb "yellow" lw 4
set style line 9 lt 1 lc rgb "green" lw 4
set style line 10 lt 1 lc rgb "cyan" lw 4
set style line 11 lt 1 lc rgb "blue" lw 4
set style line 12 lt 1 lc rgb "violet" lw 4
show style line

offset = 0
starting_time = 37824
t0(x)=(offset=($0==0) ? x : offset, x - offset)

plot "/home/jpinsonault/workflow_profiler_testing/workflow_stats_parser/_test_multi_stage/2014-05-02_15.27.29_sar.csv" using (t0(timecolumn(1))/3600):2 every ::3 ls 1 t "Stage1" with lines, \
  '' using ((timecolumn(3)-offset)/3600):4 every ::3 ls 2 t "Stage2" with lines, \
  '' using ((timecolumn(5)-offset)/3600):6 every ::3 ls 3 t "Stage3" with lines