clear
reset
print "averagie cpu utilization"
set terminal pngcairo transparent enhanced font "arial,25" fontscale 1.0 size 1920, 1080
set key outside bottom center box title "Pipeline Phase(s)" enhanced
set key maxrows 4
set key font ",25" spacing 1 samplen 2.9 width 2 height 1
set xlabel "Time (hours)" font ",25"
set ylabel "Utilization (%)" font ",25"

set output "/mnt/app_hdd/aprabh2/WP_test/test1/sim1M_pairs_ohsu_16t_30s_2014-04-15_09-03-08/post_processed_stats/output_average_cpu_utilization_plot.png"
set title "Average CPU Utilization (%) per Phase\n{/*0.5 sim1M\\_pairs\\_ohsu\\_16t\\_30s\\_2014-04-15\\_09-03-08}" font ",35"
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

plot "/mnt/app_hdd/aprabh2/WP_test/test1/sim1M_pairs_ohsu_16t_30s_2014-04-15_09-03-08/post_processed_stats/2014-04-15_09.41.06_sar.csv" using (t0(timecolumn(1))/3600):2 every ::3 ls 1 t "bwa mem" with lines, \
  '' using ((timecolumn(3)-offset)/3600):4 every ::3 ls 2 t "samtools view" with lines, \
  '' using ((timecolumn(5)-offset)/3600):6 every ::3 ls 3 t "samtools sort" with lines, \
  '' using ((timecolumn(7)-offset)/3600):8 every ::3 ls 4 t "MarkDuplicates" with lines, \
  '' using ((timecolumn(9)-offset)/3600):10 every ::3 ls 5 t "RealignerTargetCreator" with lines, \
  '' using ((timecolumn(11)-offset)/3600):12 every ::3 ls 6 t "IndelRealigner" with lines, \
  '' using ((timecolumn(13)-offset)/3600):14 every ::3 ls 7 t "BaseRecalibrator" with lines, \
  '' using ((timecolumn(15)-offset)/3600):16 every ::3 ls 8 t "PrintReads" with lines