# Command-line args: fin, fout, titin

set terminal png
set output fout

set xrange [0:100]
set yrange [0:system('wc -l < '.fin)+1]
set ylabel 'Couples'
set xlabel 'Percent blocks recovered'
set title titin
set key off

binw = 2
set boxwidth binw
set style fill solid 1.0
bin(x, width) = width * floor(x / width) + width / 2.0

plot fin using (bin($1, binw)):(1.0) smooth freq with boxes
