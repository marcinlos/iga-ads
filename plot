
unset key

set xrange [0:1]
set yrange [0:1]
set cbrange [0:value_range]

# set terminal pngcairo size 500, 450
set term png
set output output_file

plot input_file with image
