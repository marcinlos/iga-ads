#!/usr/bin/env gnuplot
# "USAGE: gnuplot -e "start=0; end=300; step=10" make_gif.gnuplot"

# default values if not provided from command line
if (!exists("start"))   start = 0
if (!exists("end"))     end = 9900
if (!exists("step"))    step = 10

set terminal gif animate delay 20 optimize size 800,600
set output 'anim.gif'

unset key
# unset colorbox
set cbrange [200:*] # sets lower bound of colorbox to 200 while allowing for dynamical scaling of the upper bound


do for [i=start:end:step] {
    filename = sprintf("out_%d.data", i)
    if (system(sprintf("ls %s 2>/dev/null", filename)) eq "") {
        print sprintf("Skipping missing file: %s", filename)
        continue
    }
    print sprintf("Plotting %s", filename)
    plot filename with image
}

unset output