set view map
set xrange [0:1]
set yrange [0:1]
set cbrange [0:1]
set term png

do for [file in system('ls -1B *.data')] {
	plot file with image
	set output file.".png"
	replot
}
