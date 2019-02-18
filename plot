set view map
set xrange [0:1]
set yrange [0:1]
set term png

do for [file in system("find . -maxdepth 1 -type f -name '*.data'")] {
	low = system("cat ".file." | sed -e's/  */ /g' | cut -d' ' -f 4 | sort -r | tail -n 1")
	high = system("cat ".file." | sed -e's/  */ /g' | cut -d' ' -f 4 | sort | tail -n 1")

	print file." range ".low."/".high

	set cbrange [low:high]

	plot file with image
	set output file.".png"
	replot
}