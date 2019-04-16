kabsmake: kabsDImaxtree2.c
	gcc -O3 -I../Dist kabsDImaxtree2.c -o kabs -lfreeimage -lm
highpassmake: gaussianblur.c
	gcc -O3 -I../Dist highPass.c -o highPass -lfreeimage -lm
unsharpmake: unsharpmasking.c
	gcc -O3 -I../Dist unsharpmasking.c -o usm -lfreeimage
peakmake: kabsDImaxtree2Peak.c
	gcc -O3 -I ../Dist kabsDImaxtree2Peak.c -o kabsPeak -lfreeimage -lm
all: kabsmake highpassmake unsharpmake peakmake