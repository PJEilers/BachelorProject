default: all
kabsmake: kabsDImaxtree2.c
	gcc -O3 -I../Dist kabsDImaxtree2.c -o kabs -lfreeimage 
highpassmake: highPass.c
	gcc -O3 -I../Dist highPass.c -o highPass -lfreeimage -lm
unsharpmake: unsharpMasking.c
	gcc -O3 -I../Dist unsharpMasking.c -o usm -lfreeimage
all: kabsmake highpassmake unsharpmake