 
Program for adaptive unsharp masking using max-tree filtering, need to have freeimage.h with the "Dist" in the parent folder.

Usage:

make
./highPass <image> <ksize> [beta]
./kabs <original inputimage> <processed input image> <outputfile> <lambda> [k]
./usm <original inputimage> <beta>

Example:

make
./highPass image.tif 11
./kabs pos.tif pos.tif posf.tif 7500 2000
./kabs neg.tif neg.tif negf.tif 7500 2000
./usm image.tif posf.tif negf.tif 25

Also included is a script that runs the max-tree filtering a number of times to get closer to idempotency.

Example:

make
./highPass image.tif 11
./run.sh pos.tif posf.tif 7500 2000 x.tif
./run.sh neg.tif negf.tif 7500 2000 x.tif
./usm image.tif posf.tif negf.tif 25