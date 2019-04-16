

./kabs $1 $1 $2 $3 $4
for i in 1 2 3 4 5
do
./kabs $2 $2 $5 $3 $4
./kabs $5 $5 $2 $3 $4
done
compare -compose src $1 $2 dif.tif