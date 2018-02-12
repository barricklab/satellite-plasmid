#global analysis - hard to parse
#nucmer --maxmatch -g 5 -l 8 -c 8 --forward --nosimplify --prefix=output pBTK555.fasta pBTK555.fasta
#delta-filter -i 8 output.delta > output.fdelta
#show-coords -r output.fdelta > output.coords

# Note that you actually have to make -g negative to avoid clustering overlapping exact matches!!
NUCMER_OPTIONS="--maxmatch -f -g -100 --noextend -d 0 -l 7 -c 7 --forward --nosimplify"

nucmer $NUCMER_OPTIONS --prefix=satellite satellite-left.fa satellite-right.fa
delta-filter -i 8 satellite.delta > satellite.fdelta
show-coords -T -r satellite.fdelta > satellite.coords.tab

nucmer $NUCMER_OPTIONS --prefix=gfp-del-1 gfp-del-1-left.fa gfp-del-1-right.fa
delta-filter -i 8 gfp-del-1.delta > gfp-del-1.fdelta
show-coords -T -r gfp-del-1.fdelta > gfp-del-1.coords.tab

nucmer $NUCMER_OPTIONS --prefix=gfp-del-2 gfp-del-2-left.fa gfp-del-2-right.fa
delta-filter -i 8 gfp-del-2.delta > gfp-del-2.fdelta
show-coords -T -r gfp-del-2.fdelta > gfp-del-2.coords.tab
