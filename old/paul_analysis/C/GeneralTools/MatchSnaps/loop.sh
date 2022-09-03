minnum=0
maxnum=135

for ((i=$minnum; i<=$maxnum; i++))
do
./MatchSnaps param.txt $i
done
