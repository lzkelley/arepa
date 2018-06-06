minnum=0
maxnum=134

for ((i=$minnum; i<=$maxnum; i++))
do
./B-BaseTree param_run.txt $i
done
