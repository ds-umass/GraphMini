dir="${1}"
fname="${dir}/snap.txt"

nedge=1000
perorig=80
nba=20
#original graph
bash convert.sh ${1}
#labels
bin/labelize ${fname}
#make daf
bin/outtoDAF ${fname}
#split for continuous/batch
bin/split_tool $fname 80 $nedge $nba
echo bin/split_tool $fname 80 $nedge $nba
#convert start
bash convertNoRelabel.sh ${fname}.start
#make turb 
bin/unlabelled_turb ${fname}.start
#make the updates
for i in $(seq 0 $((${nba}-1)))
do
    echo $nedge > ${fname}.add.$i.upd
    cat ${fname}.add.$i >> ${fname}.add.${i}.upd
    bin/updateconvert ${fname}.add.${i}
done
