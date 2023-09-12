export OMP_NUM_THREADS=8
DAF=/home/swreinehr/newdaf/daf_parallel_10min
all_pats="size4_0 size4_1 size5_0 size5_1 size6_0 size6_1 size7_0 size7_1 size8_0 size8_1"
#graphs="bio-HS-CX econ-mahindas ia-dnc fb-messages"
graphs="mico cit-Patents com-livejournal"
#"email-Enron com-amazon com-dblp wiki-Vote" #"mico cite" #"wiki-Vote email-Enron com-dblp com-amazon"
#graph_dirs="/home/swreinehr/AutoMine/graph_converter/"
graph_dirs=/home/swreinehr/graphs
for pattern in $all_pats
do
    make bin/l${pattern}_l_best    
    for gd in $graphs; do
	
	graph_dir=${graph_dirs}/${gd}
	output=labelled_${gd}_$(echo $pattern | rev | cut -d'/' -f 1 | rev)
	echo $graph_dir
	# See the documentation on how `graphflow-server` is run,
	# but in essence, it's:
	# graphflow-server [start] [pattern] [additions] [deletions].
        echo automine >> $output
	env time -v timeout 10m bin/l${pattern}_l_best $graph_dir/snap.txt 2>&1 >> $output 
	echo daf >> $output
	env time -v timeout 10m $DAF -d ${graph_dir}/snap.txt.daf -q dafpat/$pattern -n 1 -m 100000000000000 -h 8 2>&1 >> $output
    done
done
