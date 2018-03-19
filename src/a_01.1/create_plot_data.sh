# the name of the .csv file
FILENAME='plot_data.csv'

# get the maximum number of processes to use
MAXPROCS=24
if ! [ -z $1 ]
then
    MAXPROCS=$1
fi

# remove existing .csv file
if [ -f $FILENAME ]
then
    rm plot_data.csv
fi

# create .csv file containing the first line
echo 'type,size,runtime' > $FILENAME

# iterate over number of procs
for PROCS in `seq 1 $MAXPROCS`
do
    # iterate over each scheduling strategy
    for STRATEGY in "static" "static,1" "dynamic" "guided"
    do
	echo "PROCS=$PROCS, STRATEGY=$STRATEGY"
	# set the Scheduling Strategy on runtime
	export OMP_SCHEDULE="$STRATEGY"
	export OMP_NUM_THREADS=$PROCS
	# call the program and append outputs to .csv file
	echo "`echo "$STRATEGY" | tr ',' '-'`,`srun mandelseq2 -w 4096 -h 4096 -x -.59 -.54 -.58 -.53 -i 1024`" >> $FILENAME
    done
done
