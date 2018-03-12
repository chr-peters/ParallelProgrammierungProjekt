# the name of the .csv file
FILENAME='plot_data.csv'

# get the maximum number of processes to use
MAXPROCS=48
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
echo 'size, rank, type, runtime, calctime, mpitime, waittime' > $FILENAME

# iterate over number of procs
for PROCS in `seq 1 $MAXPROCS`
do
    # iterate over each type
    for TYPE in 0 1 2
    do
	# call the program and append outputs to .csv file
	mpirun -n $PROCS mandelseq -v -t $TYPE -w 4096 -h 4096 -x -.59 -.54 -.58 -.53 -i 1024 >> $FILENAME
    done
done
