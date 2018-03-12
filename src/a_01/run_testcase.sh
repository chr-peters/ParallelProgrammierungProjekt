TYPE=$1

if [ -z $1 ]
then
    TYPE=0
fi
    
mpirun mandelseq -t $TYPE -v -w 800 -h 800 -x -.59 -.54 -.58 -.53 -i 1024
