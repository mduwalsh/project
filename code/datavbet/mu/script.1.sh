#PBS -l nodes=1
#PBS -N vbet.b14g100
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR
./vbet 14 1474447344 100 1
