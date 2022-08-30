export OMP_NUM_THREADS=8
echo
echo "====================================================================="
echo "I. Choose execution option:"
echo "====================================================================="
echo "1 - Compile"
echo "2 - Compile and Run"
echo "3 - Run"
echo "4 - Compile, Run and Plot"
echo "5 - Run and Plot"
echo "6 - Plot"
echo "0 - EXIT"
echo "====================================================================="
printf "Option: "
read option
echo "====================================================================="
echo

if [[ $option -eq 1 ]]
then 
    icc -std=c17 -qopenmp -o CHAOS -O3 src/main.c src/libs/edosystems.c src/libs/interface.c src/libs/iofiles.c src/libs/nldyn.c src/modules/time_series.c src/modules/poinc_map.c src/modules/lyap_exp_wolf.c src/modules/ftime_series.c src/modules/bifurcation.c src/modules/fbifurcation.c src/modules/dyndiag.c src/modules/epbasin.c src/modules/forcedbasin.c  
fi