export OMP_NUM_THREADS=8
echo "====================================================================="
echo "I. Choose execution option"
echo "====================================================================="
echo "1 - Compile"
echo "2 - Compile and Run"
echo "3 - Run"
echo "4 - Plot"
echo "0 - EXIT"
echo "====================================================================="
printf "Option: "
read option
echo "====================================================================="
if [ $option -eq 1 ];
then 
    icc -std=c17 -no-multibyte-chars -qopenmp -o CHAOS -O3 src/main.c src/libs/odesystems.c src/libs/interface.c src/libs/iofiles.c src/libs/nldyn.c src/libs/nlosc.c src/libs/customcalc.c src/modules/convergence_test.c src/modules/time_series.c src/modules/poinc_map.c src/modules/lyap_exp_wolf.c src/modules/ftime_series.c src/modules/bifurcation.c src/modules/fbifurcation.c src/modules/dyndiag.c src/modules/fdyndiag.c src/modules/epbasin.c src/modules/forcedbasin.c src/modules/OS_time_series.c src/modules/OS_ftime_series.c src/modules/OS_bifurcation.c src/modules/OS_fbifurcation.c src/modules/OS_dyndiag.c src/modules/OS_fdyndiag.c src/modules/OS_fforcedbasin.c
elif [ $option -eq 2 ]
then 
    icc -std=c17 -no-multibyte-chars -qopenmp -o CHAOS -O3 src/main.c src/libs/odesystems.c src/libs/interface.c src/libs/iofiles.c src/libs/nldyn.c src/libs/nlosc.c src/libs/customcalc.c src/modules/convergence_test.c src/modules/time_series.c src/modules/poinc_map.c src/modules/lyap_exp_wolf.c src/modules/ftime_series.c src/modules/bifurcation.c src/modules/fbifurcation.c src/modules/dyndiag.c src/modules/fdyndiag.c src/modules/epbasin.c src/modules/forcedbasin.c src/modules/OS_time_series.c src/modules/OS_ftime_series.c src/modules/OS_bifurcation.c src/modules/OS_fbifurcation.c src/modules/OS_dyndiag.c src/modules/OS_fdyndiag.c src/modules/OS_fforcedbasin.c
    ./CHAOS
elif [ $option -eq 3 ]
then 
    ./CHAOS
elif [ $option -eq 4 ]
then 
    echo "====================================================================="
    echo "II. Choose the module to plot:"
    echo "====================================================================="
    echo "1 - Time Series"
    echo "2 - Poincaré Map"
    echo "3 - Lyapunov Exponents (Wolf Method)"
    echo "4 - Full Time Series (Integrator + Poincaré + Lyapunov Exponents)"
    echo "5 - Bifurcation Diagram"
    echo "6 - Full Bifurcation Diagram"
    echo "7 - Dynamical Diagram"
    echo "8 - Full Dynamical Diagram"
    echo "9 - Basin of Attraction (Equilibrium Point)"
    echo "10 - Basin of Attraction (Forced)"
    echo "0 - EXIT"
    echo "====================================================================="
    printf "Plot Module: "
    read plt
    echo "====================================================================="
    if [ $plt -eq 1 ];
    then 
        python -B -m TimeSeries.out.plot_rk4
    elif [ $plt -eq 2 ]
    then
        python -B -m PoincareMap.out.plot_poinc
    elif [ $plt -eq 3 ]
    then
        python -B -m LyapunovExp.out.plot_lyap
    elif [ $plt -eq 4 ]
    then
        python -B -m FTimeSeries.out.plot_ftimeseries
    elif [ $plt -eq 5 ]
    then
        python -B -m Bifurcation.out.plot_bifurc
    elif [ $plt -eq 6 ]
    then
        python -B -m FBifurcation.out.plot_fbifurc
    elif [ $plt -eq 7 ]
    then
        python -B -m DynDiagram.out.plot_dyndiag
    elif [ $plt -eq 8 ]
    then
        python -B -m FDynDiagram.out.plot_fdyndiag
    elif [ $plt -eq 9 ]
    then
        python -B -m EPBasin.out.plot_epbasin
    elif [ $plt -eq 10 ]
    then
        python -B -m FForcBasin.out.plot_forcedbasin
    elif [ $plt -eq 0 ]
    then
        exit 0
    else
        echo "Invalid Option. Exiting script..."
        exit 0
    fi
elif [ $option -eq 0 ]
then
    exit 0
else
    echo "Invalid Option. Exiting script..."
    exit 0
fi