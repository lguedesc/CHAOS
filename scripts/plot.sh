echo "====================================================================="
echo "Choose the module to plot:"
echo "====================================================================="
echo "1  - Time Series"
echo "2  - Poincaré Map"
echo "3  - Lyapunov Exponents (Wolf Method)"
echo "4  - Full Time Series (Integrator + Poincaré + Lyapunov Exponents)"
echo "5  - Bifurcation Diagram"
echo "6  - Full Bifurcation Diagram"
echo "7  - Dynamical Diagram"
echo "8  - Full Dynamical Diagram"
echo "9  - Basin of Attraction (Equilibrium Point)"
echo "10 - Basin of Attraction (Forced)"
echo "11 - Stability Analysis"
echo "0  - EXIT"
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
elif [ $plt -eq 11 ]
then
    python -B -m StabilityAnalysis.Potential_and_EP-only_mechanical
elif [ $plt -eq 0 ]
then
    exit 0
else
    echo "Invalid Option. Exiting script..."
    exit 0
fi