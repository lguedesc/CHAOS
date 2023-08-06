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
echo "12 - Full Time Series (angle test)"
echo "13 - Full Bifurcation Diagram (angle test)"
echo "0  - EXIT"
echo "====================================================================="
printf "Plot Module: "
read plt
echo "====================================================================="
if [ $plt -eq 1 ];
then 
    python plots/plot_timeseries.py
elif [ $plt -eq 2 ]
then
    python plots/plot_poinc.py
elif [ $plt -eq 3 ]
then
    python plots/plot_lyap.py
elif [ $plt -eq 4 ]
then
    python plots/plot_ftimeseries.py
elif [ $plt -eq 5 ]
then
    python plots/plot_bifurc.py
elif [ $plt -eq 6 ]
then
    python plots/plot_fbifurc.py
elif [ $plt -eq 7 ]
then
    python plots/plot_dyndiag.py
elif [ $plt -eq 8 ]
then
    python plots/plot_fdyndiag.py
elif [ $plt -eq 9 ]
then
    python plots/plot_epbasin.py
elif [ $plt -eq 10 ]
then
    python plots/plot_forcedbasin.py
elif [ $plt -eq 11 ]
then
    python -B -m StabilityAnalysis.Potential_and_EP-only_mechanical
elif [ $plt -eq 12 ]
then
    python plots/plot_ftimeseries_angles.py
elif [ $plt -eq 13 ]
then
    python plots/plot_fbifurc_angles.py
elif [ $plt -eq 0 ]
then
    exit 0
else
    echo "Invalid Option. Exiting script..."
    exit 0
fi