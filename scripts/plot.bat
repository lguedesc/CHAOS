@echo off
set pythonpackage=miniconda3
echo =====================================================================
echo Choose the module to plot:
echo =====================================================================
echo 1  - Time Series
echo 2  - Poincare Map
echo 3  - Lyapunov Exponents (Wolf Method)
echo 4  - Full Time Series (Integrator + Poincare Map + Lyapunov Exponents)
echo 5  - Bifurcation Diagram
echo 6  - Full Bifurcation Diagram
echo 7  - Dynamical Diagram
echo 8  - Full Dynamical Diagram
echo 9  - Basin of Attraction (Equilibrium Point)
echo 10 - Basin of Attraction (Forced)
echo 11 - Stability Analysis
echo 0  - EXIT
echo =====================================================================
set /p plt="Plot Module: "

if %plt%==1 goto plot_timeseries
if %plt%==2 goto plot_poincare
if %plt%==3 goto plot_lyapunov
if %plt%==4 goto plot_ftimeseries
if %plt%==5 goto plot_bifurcation
if %plt%==6 goto plot_bifurcation
if %plt%==7 goto plot_dyndiag
if %plt%==8 goto plot_fdyndiag
if %plt%==9 goto plot_epbasin
if %plt%==10 goto plot_fforcbasin
if %plt%==11 goto plot_stability
if %plt%==0 goto exit_script
if %plt% LSS 0 goto invalid
if %plt% GTR 11 goto invalid

:plot_timeseries
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_timeseries.py
    goto exit_script
:plot_poincare
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_poinc.py
    goto exit_script
:plot_lyapunov
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_lyap.py
    goto exit_script
:plot_ftimeseries
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_ftimeseries.py
    goto exit_script
:plot_bifurcation
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_bifurc.py
    goto exit_script
:plot_fbifurcation
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_fbifurc.py
    goto exit_script
:plot_dyndiag
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_dyndiag.py
    goto exit_script
:plot_fdyndiag
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_dyndiag.py
    goto exit_script
:plot_epbasin
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_epbasin.py
    goto exit_script
:plot_fforcbasin
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_forcedbasin.py
    goto exit_script
:plot_stability
    %USERPROFILE%/%pythonpackage%/python.exe plots\plot_stability.py
    goto exit_script
:invalid
    echo Invalid Option. Exiting script...
    goto exit_script
:exit_script
    exit /b

echo =====================================================================