#!/bin/bash

# Define the number of times you want the loop to run
num_interations=15
# Define the maximum CPU utilization percentage (e.g., 95%)
max_cpu=95

for ((i=1; i<=num_interations; i++))
do
    # Check the CPU utilization
    cpu=$(echo $[100-$(vmstat 1 2|tail -1|awk '{print $15}')])

    # If the CPU utilization is greater than the maximum, wait
    while (( cpu > max_cpu ))
    do
        sleep 5
        cpu=$(echo $[100-$(vmstat 1 2|tail -1|awk '{print $15}')])
    done

    
    # Start the program in a new terminal window
    echo "Running program $i..."
    gnome-terminal -- bash -c "source /opt/intel/oneapi/setvars.sh; echo '2 23 7 data/FBifurcation/in/MHEH_cr_in-$i.txt' | ./bin/CHAOS; exec bash"
done
