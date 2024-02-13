#!/bin/bash

# Define the number of times you want the loop to run
num_interations=30

for ((i=1; i<=num_interations; i++))
do
    # Start the program in the background and get its PID
    echo "2 23 8 data/DynDiagram/in/MHEH_cr_in-$i.txt" | ./bin/CHAOS &
    pid=$!

    # Wait for the program to finish
    wait $pid
done 