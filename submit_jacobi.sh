#!/bin/bash
mkdir run
cd run

cp ../bin_jacobi .
cp ../plot_out.py .
# run script

./bin_jacobi
python3 plot_out.py
