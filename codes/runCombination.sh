#!/bin/bash

# Set Rscript binary file here
R_BIN=/home/perez053/R/R3.4.3/bin/Rscript

# Set total available cores here
totalCores=70

# Set the F combination that you want, e.g. (1 2 3) or
# for sequence from 1 to 15: combination=($(seq 1 15))
#combination=($(seq 71 119))
combination=(41 42 43 51)

# Set the number of iterations (e.g. ($(seq 1 100)))
#iteration=(1 2 3 4)
iteration=($(seq 1 100))


## DO NOT EDIT ANYTHING BELOW ##

# Run combination in parallel, but only for number of cores
count=0
for combIndex in "${combination[@]}"
do
  for iterIndex in "${iteration[@]}"
  do
    $R_BIN --vanilla -e "combIndex <- $combIndex" -e "iterIndex <- $iterIndex" -e "source(\"algorithm.R\")" &
    (( count ++ ))
    if [[ $count -eq $totalCores  ]]; then
        wait
        count=0
    fi
  done
done

wait

echo All done!
