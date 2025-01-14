#!/bin/bash
# Specify a partition
#SBATCH --partition=bluemoon
# Request nodes 
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1
# Maximum runtime
#SBATCH --time=0:10:00
# Name of job
#SBATCH --job-name=L_BATCH
# Output of this job, stderr and stdout are joined by default
# %x=job-name %j=jobid
#SBATCH --output=outputs/%x_%j.out

data_files=("${@:1:10}")

source_node1=${11}
source_node2=${12}
output_node=${13}

visualize=${14}
amplitude=${15}

method=${16}
outputs=${17}

cases=("${@: -10}")


for ((i=0; i<${#data_files[@]}; i++)); do
    id=$(echo "${data_files[i]}" | grep -oP 'data_id\K\d+')
    seed=$(echo "${data_files[i]}" | grep -oP 'seed\K\d+')
    period1=0.0666
    #$(echo "${data_files[i]}" | grep -oP 'period\K\d+' | head -n 1)
    period2=0.05
    #$(echo "${data_files[i]}" | grep -oP 'period\K\d+' | tail -n 1)
    IFS=',' read -ra case_values <<< "${cases[i]}"

    for case in "${case_values[@]}"; do
        echo "${period1}"
        echo "${period2}"

        ./lmp_serial -in in.mat_simulator_varied -var data_file ${data_files[i]} \
        -var source_node1 ${source_node1} -var source_node2 ${source_node2} \
        -var output_node ${output_node} -var id ${id} -var visualize ${visualize} \
        -var seed ${seed} -var amplitude ${amplitude} -var period1 ${period1} \
        -var period2 ${period2}  -var method ${method} -var case ${case} \
        -var dump_folder ../dumps/period${period1}_period${period2}/ \
        >../outputs/${outputs}
    done
done