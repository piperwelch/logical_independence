#!/bin/bash
# Specify a partition
#SBATCH --partition=bluemoon
# Request nodes 
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1
# Maximum runtime
#SBATCH --time=0:15:00
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

poly=${16}
outputs=${17}

cases=("01" "10" "11")
#sbatch batch10_run.sh 

#{data_files} 1-10
#{self.source1} 11 
#{self.source2} 12
#{self.output}  13
#False 14
#{self.amplitude} 15
#False  16
#{current_datetime} 17

# Loop through each combination of materials, periods, and cases

for data_file in "${data_files[@]}"; do
    id=$(echo "$data_file" | grep -oP 'data_id\K\d+')
    seed=$(echo "$data_file" | grep -oP 'seed\K\d+')
    period1=$(echo "$data_file" | grep -oP 'period\K\d+' | head -n 1)
    period2=$(echo "$data_file" | grep -oP 'period\K\d+' | tail -n 1)
    
    for case in "${cases[@]}"; do # for each case, do the 2 periods 

        ./lmp_serial -in in.mat_simulator -var data_file ${data_file} -var source_node1 ${source_node1} \
        -var source_node2 ${source_node2} -var output_node ${output_node} -var id ${id} -var visualize ${visualize} \
        -var seed ${seed} -var amplitude ${amplitude} -var period1 ${period1} -var period2 ${period2} \
        -var poly ${poly} -var case ${case} -var dump_folder ../dumps/period${period1}_period${period2}/ >../outputs/${outputs}

    done
done