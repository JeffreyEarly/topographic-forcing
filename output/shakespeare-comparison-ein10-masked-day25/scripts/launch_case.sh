#!/bin/zsh
set -eu

case_name="$1"
run_directory="/Users/jearly/Documents/OceanKitRepositories/topographic-forcing/output/shakespeare-comparison-ein10-masked-day25"
script_directory="${run_directory}/scripts"
log_path="${run_directory}/logs/${case_name}.log"
matlab_binary="/Applications/MATLAB_R2025b.app/bin/matlab"

exec "${matlab_binary}" -batch "addpath('${script_directory}'); run_masked_case('${case_name}')" > "${log_path}" 2>&1
