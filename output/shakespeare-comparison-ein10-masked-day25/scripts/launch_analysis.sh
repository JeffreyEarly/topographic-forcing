#!/bin/zsh
set -eu

run_directory="/Users/jearly/Documents/OceanKitRepositories/topographic-forcing/output/shakespeare-comparison-ein10-masked-day25"
script_directory="${run_directory}/scripts"
log_path="${run_directory}/logs/analysis.log"
matlab_binary="/Applications/MATLAB_R2025b.app/bin/matlab"

exec "${matlab_binary}" -batch "addpath('${script_directory}'); analyze_masked_validation" > "${log_path}" 2>&1
