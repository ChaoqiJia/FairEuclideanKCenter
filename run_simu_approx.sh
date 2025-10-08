#!/bin/bash

# ----------------------------
# Basic configuration
# ----------------------------
PYTHON_SCRIPT="simulate_approx.py"  
ITER_NUM=10                         
OUTPUT_DIR="output"
LOG_FILE="${OUTPUT_DIR}/simu_approx_$(date +%Y%m%d_%H%M%S).log"


# virtual environment
# source ~/venvs/fkc_env/bin/activate
# echo "Activated virtual environment."

mkdir -p "${OUTPUT_DIR}"

# ----------------------------
# Run experiment
# ----------------------------
echo "Starting simulation with ${ITER_NUM} iterations..."
echo "Logs will be saved to ${LOG_FILE}"

python3 "${PYTHON_SCRIPT}" "${ITER_NUM}" > "${LOG_FILE}" 2>&1

# ----------------------------
# Post-run summary
# ----------------------------
if [ $? -eq 0 ]; then
    echo "completed successfully!"
else
    echo " an error. Check the log: ${LOG_FILE}"
fi
