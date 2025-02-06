#!/bin/bash

# Load the necessary module, e.g., Python

#SBATCH --job-name=data_processor_job   # Specify the job name
#SBATCH --output=data_processor_output_%j.txt  # Output file name, where %j is the job ID
#SBATCH --ntasks=1                      # Total number of tasks (change if needed)
#SBATCH --cpus-per-task=1               # Number of CPUs per task


# Execute the Python script
python data_processor.py --start_location 1 --end_position 60 --total_cycle 12 --num_threads 6

