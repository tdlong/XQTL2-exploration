#!/bin/bash

# Identify and Manage Adaptive Window Slurm Jobs
# This script helps find, cancel, and restart affected jobs

echo "=== Adaptive Window Job Management ==="
echo ""

# Check current jobs
echo "Current slurm jobs for user $USER:"
squeue -u $USER
echo ""

# Look for adaptive window jobs specifically
echo "Looking for adaptive window related jobs..."
squeue -u $USER | grep -E "(REFALT2haps|adaptive|haplotype)" || echo "No obvious adaptive window jobs found"
echo ""

# Check job history for recent adaptive window jobs
echo "Recent job history (last 30 days):"
sacct -u $USER --starttime=$(date -d '30 days ago' +%Y-%m-%d) | grep -E "(REFALT2haps|adaptive|haplotype)" || echo "No recent adaptive window jobs in history"
echo ""

# Find job output files
echo "Looking for slurm output files..."
find . -maxdepth 2 -name "slurm-*.out" -o -name "slurm-*.err" | head -20
echo ""

# Interactive job management
echo "=== Job Management Options ==="
echo "1. Cancel specific job: scancel <JOB_ID>"
echo "2. Cancel all user jobs: scancel -u $USER"
echo "3. Check specific job details: scontrol show job <JOB_ID>"
echo "4. Check job array status: scontrol show job <ARRAY_JOB_ID>"
echo ""

# Check if there are any array jobs
echo "Checking for array jobs..."
squeue -u $USER -o "%.18i %.9P %.20j %.8u %.2t %.10M %.6D %R" | grep -E "\[" || echo "No array jobs found"
echo ""

echo "=== Next Steps ==="
echo "1. Identify the job IDs that need to be cancelled"
echo "2. Cancel them with: scancel <JOB_ID>"
echo "3. Clean up output directories"
echo "4. Restart with: sbatch <your_script>.sh"
echo ""

# Optional: Check specific directories for job outputs
if [ -d "slurm_logs" ]; then
    echo "Found slurm_logs directory:"
    ls -la slurm_logs/ | head -10
fi

if [ -d "jobs" ]; then
    echo "Found jobs directory:"
    ls -la jobs/ | head -10
fi
