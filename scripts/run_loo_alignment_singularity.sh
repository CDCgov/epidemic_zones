#!/bin/bash
# Runs loo_analysis.R inside the influenza-clustering singularity container.
# Binds this script's directory into the container (same path on both sides)
# and sets it as the container's working directory, so loo_analysis.R's
# relative paths (loo_input_bundle.RData, loo_checkpoints/, loo_analysis_log.txt,
# loo_analysis_final.RData) all resolve to files here on the host.
#
# Usage:
#   Foreground (quick test):   bash run_loo_singularity.sh
#   Detached, no scheduler:    nohup bash run_loo_singularity.sh &
#   SLURM:                     sbatch run_loo_singularity.sh
#     (edit/uncomment the #SBATCH lines below first)

#SBATCH --job-name=loo_analysis
#SBATCH --output=loo_analysis_%j.out
#SBATCH --error=loo_analysis_%j.err
#SBATCH --time=24:00:00          # <-- adjust to expected runtime
#SBATCH --mem=16G                # <-- adjust to expected memory needs
#SBATCH --cpus-per-task=1
# #SBATCH --partition=YOUR_PARTITION   # <-- uncomment and set if required
# #SBATCH --account=YOUR_ACCOUNT       # <-- uncomment and set if required

set -euo pipefail

SIF="/scicomp/groups-pure/ID-VSDB/GAT/sif_files/influenza-clustering.sif"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$SCRIPT_DIR"

singularity exec \
  --bind "$SCRIPT_DIR:$SCRIPT_DIR" \
  --pwd "$SCRIPT_DIR" \
  "$SIF" \
  Rscript loo_analysis_alignment.R
