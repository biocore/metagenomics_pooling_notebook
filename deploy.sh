#!/bin/bash

#########################################################################
# JupyterHub Environment and Kernel Deployment Script
#
# Purpose:
#   Deploy a github repo to a new conda environment and JupyterHub kernel
#   - Implements verification and logging capabilities
#
# Usage:
# bash  ./deploy.sh <repo_name> <tag_name> [kernel_prefix] [--dry-run]
##########################################################################

# ensure error codes from upstream calls are passed through pipes
set -euo pipefail

# Configuration
LOG_FILE="deployment_$(date +%Y%m%d_%H%M%S).log"

# Function to show usage instructions
show_usage() {
  cat << EOF
Usage: $0 <repo_name> <tag_name> [kernel_prefix] [OPTIONS]

Deploy a github repo to a new conda environment and JupyterHub kernel.

Positional arguments:
  repo_name          GitHub repository name (e.g., MyUser/my-repo)
  tag_name           Git tag or branch to deploy (e.g., 2025.05.1+testdeploy)
  kernel_prefix      (Optional) Path to install the Jupyter kernel spec
                     (e.g., /shared/local).  Necessary for installing kernels
                     to a system-wide shared location.

Options:
  --dry-run          Show what would happen without making changes
  --help             Show this help message and exit

Example:
  $0 MyUser/my-repo 2025.05.1+testdeploy
  $0 MyUser/my-repo 2025.05.1+testdeploy /shared/local --dry-run
EOF
  exit 1
}

# Parse command line arguments
parse_args() {
  # Check if at least repo and deploy tag are provided
  if [ $# -lt 2 ]; then
    show_usage
  fi

  GITHUB_REPO=$1
  shift

  DEPLOY_TAG=$1
  shift

    # Check if next arg is the optional kernel prefix (i.e., not a flag like --dry-run)
  if [[ $# -gt 0 && "$1" != --* ]]; then
    KERNEL_PREFIX=$1
    KERNEL_PREFIX=$(echo "$KERNEL_PREFIX" | sed 's:/*$::') # Remove trailing slashes
    shift
  else
    KERNEL_PREFIX=""
  fi

  
  # Default values
  DRY_RUN=false
  
  # Parse optional arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --dry-run)
        DRY_RUN=true
        log "INFO" "Dry run mode enabled - no changes will be made"
        ;;
      --help)
        show_usage
        ;;
      *)
        echo "Unrecognized input: $1"
        exit 1
        ;;
    esac
    shift
  done
}

# Log message with level
log() {
  local level=$1
  local message=$2
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message"
}

# Initialize logging
setup_logging() {
  exec > >(tee -a "$LOG_FILE") 2>&1
  log "INFO" "Starting deployment script..."
}

# Exit with error message
error_exit() {
  log "ERROR" "$1"
  exit 1
}

# Check dependencies
check_dependencies() {
  log "INFO" "Checking dependencies..."
  
  # NB: these are just the dependencies to run this script, not the dependencies for the lab notebooks
  for cmd in conda git python; do
    if ! command -v $cmd &> /dev/null; then
      error_exit "Required command not found: $cmd"
    fi
  done
  
  log "INFO" "All dependencies are available"
}

# Format kernels directory path
format_kernels_dir() {
  local a_kernel_prefix="$1"
  local kernels_dir=""
  if [[ -n "$a_kernel_prefix" ]]; then
    kernels_dir="$a_kernel_prefix/share/jupyter/kernels/"
  fi
  echo "$kernels_dir"
}


# Check if kernel exists
kernel_exists() {
  local kernel_name=$1
  log "INFO" "Checking if kernel '$kernel_name' exists for prefix '$KERNEL_PREFIX'..."

  # Check if the kernels directory exists
  local formatted_kernel_dir=""
  formatted_kernel_dir=$(format_kernels_dir "$KERNEL_PREFIX")
  if [[ -n "$formatted_kernel_dir" ]]; then
    log "INFO" "Checking for existing kernel in $formatted_kernel_dir"
    if [ ! -d "$formatted_kernel_dir" ]; then
        log "ERROR" "Jupyter kernels directory not found at $formatted_kernel_dir"
        return 1
    fi

    # Get all directories in the kernels directory
    for dir in "$formatted_kernel_dir"/*; do
        if [ -d "$dir" ]; then
            # Extract just the kernel name (basename)
            local name=""
            name=$(basename "$dir")

            # Check if it matches the input
            if [ "$name" = "$kernel_name" ]; then
                echo 1 # Kernel exists
                return 0  # Function succeeded
            fi
        fi
    done
  fi

  echo 0 # Kernel does not exist
  return 0  # Function succeeded
}

# Create and set up new environment
setup_new_environment() {
  # Create environment name based on deploy type
  log "INFO" "Setting up new environment '$DEPLOY_NAME'..."

  if [ "$DRY_RUN" = true ]; then
    log "INFO" "DRY RUN: Would create conda environment '$DEPLOY_NAME', then install requirements and repo '$GITHUB_REPO'"
    log "INFO" "DRY_RUN: Would install kernel '$DEPLOY_NAME'"
    return
  fi

  # Clone the repository to get requirements
  TEMP_DIR=$(mktemp -d)
  log "INFO" "Cloning repository to get requirements..."

  # Note that lightweight cloning (e.g. --depth 1) that leaves out full history only works for lightweight (not annotated) tags
  GITHUB_URL="https://github.com/$GITHUB_REPO"
  git clone --depth 1 --branch "$DEPLOY_TAG" "$GITHUB_URL" "$TEMP_DIR" || report "Failed to clone repository"

  # Create new conda environment from environment.yml
  if [ -f "$TEMP_DIR/environment.yml" ]; then
    log "INFO" "Found environment.yml, installing conda environment and dependencies..."
    conda env create --file "$TEMP_DIR/environment.yml" --name "$DEPLOY_NAME"  || report "Failed to install from environment.yml"
  else
    report "Could not find environment.yml"
  fi

  # Install the repo
  log "INFO" "Installing repo $GITHUB_REPO"
  conda run -n "$DEPLOY_NAME" pip install "git+$GITHUB_URL@$DEPLOY_TAG" || rollback "Failed to install repo" "$DEPLOY_NAME"

  # Install the kernel; send to user-specified directory iff KERNEL_PREFIX is set else to new conda env
  if [ -z "$KERNEL_PREFIX" ]; then
    KERNEL_PREFIX=$(conda run -n "$DEPLOY_NAME" python -c 'import sys; print(sys.prefix)')
  fi

  log "INFO" "Installing kernel $DEPLOY_NAME pointing to environment $DEPLOY_NAME in $KERNEL_PREFIX ..."
  conda run -n "$DEPLOY_NAME" python -m ipykernel install --name="$DEPLOY_NAME" --display-name="$DEPLOY_NAME" --prefix="$KERNEL_PREFIX" || rollback "Failed to install kernel" "$DEPLOY_NAME"

  # Clean up
  rm -rf "$TEMP_DIR"
}

# Verify a newly installed environment
verify_environment() {
  local env_name=$1
  local kernel_name=$2
  
  log "INFO" "Verifying environment '$env_name' and kernel '$kernel_name'..."
  
  # Check if environment exists
  if ! conda info --envs | awk '{print $1}' | grep -Fxq "$env_name"; then
    log "ERROR" "Conda environment '$env_name' not found"
    return 1
  fi
  
  # Check if kernel exists
  exists=$(kernel_exists "$kernel_name")
  if [ $? -ne 0 ]; then
    log "ERROR" "Error checking kernel existence"
    return 1
  elif [ "$exists" -eq 0 ]; then
    log "ERROR" "Kernel '$kernel_name' not found"
    return 1
  fi
  
  # Verify kernel can start
  log "INFO" "Testing if kernel can start..."
  
  TEMP_DIR=$(mktemp -d)
  TEMP_NOTEBOOK="$TEMP_DIR/deploy_test.ipynb"
  cat > "$TEMP_NOTEBOOK" << EOF
{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c59cc569-40ad-4881-acde-f4099e79edbf",
   "metadata": {},
   "outputs": [],
   "source": [
    "print('Kernel verification successful')"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "$kernel_name",
   "language": "python",
   "name": "$kernel_name"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
EOF

  log "INFO" "Executing temporary notebook '$TEMP_NOTEBOOK'..."

  # Activate the conda environment and run the notebook
  eval "$(conda shell.bash hook)"
  conda activate "$DEPLOY_NAME"

  local return_value
  if ! jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=60 "$TEMP_NOTEBOOK"; then
    log "ERROR" "Kernel verification failed - kernel could not execute notebook"

    log "WARNING" "Removing kernel '$kernel_name' due to error ..."
    jupyter kernelspec remove -f "$kernel_name" || log "ERROR" "Failed to clean up kernel '$kernel_name'"
    return_value=1
  else
    log "INFO" "Environment and kernel verification successful"
    return_value=0
  fi

  conda deactivate
  rm -r "$TEMP_DIR"
  return $return_value
}

# Report if deployment fails
report() {
  local message=$1
  error_exit "Deployment failed: $message"
}

# Undo creation of environment if downstream steps fail
rollback() {
    local message=$1
    local deploy_name=$2

    log "WARNING" "Removing conda environment '$deploy_name' due to error..."
    conda remove -n "$deploy_name" --all -y
    report $message
}

# Main function
main() {
  parse_args "$@"

  setup_logging
  check_dependencies
  
  log "INFO" "Starting deployment for tag '$DEPLOY_TAG'..."

  # Replace literal periods (.) and plus signs (+) in the tag name with underscores (_)
  DEPLOY_NAME=$(echo "$DEPLOY_TAG" | sed 's/[.+]/_/g')
  
  # Check for existing kernel iff a kernel prefix is provided
  if [[ -n "$KERNEL_PREFIX" ]]; then
    log "INFO" "Checking for pre-existing kernel '$DEPLOY_NAME'..."
    exists=$(kernel_exists "$DEPLOY_NAME")
    if [ $? -ne 0 ]; then
      error_exit "Error checking kernel existence"
    elif [ "$exists" -eq 1 ]; then
      error_exit "Kernel '$DEPLOY_NAME' already exists"
    fi
  fi
  
  # Set up new environment
  setup_new_environment
  
  # Verify the new environment and kernel
  if [ "$DRY_RUN" = false ]; then
    log "INFO" "Verifying new environment..."
    # NB: double use of "$DEPLOY_NAME" is NOT a typo :)
    if ! verify_environment "$DEPLOY_NAME" "$DEPLOY_NAME"; then
      rollback "Environment verification failed" "$DEPLOY_NAME"
    fi
  else
    log "INFO" "DRY RUN: Would verify environment '$DEPLOY_NAME' and kernel '$DEPLOY_NAME'"
  fi
  
  log "INFO" "Deployment successful!"
  log "INFO" "New kernel '$DEPLOY_NAME' is using conda environment '$DEPLOY_NAME'"
  log "INFO" "Log file: $LOG_FILE"
  exit 0
}

# Execute main function with all arguments
main "$@"
