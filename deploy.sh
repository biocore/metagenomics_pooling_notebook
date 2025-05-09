#!/bin/bash

#############################################################
# JupyterHub Environment and Kernel Deployment Script
#
# Purpose:
#   Deploy a github repo to a JupyterHub kernel environment
#   - Implements verification and logging capabilities
#
# Usage:
# bash  ./deploy.sh <repo_name> <tag_name> [--dry-run]
#############################################################

# ensure error codes from upstream calls are passed through pipes
set -euo pipefail

# Configuration
LOG_FILE="deployment_$(date +%Y%m%d_%H%M%S).log"

# Initialize logging
setup_logging() {
  exec > >(tee -a "$LOG_FILE") 2>&1
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] Starting deployment script"
}

# Log message with level
log() {
  local level=$1
  local message=$2
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message"
}

# Function to show usage instructions
show_usage() {
  cat << EOF
Usage: $0 <repo_name> <tag_name> [OPTIONS]

Deploy a github repo to a JupyterHub kernel environment.

Options:
  --dry-run          Show what would happen without making changes
  --help             Show this help message and exit

Example:
  $0 MyUser/my-repo 2025.05.1+testdeploy
EOF
  exit 1
}

# Exit with error message
error_exit() {
  log "ERROR" "$1"
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
        error_exit "Unknown option: $1"
        ;;
    esac
    shift
  done
}

# Check dependencies
check_dependencies() {
  log "INFO" "Checking dependencies..."
  
  # NB: these are just the dependencies to run this script, not the dependencies for the lab notebooks
  for cmd in conda jupyter jq curl git python; do
    if ! command -v $cmd &> /dev/null; then
      error_exit "Required command not found: $cmd"
    fi
  done
  
  log "INFO" "All dependencies are available"
}


# Check if kernel exists
kernel_exists() {
  local kernel_name=$1
  log "INFO" "Checking if kernel '$kernel_name' exists"
  jupyter kernelspec list --json | jq -e --arg name "$kernel_name" '.kernelspecs[$name]' > /dev/null
  return $?
}


# Verify a newly installed environment
verify_environment() {
  local env_name=$1
  local kernel_name=$2
  
  log "INFO" "Verifying environment '$env_name' and kernel '$kernel_name'"

  # TODO: remove fake fail inserted for testing
  return 1
  
  # Check if environment exists
  if ! conda info --envs | awk '{print $1}' | grep -Fxq "$env_name"; then
    log "ERROR" "Conda environment '$env_name' not found"
    return 1
  fi
  
  # Check if kernel exists
  if ! kernel_exists "$kernel_name"; then
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

  log "INFO" "Execute temporary notebook '$TEMP_NOTEBOOK'"
  
  # Execute the notebook
  if ! jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=60 "$TEMP_NOTEBOOK"; then
    log "ERROR" "Kernel verification failed - kernel could not execute notebook"

    log "INFO" "Removing kernel '$DEPLOY_NAME'"
    jupyter kernelspec remove -f "$DEPLOY_NAME"
    rm -r "$TEMP_DIR"
    return 1
  fi
   
  rm -r "$TEMP_DIR"
  log "INFO" "Environment and kernel verification successful"
  return 0
}

# Report if deployment fails
report() {
  local message=$1
  
  log "ERROR" "Deployment failed: $message"
  exit 1
}

# Undo creation of environment if downstream steps fail
rollback() {
    local message=$1
    local deploy_name=$2

    log "WARNING" "removing conda environment $deploy_name"
    conda remove -n "$deploy_name" --all -y
    report $message
}


# Create and set up new environment
setup_new_environment() {
  
  # Create environment name based on deploy type
  log "INFO" "Setting up new environment: $DEPLOY_NAME"
  
  if [ "$DRY_RUN" = true ]; then
    log "INFO" "DRY RUN: Would create conda environment $DEPLOY_NAME and install requirements and repo $GITHUB_REPO"
    log "INFO" "DRY RUN: Would create kernel $DEPLOY_NAME"
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
    log "INFO" "Found environment.yml, installing conda environment and dependencies ..."
    conda env create --file "$TEMP_DIR/environment.yml" --name "$DEPLOY_NAME"  || report "Failed to install from environment.yml"
  else
    report "Could not find environment.yml"
  fi

  # Install the repo
  log "INFO" "Installing repo $GITHUB_REPO"
  conda run -n "$DEPLOY_NAME" pip install "git+$GITHUB_URL@$DEPLOY_TAG" || rollback "Failed to install repo" "$DEPLOY_NAME"

  # Install the kernel
  log "INFO" "Installing kernel $DEPLOY_NAME pointing to environment $DEPLOY_NAME..."
  conda run -n "$DEPLOY_NAME" python -m ipykernel install --user --name="$DEPLOY_NAME" --display-name="$DEPLOY_NAME" || rollback "Failed to install kernel" "$DEPLOY_NAME"
  
  # Clean up
  rm -rf "$TEMP_DIR"
}

# Main function
main() {
  parse_args "$@"

  setup_logging
  check_dependencies
  
  log "INFO" "Starting deployment for tag '$DEPLOY_TAG'"
  
  DEPLOY_NAME=$(echo "$DEPLOY_TAG" | sed 's/[.+]/_/g')
  
  # Check for existing kernel
  if kernel_exists $DEPLOY_NAME; then
    log "ERROR" "Kernel '$DEPLOY_NAME' already exists"
    exit 1
  fi
  
  # Setup new environment
  setup_new_environment
  
  # Verify the new environment and kernel
  if [ "$DRY_RUN" = false ]; then
    log "INFO" "Verifying new environment..."
    # NB: double use of "$DEPLOY_NAME" is NOT a typo :)
    if ! verify_environment "$DEPLOY_NAME" "$DEPLOY_NAME"; then
      rollback "Environment verification failed" "$DEPLOY_NAME"
    fi
  else
    log "INFO" "DRY RUN: Would verify environment $DEPLOY_NAME and kernel $DEPLOY_NAME"
  fi
  
  log "INFO" "Deployment successful!"
  log "INFO" "New kernel '$DEPLOY_NAME' is using conda environment '$DEPLOY_NAME'"
  log "INFO" "Log file: $LOG_FILE"
  exit 0
}

# Execute main function with all arguments
main "$@"
