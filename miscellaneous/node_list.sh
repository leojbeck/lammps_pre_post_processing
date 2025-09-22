#!/bin/bash

# Check if partition name was provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <partition>"
    exit 1
fi

PARTITION="$1"

# ANSI colors
GREEN='\033[1;32m'
RESET='\033[0m'

# Get nodes used by current user's jobs
used_nodes=$(squeue -u "$USER" -h -o "%R" | sort | uniq)

# Header
echo "Partition, Node #, Free CPUs, Used CPUs"

# Get node info
sinfo -p "$PARTITION" -N -o "%n %C" | while read -r node cpuinfo; do
    node_number=$(echo "$node" | sed 's/[^0-9]*//g')  # Extract numeric part
    used_cpus=$(echo "$cpuinfo" | cut -d'/' -f1)     # Used CPUs
    free_cpus=$(echo "$cpuinfo" | cut -d'/' -f2)     # freeable CPUs

    # Check if current node is in the list of used nodes
    if echo "$used_nodes" | grep -q "$node"; then
        echo -e "${GREEN}$PARTITION, $node_number, $free_cpus, $used_cpus${RESET}"
    else
        echo "$PARTITION, $node_number, $free_cpus, $used_cpus"
    fi

done
