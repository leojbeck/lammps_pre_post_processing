#!/bin/bash

# ── argument parsing ───────────────────────────────────────────────────────────

SUMMARY=0
SUMMARY_N=0
PARTITION=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            echo "Usage: $0 <partition> [-s <N>]"
            echo ""
            echo "  <partition>       Slurm partition to query"
            echo "  -s, --summary N   Print aggregate counts instead of per-node list;"
            echo "                    N = minimum free CPUs threshold"
            echo ""
            echo "Nodes running your jobs are highlighted in green."
            exit 0
            ;;
        -s|--summary)
            SUMMARY=1
            SUMMARY_N="$2"
            shift 2
            ;;
        *)
            PARTITION="$1"
            shift
            ;;
    esac
done

# Check if partition name was provided
if [[ -z "$PARTITION" ]]; then
    echo "Usage: $0 <partition> [-s|--summary <N>]"
    exit 1
fi

# Check if summary flag is set and if the argument is a valid integer
if [[ "$SUMMARY" -eq 1 ]] && ! [[ "$SUMMARY_N" =~ ^[0-9]+$ ]]; then
    echo "Error: --summary requires an integer argument"
    exit 1
fi

# ── setup ──────────────────────────────────────────────────────────────────────

# ANSI colors
GREEN='\033[1;32m'
RESET='\033[0m'

# Get nodes used by current user's jobs
used_nodes=$(squeue -u "$USER" -h -o "%R" | sort | uniq)

# ── output ─────────────────────────────────────────────────────────────────────

if [[ "$SUMMARY" -eq 0 ]]; then
    # Detail view: one line per node
    echo "Partition, Node #, Free CPUs, Used CPUs"
    while read -r node cpuinfo; do
        node_number=$(echo "$node" | sed 's/[^0-9]*//g')
        used_cpus=$(echo "$cpuinfo" | cut -d'/' -f1)
        free_cpus=$(echo "$cpuinfo" | cut -d'/' -f2)
        # Highlight nodes the current user has jobs on
        if echo "$used_nodes" | grep -q "$node"; then
            echo -e "${GREEN}$PARTITION, $node_number, $free_cpus, $used_cpus${RESET}"
        else
            echo "$PARTITION, $node_number, $free_cpus, $used_cpus"
        fi
    done < <(sinfo -p "$PARTITION" -N -h -o "%n %C")
else
    # Summary view: aggregate counts
    total=0
    any_free=0
    at_least_n=0
    while read -r node cpuinfo; do
        free_cpus=$(echo "$cpuinfo" | cut -d'/' -f2)
        (( total++ ))
        (( free_cpus > 0          )) && (( any_free++ ))
        (( free_cpus >= SUMMARY_N )) && (( at_least_n++ ))
    done < <(sinfo -p "$PARTITION" -N -h -o "%n %C")

    echo "Partition:                    $PARTITION"
    echo "Total nodes:                  $total"
    echo "Nodes with any free CPUs:     $any_free"
    echo "Nodes with >= $SUMMARY_N free CPUs:  $at_least_n"
fi
