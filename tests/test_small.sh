#!/bin/bash
prog="../bin/tmf"

if [ ! -e "${prog}" ]; then
    echo "Error: Program file '${prog}' does not exist."
    echo "Please run 'make' in directory '../src'."
    exit
fi

# Input data file that contains the events of the temporal network.
data_file='test_data.dat'

# Input data file that contains the node types.
node_types='node_types.dat'

## Output file name (.dat will be appended).
test_output="test_small_output"

# Parameters for detecting temporal motifs.
tw=10         # Time window
motif_size=3  # Max number of events in temporal motifs to find
r=0           # Number of references to create (0 means no references)

# Find motifs up to ${motif_size} events.
${prog} ${tw} ${test_output} -m ${motif_size} -r ${r} -nf ${node_types} < ${data_file}
