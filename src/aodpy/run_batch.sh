#!/bin/bash

# Path to the config file
CONFIG_FILE="config.toml"

# Path to the permutations file
BATCH_FILE="batch.txt"

# Function to update the config file
update_config() {
    local site=$1
    local startdate=$2
    local enddate=$3

    # Use sed to update the config file
    sed -i "s/^site = .*/site = \"$site\"/" "$CONFIG_FILE"
    sed -i "s/^startdate = .*/startdate = $startdate/" "$CONFIG_FILE"
    sed -i "s/^enddate = .*/enddate = $enddate/" "$CONFIG_FILE"
}

# Read each line from the permutations file
while IFS=, read -r site startdate enddate; do
    # Update the config file with the current permutation
    update_config "$site" "$startdate" "$enddate"

    #python envcal.py -b || { exit 1; }
    python langley.py || { exit 1; }
    sed -i "s/^calfile_in = .*/calfile_in = \"lcl\"/" "$CONFIG_FILE"
    sed -i "s/^refchannel = .*/refchannel = 870/" "$CONFIG_FILE"
    python gencal.py || { exit 1; }
    sed -i "s/^calfile_in = .*/calfile_in = \"870\"/" "$CONFIG_FILE"
    sed -i "s/^refchannel = .*/refchannel = 670/" "$CONFIG_FILE"
    python gencal.py || { exit 1; }
    sed -i "s/^calfile_in = .*/calfile_in = \"670\"/" "$CONFIG_FILE"
    sed -i "s/^refchannel = .*/refchannel = 500/" "$CONFIG_FILE"
    python gencal.py || { exit 1; }
    python aod2p.py  || { exit 1; }

done < "$BATCH_FILE"
