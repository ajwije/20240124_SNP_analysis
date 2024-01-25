#!/bin/bash

# Check if no arguments are passed
if [ $# -eq 0 ]; then
    echo "No directories supplied. Please provide the directories as arguments."
    exit 1
fi

# Loop through all arguments
for dir in "$@"; do
    # Check if directory does not exist
    if [ ! -d "$dir" ]; then
        # Create directory
        mkdir -p "$dir"
        echo "Directory created: $dir"
    else
        echo "Directory already exists: $dir"
    fi
done

