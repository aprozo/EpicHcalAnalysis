#!/bin/bash

# Define the input file
input_file="tileMap.txt"

# Check if the input file exists
if [ ! -e "$input_file" ]; then
  echo "Input file '$input_file' does not exist."
  exit 1
fi

# Loop through each line in the input file
while read -r line; do
  # Split the line into two variables, assuming they are separated by a space
  PHI=$(echo "$line" | cut -d' ' -f1)
  THETA=$(echo "$line" | cut -d' ' -f2)

  # Print or use X and Y as needed
  echo "PHI: $PHI, THETA: $THETA"

  # You can perform any actions with X and Y here
done < "$input_file"
