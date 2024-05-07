#!/bin/bash

# Help
Help()
{
   # Display Help
   echo "This program merges files with the same structure, mantaining the header of the first file. If the merge is successful, the input files are removed."
   echo
   echo "Syntax: merge_files.sh file1.txt file2.txt ... filen.txt output.txt"
   echo
}

# Store arguments in array
Help
input=("$@")
output=${@: -1}

unset "input[${#input[@]}-1]"

# Check if output exists
if test -f "$output"
then
    echo "Output file $output exist, specify a new file"
    exit 1
fi

# Merge files
echo "Files to merge:     " ${input[@]}
echo "Creating output in: " $output
echo

awk '
    NR == 1 {print}
    FNR == 1 {next}
    {print}
' ${input[@]} > $output

# Check number of lines and remove original files
lines_input=`wc -l ${input[@]} | tail -n 1 | grep -Eo '[0-9]+'`
lines_output=`wc -l $output | grep -Eo '^[0-9]+'`

if (( $lines_input - $# == $lines_output - 2 ))
then
    echo "Merging successful, removing input files"
    rm ${input[@]}
else
    echo "Merging not successful"
fi
