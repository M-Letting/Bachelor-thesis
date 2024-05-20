#!/bin/bash

# Set working directory in Ucloud to / before running (use "cd .."")
# Path to FragPipe, path to workflow FragPipe workflow file, 
# and FragPipe output directory should be set in FP_command_template
# Protein database file (FASTA) is set in workflow file

#Specify the directory where you want .txt files to be stored
text_file_storage="work/Bachelor_project/real_results/exp_text_storage"
mkdir -p "$text_file_storage"

# Specify the directory where you want to start the search
search_directory="work/Folder/exploris_rawfiles"

# Specify the keyword to search for in file names
keyword="HeLa"

# Specify the output file for saving the results
match_output_file="$text_file_storage/exploris_match_locations.txt"

# Use find to search for files containing the keyword in their names
find "$search_directory" -type f -iname "*$keyword*" 2>/dev/null | grep -i "$keyword" | grep -iv "MS1" | grep -iv "mgf" | tee "$match_output_file"

echo "Search complete. Results saved to $match_output_file"

# Following part of script creates manifest files for FragPipe

# Define path to input file containing paths to .raw data files
mani_input_file="$match_output_file"

# Make output directory for manifest files
mani_output_dir="work/Bachelor_project/manifests/exploris_manis"
mkdir -p "$mani_output_dir"

echo "Manifest output directory created $mani_output_dir"

# Counter for file names
count=1

# Read each line from the input file
while IFS= read -r path; do
    # Create a new file for every 10 paths
    if ((count % 10 == 1)); then
        output_file="$mani_output_dir/exp_mani_$((count / 10 + 1)).fp-manifest"
        touch "$output_file"
    fi

    # Extract the base file name
    base_name=$(basename "$path")

    # Append the path to the output file with the specified format
    echo -e "/$path\t$base_name\t\tDDA" >> "$output_file"

    # Increment the counter
    ((count++))
done < "$mani_input_file"

echo "Clusters created in $mani_output_dir"

# Set output directory for .txt file containing paths of all manifest files
manifest_list_output_dir="$text_file_storage/list_of_exploris_manifest.txt"

# Creates a .txt file containing the paths of all manifest files
find "$mani_output_dir" -type f > "$manifest_list_output_dir"

echo "File paths of clusters have been appended to $manifest_list_output_dir"

#All off the above requires specific paths to written where specified
#all of the above works as intended

#Followign part runs FragPipe and saves results in desired folders
# Specifies paths and command template to run FragPipe
FP_command_template=""opt/fragpipe/bin/fragpipe" --headless \
--workflow "work/Bachelor_project/important_stuff/WFF/exploris_WF.workflow" \
--manifest "MANIFEST_FILE" \
--workdir "work/Bachelor_project/real_results/exploris_results/ExpRres_XOXO" \
--ram 191 \
--threads 31"

# Set variable used for naming output directories for FragPipe runs
naming_var=1
# Loop through each path in the list_of_manifest_files.txt file and execute the command
while IFS= read -r path; do
  # Replace "MANIFEST_FILE" with the current path in the list_of_manifest_files.txt
  mkdir -p "work/Bachelor_project/real_results/exploris_results/ExpRres_$naming_var"
  escaped_path=$(echo "$path" | sed 's/\//\\\//g')
  escaped_naming_var=$(echo "$naming_var" | sed 's/\//\\\//g')
  full_command=$(echo "$FP_command_template" | sed "s/MANIFEST_FILE/$escaped_path/g" | sed "s/XOXO/$escaped_naming_var/g")
  echo "Executing: $full_command"
  eval "$full_command"
  ((naming_var++))
done < "$manifest_list_output_dir"

echo "Will now go through missed raw files and analyse them if possible"

# Set the directory path to results folder
directory_path="work/Bachelor_project/real_results/exploris_results"

# Set the output file for list of missed folders 
missed_folders_output_file="$text_file_storage/folders_under_1MiB.txt"

# Find immediate folders under 1 MiB and write their paths to a text file
find "$directory_path" -maxdepth 1 -type d -exec du -s {} + | awk '$1 < 1024' | cut -f2- > "$missed_folders_output_file"

echo "Immediate folders under 1 MiB have been listed in $missed_folders_output_file"

# Set the output file for subfolder names
subfolder_file="$text_file_storage/subfolder_names.txt"

# Extract subfolder names, replace the last underscore with dot, and write them to a new text file
while IFS= read -r folder_path; do
    find "$folder_path" -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | sed 's/_[^_]*$/./' >> "$subfolder_file"
done < "$missed_folders_output_file"

echo "Subfolder names (with replaced underscore) have been listed in $subfolder_file"

#Set output file for matching file paths
matching_files="$text_file_storage/matching_missed_files.txt"

#Search the specified folder for the filenames from subfolder_file
while IFS= read -r filename; do
    matching_path=$(find "work/Folder/exploris_rawfiles" -type f -name "$filename")
    if [ -n "$matching_path" ]; then
        echo "$matching_path" >> "$matching_files"
    fi
done < "$subfolder_file"

echo "Matching file paths have been listed in $matching_files"

#Make new manifest files from paths of missed files
# Counter for file names
count=1

# Read each line from the input file
while IFS= read -r path; do
    output_file="$mani_output_dir/exp_mani_missed_$count.fp-manifest"
    touch "$output_file"

    # Extract the base file name
    base_name=$(basename "$path")

    # Append the path to the output file with the specified format
    echo -e "/$path\t$base_name\t\tDDA" >> "$output_file"

    # Increment the counter
    ((count++))
done < "$matching_files"

echo "Manifests created in $mani_output_dir"

# Set output directory for .txt file containing paths of all missed manifest files
missed_manifest_list_output_dir="$text_file_storage/list_of_missed_exp_manifest.txt"

# Creates a .txt file containing the paths of all manifest files
find "$mani_output_dir" -type f -name "*missed*" > "$missed_manifest_list_output_dir"

echo "File paths of missed clusters have been appended to $missed_manifest_list_output_dir"

# Set variable used for naming output directories for FragPipe runs
naming_var=1
# Loop through each path in the list_of_manifest_files.txt file and execute the command
while IFS= read -r path; do
  # Replace "MANIFEST_FILE" with the current path in the list_of_manifest_files.txt
  mkdir -p "work/Bachelor_project/real_results/exploris_results/ExpRes_miss_$naming_var"
  escaped_path=$(echo "$path" | sed 's/\//\\\//g')
  escaped_naming_var=$(echo "$naming_var" | sed 's/\//\\\//g')
  full_command=$(echo "$FP_command_template" | sed "s/MANIFEST_FILE/$escaped_path/g" | sed "s/XOXO/miss_$escaped_naming_var/g")
  echo "Executing: $full_command"
  eval "$full_command"
  ((naming_var++))
done < "$missed_manifest_list_output_dir"
