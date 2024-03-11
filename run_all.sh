#!/bin/bash

# # Run DataCollection
# python3 DataCollection/import_fits_images.py "$1" "$2" "$3"
 folder_name="DataCollection/FITSImages_${1}_${2}"
# echo "Images retrieved from CADC placed into $folder_name."
# # Check if the folder exists and is empty
if [ -d "$folder_name" ] && [ -z "$(ls -A "$folder_name")" ]; then
    echo "The folder '$folder_name' is empty. Stopping the script."
    exit 1
fi

# Prepare folder structure for cleaner
python3 n1fits/bin/launch_fits.py
mkdir -p "$n1fits/mission/image/fits_processor/log"
chmod +rw "$n1fits/mission/image/fits_processor/log"
mkdir -p "$n1fits/mission/image/fits_processor/incoming"
chmod +rw "$n1fits/mission/image/fits_processor/incoming"
mkdir -p "$n1fits/mission/image/fits_processor/outgoing"
chmod +rw "$n1fits/mission/image/fits_processor/outgoing"
mkdir -p "$n1fits/mission/image/fits_processor/outgoing/ASTRO"
chmod +rw "$n1fits/mission/image/fits_processor/outgoing/ASTRO"

# Move the images retrieved to the working directory of the cleaner 
working_dir="n1fits/mission/image/fitsprocessor/outgoing/ASTRO/"

# Find all .fits files in the subfolders of the specified folder
# and duplicate them into the working directory's folder
find "$folder_name" -type f -name "*.fits" -exec cp {} "$working_dir" \;
echo "All .fits files have been duplicated into $new_folder."

# Run cleaner
python3 n1fits/bin/launch_fits.py # places cleaned images in **../mission/image/outgoing/ASTRO/**

# Run ObjectDetection
python3 ObjectDetection/detect_objects.py