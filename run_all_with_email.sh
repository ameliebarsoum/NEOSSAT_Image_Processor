#!/bin/bash

# TODO: Activate the virtual environment
# source  PATH/TO/VENV/activate

cd /Users/rebeccamizrahi/Documents/McGill/Capstone/NEOSSAT_Image_Processor

# Remove intermediary folders if they exist. Will do nothing if not.
    # From cleaner
rm -rf "n1fits/mission/image/fits_processor"
rm -rf n1fits/mission/image/outgoing/ASTRO/*
rm -rf "n1fits/mission/image/fits_processor/incoming"
rm -rf "n1fits/mission/image/fits_processor/outgoing"
rm -rf "n1fits/mission/image/fits_processor/outgoing/ASTRO"

    # From ObjectDetection
rm -rf "ObjectDetection/aligned/"
rm -rf "ObjectDetection/differenced/"
rm -rf "ObjectDetection/stacked/"
rm -rf "ObjectDetection/flagged/"
rm -rf "ObjectDetection/cosmic_ray_hits/"
rm -rf "ObjectDetection/flagged_original/"
rm -rf "ObjectDetection/flagged_original_fits/"
rm -rf "ObjectDetection/filtered_flagged_original/"

find . -name ".DS_Store" -delete

cd DataCollection
paths=$(python3 import_fits_automated.py)
cd ..

# Append DataCollection/ to the beginning of each path in paths
paths=$(echo "$paths" | sed "s/^/DataCollection\//")
# Echo paths
echo "$paths"

# Iterate through each path of FITS image sets
for path in $paths; do

    # Check if the path contains the expected folder name
    if [[ $path != *"FITSImages_"* ]]; then
        continue
    fi

    echo "Processing $path..."

    # Confirm the folder exists and is not empty
    if [ -d "$path" ] && [ -z "$(ls -A "$path")" ]; then
        echo "The folder '$path' is empty. Stopping the script."
        exit 1
    fi  

    # Prepare folder structure for cleaner
    echo "Preparing folder structure for cleaner..."
    python3 n1fits/bin/launch_fits.py
    mkdir -p "n1fits/mission/image/fits_processor/log"
    chmod +rwx "n1fits/mission/image/fits_processor/log"
    mkdir -p "n1fits/mission/image/fits_processor/incoming"
    chmod +rwx "n1fits/mission/image/fits_processor/incoming"
    mkdir -p "n1fits/mission/image/fits_processor/outgoing"
    chmod +rwx "n1fits/mission/image/fits_processor/outgoing"
    mkdir -p "n1fits/mission/image/fits_processor/outgoing/ASTRO"
    chmod +rwx "n1fits/mission/image/fits_processor/outgoing/ASTRO"
    
    # Move the images retrieved to the working directory of the cleaner 
    working_dir="n1fits/mission/image/fits_processor/outgoing/ASTRO/"

    # Find all .fits files in the subfolders of the specified folder
    # and duplicate them into the working directory's folder
    find "$path" -type f -name "*.fits" -exec cp {} "$working_dir" \;
    echo "All .fits files have been duplicated into $new_folder. Launching cleaner..."

    # Run cleaner
    python3 n1fits/bin/launch_fits.py # places cleaned images in **../mission/image/outgoing/ASTRO/**

    # Run ObjectDetection
    cd ObjectDetection
    python3 detect_objects.py
    cd ..

    # Run ObjectClassification from ObjectClassification/classify_object.py
    cd ObjectClassification
    python3 classify_object.py
    cd ..

    python3 upload_folders_to_drive.py "$path"

    # Run email trigger 
    cd Notification
    python3 notify.py
    cd ..

    # Remove intermediary folders
        # From cleaner
    rm -rf "n1fits/mission/image/fits_processor"
    rm -rf n1fits/mission/image/outgoing/ASTRO/*
        # From ObjectDetection
    rm -rf "ObjectDetection/aligned/"
    rm -rf "ObjectDetection/differenced/"
    rm -rf "ObjectDetection/stacked/"
    rm -rf "ObjectDetection/flagged/"
    rm -rf "ObjectDetection/cosmic_ray_hits/"
    rm -rf "ObjectDetection/flagged_original/"
    rm -rf "ObjectDetection/flagged_original_fits/"
    rm -rf "ObjectDetection/filtered_flagged_original/"
    
    rm -rf "$path"

    echo "* * * *"
    echo " "

done

echo " "
echo " . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ."
echo " "

# Deactivate the virtual environment
# deactivate