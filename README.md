## NEOSSAT_Image_Processor
### Current steps to run entire object detection pipeline. Will be automated to run in sequence.
Step 1: navigate to DataCollection and run import_fits_images.py. When prompted, you can enter any coordinates. (14, -13) works. It will generate a folder and put the images there.
Step 2: Follow the README in the cleaner if u want to run the cleaner. Not necessary if you are just testing on the output of object detection. 
Step 3: Take the set of images (from the generated folder in DataCollection, or ../mission/image/outgoing/ASTRO/ from the cleaner) and put them in the folder ObjectDetection/mission_images/. Run detect_objects.py.