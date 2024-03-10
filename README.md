## NEOSSAT_Image_Processor
### Current steps to run entire object detection pipeline. Will be automated to run in sequence.

_Step 1:_ navigate to DataCollection and run import_fits_images.py. When prompted, you can enter any coordinates. (14, -13) works. It will generate a folder and put the images there.

_Step 2:_ Follow the README in the cleaner if u want to run the cleaner. Not necessary if you are just testing on the output of object detection. 

_Step 3:_ Take the set of images (from the generated folder in DataCollection, or ../mission/image/outgoing/ASTRO/ from the cleaner) and put them in the folder ObjectDetection/mission_images/. Run detect_objects.py.
