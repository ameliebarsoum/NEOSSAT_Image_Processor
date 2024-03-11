## NEOSSAT_Image_Processor
### Current steps to run entire object detection pipeline. Will be automated to run in sequence.

_Step 0:_ When running the code for the first time, make sure you have astroquery installed. You can do this by running `pip install astroquery` on your terminal

_Step 1:_ navigate to DataCollection and run import_fits_images.py. When prompted, you can enter any coordinates. (RA: 290, Dec: -23, radius: 1 (default)) works. It will generate a folder and put the images there.

_Step 2:_ Follow the README in the cleaner if you want to run the cleaner. Not necessary if you are just testing on the output of object detection. 

_Step 3:_ 

Run

```
$ cd ../../../../../
$ cd ObjectDetection
$ mkdir mission_images
```

Take the set of images (from the generated folder in DataCollection, or ../mission/image/outgoing/ASTRO/ from the cleaner) and put them in the folder ObjectDetection/mission_images/. Run detect_objects.py.
