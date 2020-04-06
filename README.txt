Installation:
For Chrome users, download the dendritic_segmentation_Integration.mlapp file.
For other browsers, download dendritic_segmentation_Integration_exported.m
Log into MATLAB Online on a Chrome browser and upload the file to your MATLAB online repository.
Upload the images and their associated Selection.csv files to your MATLAB online repository.

Usage:
Open the folder in your MATLAB repository in which you wish to save the generated MAT files to.
Run the app be right clicking the dendritic_segmentation_Integration.mlapp or dendritic_segmentation_Integration_exported.m file and clicking run.
Click the "Select TIF" button and choose the desired TIF to be analyzed.
Click the "Select CSV" button and choose the CSV corresponding to the TIF.
Click the "Display TIF Channels" button to see each channel in the TIF.
Set the "Nuclei Channel" and "DC Channel" spinners to the corresponding channels as seen in the TIF channels.
(Blue should be the Nuclei Channel, Green should be the DC Channel, Red should be Epithelium the Channel)
Click "Run" and note the figure generated is an overlay of the DC Channel (Green), Nuclei Channel (Blue), and located DC nuclei (Red).
Ensure that a MAT file has been generated in your "Current Folder".
Checking "debug" will generate 10 figures of various parts of the image analysis when the app is run again.
If necessary adjust the "Threshold" slider to more or less severe thresholding to get the correct nuclei.
Threshold values are between 0 and 1.
Note figures can be closed individually, or all at once by clicking "Close Figures".

Use this URL for a video tutorial on the MATLAB app: https://www.youtube.com/watch?v=1bb0FJu9tZc&feature=youtu.be
Use this URL for a video overview of the Statistics Pipeline: https://www.youtube.com/watch?v=sl2b1jzobso&feature=youtu.be


