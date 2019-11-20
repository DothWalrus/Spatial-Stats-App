# Spatial-Stats-App
App used to locate dendritic cell nuclei in images of gastric mucosa.

Installation and Usage:
  Download the dendritic_segmentation_Integration.mlapp file.
  Log into MATLAB Online on a Chrome browser and upload the file to your MATLAB online repository.
  Upload the images and their associated Selection.csv files to your MATLAB online repository.
  Open the folder in your MATLAB repository in which you wish to save the generated MAT files to.
  Run the app be right clicking the dendritic_segmentation_Integration.mlapp file and clicking run.
  Click the "Select TIF" button and choose the desired TIF to be analyzed.
  Click the "Select CSV" button and choose the CSV corresponding to the TIF.
  Click the "Display TIF Channels" button to see each channel in the TIF.
  Set the "Nuclei Channel" and "DC Channel" spinners to the corresponding channels as seen in the TIF channels.
  Click "Run" and note the figure generated is a visualization of the image analysis pipeline.
  Ensure that a MAT file has been generated in your "Current Folder".
  Checking "debug" will generate 11 figures of various parts of the image analysis.
  If necessary adjust the "Threshold" slider to more or less severe thresholding to get the correct nuclei.
