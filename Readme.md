# MIsTic: image tSNE visualizer

This is a Python tool using the Bokeh library to view multiple multiplex images simultaneously. Currently the code supports 7-panel Vectra images.

A sample Bokeh GUI with user inputs is shown below:

<img src=/figures/GUI_Mistic.png width="80%"></img>

## Features of MIsTic
* Two canvases: 
  *   still canvas with the image tSNE rendering 
  *   live canvases with tSNE scatter plots for image metadata rendering
* Checkbox to choose the markers to be visualised at once
* Option to place a border around each image based on image metadata
* Option to use a pre-defined tSNE or generate a new set of tSNE co-ordinates
* Option to shuffle images with the tSNE co-ordinates
* Option to render multiple tSNE scatter plots based on image metadata
* Hover functionality available on the tSNE scatter plot to get more information of each image
* Save, zoom, etc each of the Bokeh canvases

## Requirements

* Python 3.6 (may work on other versions but has not been tested)
  * Install Python from here: https://www.python.org/downloads/
<!---* bokeh 0.12.16 -->
  <!---* For installation, open terminal and type: ``` pip install bokeh ```-->
  <!---* To verify a successful installation, type: ``` bokeh info ```-->
  <!---* Further information: https://docs.bokeh.org/en/latest/docs/first_steps/installation.html -->
* Open a command prompt (or the Terminal application) to download ``` pip ```. Type: 
  * ``` curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py ```
  * ``` python3 get-pip.py ``` and wait for the installation to complete
  * Verify the ``` pip ``` installation by typing ``` python3 get-pip.py ```

## Run the application

* In the Mistic folder, navigate to /user_inputs folder to upload input files:
  * ```Mistic_code/bokeh_GUI/bokeh-branch-2.3/examples/app/user_inputs/```
  * Use the /figures folder to upload the multiplexed images
  * Use the /metadata folder to upload the imaging markers of interest as Markers_ids.csv
  * (Optional) Use the /metadata folder to 
    * Upload image tSNE co-ordinates as X_imagetSNE.csv
    * Upload image metadata such as 
      * Cluster labels as Cluster_categories.csv
      * Patient_ids as Patient_ids.csv
      * Treatments as Treatment_catgories.csv
      * Patent response as Response_categories.csv 
    * Sample metadata files are provided for reference 

* The main application is located in the /app folder:
  * ```Mistic_code/bokeh_GUI/bokeh-branch-2.3/examples/app/```
* Open a command prompt (or the Terminal application), change to the directory containing /app and type
  *  ```pip install -r requirements.txt``` 
* Upon successful completion of the previous step, Mistic can be run. To run the application, at the command prompt pointing to /app, type
  * ```bokeh serve --port 5098 --show image_tSNE_GUI```
  * This runs a bokeh server locally and will automatically open the interactive dashboard in your browser at http://localhost:5098/image_tSNE_GUI
