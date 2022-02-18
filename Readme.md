# Mistic: image tSNE visualizer

This is a Python tool using the Bokeh library to view multiple multiplex images simultaneously. The code has been tested on 7-panel Vectra TIFF, 32- & 64-panel CODEX TIFF, 16-panel CODEX QPTIFF and 44-panel t-CyCIF TIFF images.

Mistic's GUI with user inputs is shown below:

<img src=/fig_readme/Figure_2.jpg width="80%"></img>

Figure description: A sample Mistic GUI with user inputs is shown. **A.** User-input panel where imaging technique choice, stack montage option or markers can be selected, images borders can be added, new or pre-defined image display coordinates can be chosen, and a theme for the canvases can be selected. **B.** Static canvas showing the image t-SNE colored and arranged as per user inputs. **C.** Live canvas showing the corresponding t-SNE scatter plot where each image is represented as a dot. The live canvas has tabs for displaying additional information per image. Metadata for each image can be obtained by hovering over each dot.


## Features of Mistic
* Two canvases: 
  *   still canvas with the image tSNE rendering 
  *   live canvases with tSNE scatter plots for image metadata rendering
* Dropdown option to select the imaging technique: Vectra, t-CyCIF, or CODEX
* Option to choose between Stack montage view or multiple multiplexed images by selecting the markers to be visualised at once
* Option to place a border around each image based on image metadata
* Option to use a pre-defined tSNE or generate a new set of tSNE co-ordinates
* Option to shuffle images with the tSNE co-ordinates
* Option to render multiple tSNE scatter plots based on image metadata
* Hover functionality available on the tSNE scatter plot to get more information of each image
* Save, zoom, etc each of the Bokeh canvases

## Requirements

* Python 3.6 (code also tested and runs on Python 3.7)
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

* Download this code repository or Open Terminal and use `git clone`

  `$ git clone https://github.com/MathOnco/Mistic.git`
 
* In the Mistic folder, navigate to /user_inputs folder to upload input files:
  * ```Mistic_code/code/user_inputs/```
  * Use the /figures folder to upload the multiplexed images
  * Use the /metadata folder to 
    * Upload the imaging markers of interest as Markers_ids.csv and markers.csv. 
       * Example files are provided in the metadata folder
       * Note: For the Stack Montage option, only the markers.csv file is required. 
    * Optional uploads:
      * Upload image tSNE co-ordinates as X_imagetSNE.csv
        * If no user-generated tSNE co-ordinates are provided, Mistic will generate a set of random coordinates to render the images
      * Upload image metadata such as 
        * Cluster labels as Cluster_categories.csv
        * Patient_ids as Patient_ids.csv
        * Treatments as Treatment_catgories.csv
        * Patient response as Response_categories.csv 
        * If any of these are unavailable, Mistic will use either the randomly-generated or user-provided tSNE points without any color coding i.e. dots are colored in gray.
    * Sample metadata files are provided for reference in the /metadata folder

* The main application is located in the /image_tSNE_GUI folder:
  * ```Mistic_code/code/image_tSNE_GUI/```
* Open a command prompt (or the Terminal application), change to the directory containing /code and type
  *  ```pip install -r requirements.txt``` 
* Upon successful completion of the previous step, Mistic can be run. To run the application, at the command prompt pointing to /code, type
  * ```bokeh serve --port 5098 --show image_tSNE_GUI```
  * This runs a bokeh server locally and will automatically open the interactive dashboard in your browser at http://localhost:5098/image_tSNE_GUI

* For instructions on how to run Mistic on the t-CyCIF data, please check: https://mistic-rtd.readthedocs.io/en/latest/vignette_example_tcycif.html
* For instructions on how to run Mistic on the toy data from our NSCLC Vectra FoVs, please check:https://mistic-rtd.readthedocs.io/en/latest/vignette_example_vectra.html


## Additional information

* Paper on bioRxiv: https://www.biorxiv.org/content/10.1101/2021.10.08.463728v1
* Documentation: https://mistic-rtd.readthedocs.io
* Code has been published at Zenodo: https://doi.org/10.5281/zenodo.5912169 (Mistic Version v1.0.1)
* Toy data is published here: https://doi.org/10.5281/zenodo.6131933
* Mistic is highlighted on Bokeh's user showcase: http://bokeh.org/


