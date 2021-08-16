### Code for Mistic 
### Image t-SNE viewer for multiplexed images
### 
### 
### Jan 29th 2021
### code author SP
### 25th June 2021
### github version

### to do 
#
# work on paths
# give 10 images vachu code
# user input area
# python installation - try on your machine 

import os
import sys
import random
import warnings

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from matplotlib.figure import Figure 
#from mpl_toolkits.axisartist.axislines import Subplot

#from tqdm import tqdm
#from itertools import chain

from matplotlib import cm

#pip install --ignore-installed --upgrade tensorflow



from scipy import ndimage

from sklearn.manifold import TSNE

import phenograph as pg
from scipy.stats import zscore

from skimage.color import rgb2gray

from matplotlib.patches import Polygon

from skimage import io

from skimage import data
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb
from skimage.io import imread, imshow, imread_collection, concatenate_images
from skimage.transform import resize
from skimage.morphology import label

import seaborn as sns

import scipy.spatial
#import pysal

#from scipy.spatial import ConvexHull, convex_hull_plot_2d

#import multiview as mv

from scipy.spatial import distance
#from multiview.mvsc import MVSC

import seaborn as sns
#from statannot import add_stat_annotation

import matplotlib.image as mpimg 
#import time as time

#import random

import matplotlib.colors as colors 

from matplotlib import colors as mcolors

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

from PIL import Image, ImageOps

from bokeh.layouts import column, row
from bokeh.models import (Select, Button)
from bokeh.palettes import Spectral5
from bokeh.plotting import curdoc, figure
from bokeh.models.widgets import Div
from bokeh.layouts import column, layout

from bokeh.models import (HoverTool, ColumnDataSource)
from bokeh.models.widgets import RadioButtonGroup
from bokeh.models import CheckboxGroup
from bokeh.models import CustomJS
from bokeh.models import BoxSelectTool
from bokeh.themes import Theme

from bokeh.models.widgets import Panel, Tabs
from bokeh.io import output_file, show


##################################
#### Section that gets populated #
####based on User uploads ########
##################################
##################################
## This section has path to images, 
## outputs, markers in the data,
## and user choices
#############################
#############################


path_wd = os.getcwd()
print(path_wd)
#path_python = '/Users/4470526/Downloads/Projects/Amer_images'

# tSNE co-ordinates input
df_Xtsne = pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE.csv'), index_col=None, header= None)
tsne = np.array(df_Xtsne )
tx, ty = tsne[:,0], tsne[:,1]
tx = (tx-np.min(tx)) / (np.max(tx) - np.min(tx))
ty = (ty-np.min(ty)) / (np.max(ty) - np.min(ty))
num_images = tsne.shape[0]
# cluster assignments input
#df_clustasgn = pd.read_csv(os.path.join(path_wd + '/user_inputs/metadatacluster_asgn_imtsne.csv'), index_col=None, header= None)
#df_clustasgn.shape
#cluster_asgn = np.array(df_clustasgn )

#LABELS_MARKERS = ['Marker 1', 'Marker 2', 'Marker 3', 'Marker 4', 'Marker 5', 'Marker 6']
df_markers = pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/Marker_ids.csv'), index_col=None, header= None)
markers_list = np.array(df_markers ).flatten()
LABELS_MARKERS = []
for j in range(markers_list.shape[0]):
    LABELS_MARKERS.append(markers_list[j]) 
   
colours_58 = ["firebrick","gold","royalblue","green","dimgray","orchid","darkviolet",
              "red", "orange", "limegreen", "blue", "purple", "seagreen","gold","darkolivegreen",
              "lightpink","thistle","mistyrose","saddlebrown","slategrey","powderblue",
            "palevioletred","mediumvioletred","yellowgreen","lemonchiffon","chocolate",
              "lightsalmon","lightcyan","lightblue", "darkorange","black","darkblue","darkgreen","paleturquoise","yellow","rosybrown",
             "steelblue","dodgerblue","darkkhaki","lime","coral","aquamarine","mediumpurple","violet","plum",
             "deeppink","navy","seagreen","teal","mediumspringgreen","cadetblue",
             "maroon","silver","sienna","crimson","slateblue","magenta","darkmagenta"]

colours_resp = ["yellow","red","green"]
colours_tx = ["orange","limegreen","violet"]

## read in images
marker_image_list = []
pat_fov_list = []
FoV_path = os.path.join(path_wd + '/user_inputs/figures/')
for fname in os.listdir(FoV_path):
    print(fname)
    marker_image_list.append(FoV_path+fname)
    pat_fov_list.append(fname)
                            
#for tooltip
#print(pat_id_mt)
pat_ind_list = np.array(pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/Patient_ids.csv'),header= None,index_col=None)).flatten()
pat_ind_list = pat_ind_list.tolist()

# collecting the response metadata
color_vec = []
resp_list= []
resp_list_1 = np.array(pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/Response_categories.csv'),header= None,index_col=None)).flatten()
uni_resp,counts_resp = np.unique(resp_list_1,return_counts=True)

for i in range(resp_list_1.shape[0]):
    row_t = np.int(np.array(np.where(resp_list_1[i]==uni_resp)).flatten())
    #if resp_list_1[i]=='Response 1':
    color_vec.append(colours_resp[row_t])
    #elif resp_list_1[i]=='Response 2':
    #    color_vec.append('red')
    resp_list.append(resp_list_1[i])

        
# collecting the treatment metadata
color_vec_tx = []
tx_list = []
tx_list_1= np.array(pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/Treatment_categories.csv'),header= None,index_col=None)).flatten()
uni_tx,counts_tx = np.unique(tx_list_1,return_counts=True)
for i in range(tx_list_1.shape[0]):
    row_t = np.int(np.array(np.where(tx_list_1[i]==uni_tx)).flatten())
    #if resp_list_1[i]=='Response 1':
    color_vec_tx.append(colours_tx[row_t])
    #if tx_list_1[i]=='Treatment 1':
    #    color_vec_tx.append('orange')
    #elif tx_list_1[i]=='Treatment 2':
    #    color_vec_tx.append('limegreen')
    tx_list.append(tx_list_1[i])

#print(tx_list)
# collecting the cluster assignments metadata
color_vec_clasgn = []
cluster_anno_list = []
clust_asgn_list = []
clust_asgn_list_1 = np.array(pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/Cluster_categories.csv'),
                                         header= None,index_col=None)).flatten()
for i in range(clust_asgn_list_1.shape[0]):
   
    color_vec_clasgn.append(colours_58[clust_asgn_list_1[i]])
    clust_asgn_list.append(clust_asgn_list_1[i])
    cluster_anno_list.append('Cluster '+ str(clust_asgn_list_1[i]))  

     

        
#########
#########

# set up widgets


RB_1 = ['Yes', 'No']

RB_2 = ['Generate new co-ordinates', 'Use pre-defined co-ordinates', 'Arrange in rows']

RB_3 = ['Yes', 'No']

TS_1 = ['black','gray', 'dark blue']

TOOLS="hover,pan,crosshair,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,save,box_select,"

x_range=(0,1)
y_range=(0,1)


#main image tSNE canvas
p = figure(tools=TOOLS,x_range=x_range, y_range=y_range,width=1000,height=1000)
tsne_points = np.zeros([1,2])

#additional tSNE scatter plot canvas



point_tSNE = figure(plot_width=350, plot_height=350,
              tools='hover,pan,wheel_zoom,box_select,reset')
point_tSNE.title.text = 'tSNE point cloud for Patient response'


point_tSNE.scatter(tsne[:,0],tsne[:,1],fill_alpha=0.6, color ='red',size=8,legend='Response')

point_tSNE.legend.location = "bottom_left"

theme_black = Theme(json={
    'attrs': {
        'Figure': {
            'background_fill_color': '#2F2F2F',
            'border_fill_color': '#2F2F2F',
            'outline_line_color': '#444444'
            },
        'Axis': {
            'axis_line_color': "white",
            'axis_label_text_color': "white",
            'major_label_text_color': "white",
            'major_tick_line_color': "white",
            'minor_tick_line_color': "white",
            'minor_tick_line_color': "white"
            },
        'Grid': {
            'grid_line_dash': [6, 4],
            'grid_line_alpha': .3
            },
        'Circle': {
            'fill_color': 'lightblue',
            'size': 10,
            },
        'Title': {
            'text_color': "white"
            }
        }
    })

theme_gray = Theme(json={
    'attrs': {
        'Figure': {
            'background_fill_color': '#555555',
            'border_fill_color': '#2F2F2F',
            'outline_line_color': '#444444'
            },
        'Axis': {
            'axis_line_color': "white",
            'axis_label_text_color': "white",
            'major_label_text_color': "white",
            'major_tick_line_color': "white",
            'minor_tick_line_color': "white",
            'minor_tick_line_color': "white"
            },
        'Grid': {
            'grid_line_dash': [6, 4],
            'grid_line_alpha': .3
            },
        'Circle': {
            'fill_color': 'lightblue',
            'size': 10,
            },
        'Title': {
            'text_color': "white"
            }
        }
    })

theme_blue = Theme(json={
    'attrs': {
        'Figure': {
            'background_fill_color': '#25256d',
            'border_fill_color': '#2F2F2F',
            'outline_line_color': '#444444'
            },
        'Axis': {
            'axis_line_color': "white",
            'axis_label_text_color': "white",
            'major_label_text_color': "white",
            'major_tick_line_color': "white",
            'minor_tick_line_color': "white",
            'minor_tick_line_color': "white"
            },
        'Grid': {
            'grid_line_dash': [6, 4],
            'grid_line_alpha': .3
            },
        'Circle': {
            'fill_color': 'lightblue',
            'size': 10,
            },
        'Title': {
            'text_color': "white"
            }
        }
    })

#########
#########


## 8th March 2021
TOOLS="hover,pan,crosshair,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"


TOOLTIPS = [
        ("index", "$index"),
        ("(x,y)", "($x, $y)"),
        ("Pat_id", "@pat_list"),
        ("response", "@res_list"),
        ("FoV","@fov_list")
    ]
    
p1 = figure(plot_width=400, plot_height=400, tooltips=TOOLTIPS,tools = TOOLS,
               title="Patient response")
'''
## added to make source available outside the ()
source = ColumnDataSource(data=dict(
    x=tsne[:,0],
    y=tsne[:,1],
    pat_list = pat_id_mt[:,0],#np.array(pat_ind_list), #[item for sublist in pat_ind_list for item in sublist]
    res_list = resp_list,#[item for sublist in resp_list for item in sublist]
    fov_list = pat_fov_list,
    color_vec_list = color_vec
    #legend_p11 = legend_p1
))
'''

def draw_tSNE_scatter(tsne1):
    tsne=np.asarray(tsne1)
    source = ColumnDataSource(data=dict(
        x=tsne[:,0],
        y=tsne[:,1],
        pat_list = pat_ind_list, #mt[:,0],#np.array(pat_ind_list), #[item for sublist in pat_ind_list for item in sublist]
        res_list = resp_list,#[item for sublist in resp_list for item in sublist]
        fov_list = pat_fov_list,
        color_vec_list = color_vec,
        tx_list = tx_list,
        color_vec_tx_list = color_vec_tx,
        clust_asgn_list = clust_asgn_list,
        color_vec_clasgn_list = color_vec_clasgn,
        cluster_anno_list = cluster_anno_list
        #legend_p11 = legend_p1
    ))
    TOOLS="hover,pan,crosshair,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"

    TOOLTIPS = [
        ("index", "$index"),
        ("(x,y)", "($x, $y)"),
        ("Pat_id", "@pat_list"),
        ("response", "@res_list"),
        ("FoV","@fov_list")
    ]




    p1 = figure(plot_width=400, plot_height=400, tooltips=TOOLTIPS,tools = TOOLS,
               title="Patient response")
    #p1.circle('x', 'y', size=10, source=source)
    #26th Feb 2021
    p1.scatter('x', 'y', size=10, source=source, legend= 'res_list', color = 'color_vec_list',fill_alpha=0.6)
    p1.legend.location = "bottom_left"

    #8th June 2021
    p2 = figure(plot_width=400, plot_height=400, tooltips=TOOLTIPS,tools = TOOLS,
               title="Treatment category")
    p2.scatter('x', 'y', size=10, source=source, legend= 'tx_list', color = 'color_vec_tx_list',fill_alpha=0.6)
    p2.legend.location = "bottom_left"
    
    #8th June 2021
    p3 = figure(plot_width=400, plot_height=400, tooltips=TOOLTIPS,tools = TOOLS,
               title="Cluster annotations")
    p3.scatter('x', 'y', size=10, source=source, legend= 'cluster_anno_list', color = 'color_vec_clasgn_list',fill_alpha=0.6)
    p3.legend.location = "bottom_left"
    #p3.legend.label_text_font_size = "x-small"
    
    
    return ([p1,p2,p3, source])


###

                            
### image-tSNE

width = 5000 #15000
height = 5000
max_dim = 500
tw = 1344
th = 1008

full_image_1 = Image.new('RGBA', (width, height))


#####code from  https://scipython.com/blog/poisson-disc-sampling-in-python/ ####



def get_cell_coords(pt,a):
    """Get the coordinates of the cell that pt = (x,y) falls in."""

    return int(pt[0] // a), int(pt[1] // a)

def get_neighbours(coords,nx,ny,cells):
    """Return the indexes of points in cells neighbouring cell at coords.
    For the cell at coords = (x,y), return the indexes of points in the cells
    with neighbouring coordinates illustrated below: ie those cells that could 
    contain points closer than r.
                                     ooo
                                    ooooo
                                    ooXoo
                                    ooooo
                                     ooo
    """

    dxdy = [(-1,-2),(0,-2),(1,-2),(-2,-1),(-1,-1),(0,-1),(1,-1),(2,-1),
            (-2,0),(-1,0),(1,0),(2,0),(-2,1),(-1,1),(0,1),(1,1),(2,1),
            (-1,2),(0,2),(1,2),(0,0)]
    neighbours = []
    for dx, dy in dxdy:
        neighbour_coords = coords[0] + dx, coords[1] + dy
        if not (0 <= neighbour_coords[0] < nx and
                0 <= neighbour_coords[1] < ny):
            # We're off the grid: no neighbours here.
            continue
        neighbour_cell = cells[neighbour_coords]
        if neighbour_cell is not None:
            # This cell is occupied: store this index of the contained point.
            neighbours.append(neighbour_cell)
    return neighbours

def point_valid(pt,a,nx,ny,cells,samples,r):
    """Is pt a valid point to emit as a sample?
    It must be no closer than r from any other point: check the cells in its
    immediate neighbourhood.
    """

    cell_coords = get_cell_coords(pt,a)
    for idx in get_neighbours(cell_coords,nx,ny,cells):
        nearby_pt = samples[idx]
        # Squared distance between or candidate point, pt, and this nearby_pt.
        distance2 = (nearby_pt[0]-pt[0])**2 + (nearby_pt[1]-pt[1])**2
        if distance2 < r**2:
            # The points are too close, so pt is not a candidate.
            return False
    # All points tested: if we're here, pt is valid
    return True

def get_point(k, refpt,r,a,nx,ny,cells,samples):
    """Try to find a candidate point relative to refpt to emit in the sample.
    We draw up to k points from the annulus of inner radius r, outer radius 2r
    around the reference point, refpt. If none of them are suitable (because
    they're too close to existing points in the sample), return False.
    Otherwise, return the pt.
    """
    i = 0
    while i < k:
        rho, theta = np.random.uniform(r, 2*r), np.random.uniform(0, 2*np.pi)
        pt = refpt[0] + rho*np.cos(theta), refpt[1] + rho*np.sin(theta)
        if not (0 <= pt[0] < width and 0 <= pt[1] < height):
            # This point falls outside the domain, so try again.
            continue
        if point_valid(pt,a,nx,ny,cells,samples,r):
            return pt
        i += 1
    # We failed to find a suitable point in the vicinity of refpt.
    return False




    
    
    
##################
    

def generate_image_tSNE(chk_box_marker,rb_val,rb_rs_val,rb_shf_val, LABELS_MARKERS):
    full_image = Image.new('RGBA', (width, height))
    size = [256,256]

    if (rb_shf_val == 0):
        #create the random tsne projections
        if (rb_rs_val==1):
            # Choose up to k points around each reference point as candidates for a new
            # sample point
            k = 10

            # Minimum distance between samples
            r = 1.7

            width_1, height_1 = 20, 22

            print('Generating random co-ordinates')

            # Cell side length
            a = r/np.sqrt(2)
            # Number of cells in the x- and y-directions of the grid
            nx, ny = int(width_1 / a) + 1, int(height_1 / a) + 1

            # A list of coordinates in the grid of cells
            coords_list = [(ix, iy) for ix in range(nx) for iy in range(ny)]
            # Initilalize the dictionary of cells: each key is a cell's coordinates, the
            # corresponding value is the index of that cell's point's coordinates in the
            # samples list (or None if the cell is empty).
            cells = {coords: None for coords in coords_list}



            # Pick a random point to start with.
            pt = (np.random.uniform(0, width_1), np.random.uniform(0, height_1))
            samples = [pt]
            # Our first sample is indexed at 0 in the samples list...
            cells[get_cell_coords(pt,a)] = 0
            # ... and it is active, in the sense that we're going to look for more points
            # in its neighbourhood.
            active = [0]

            nsamples = 1
            # As long as there are points in the active list, keep trying to find samples.
            while (nsamples < num_images): #active:
                # choose a random "reference" point from the active list.
                idx = np.random.choice(active)
                refpt = samples[idx]
                # Try to pick a new point relative to the reference point.
                pt = get_point(k, refpt,r,a,nx,ny,cells,samples)
                if pt:
                    # Point pt is valid: add it to the samples list and mark it as active
                    samples.append(pt)
                    nsamples += 1
                    active.append(len(samples)-1)
                    cells[get_cell_coords(pt,a)] = len(samples) - 1
                    print('nsamples is: ',str(nsamples))
                else:
                    # We had to give up looking for valid points near refpt, so remove it
                    # from the list of "active" points.
                    active.remove(idx)

            tsne = np.asarray(samples)
            df_Xtsne = pd.DataFrame(tsne)
            df_Xtsne.to_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE_user.csv'), header=None, index=None)


        elif(rb_rs_val==0):

            df_Xtsne = pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE.csv'), index_col=None, header= None)
            df_Xtsne.shape
            tsne = np.array(df_Xtsne )
            
        elif(rb_rs_val==2):
            
            # Choose up to k points around each reference point as candidates for a new
            # sample point
            k = 10

            # Minimum distance between samples
            r = 2

            width_1, height_1 = 20, 22

            print('Arrange images side-by-side')

            # Cell side length
            a = r/np.sqrt(2)
            # Number of cells in the x- and y-directions of the grid
            nx, ny = int(width_1 / a) + 1, int(height_1 / a) + 1

            # A list of coordinates in the grid of cells
            coords_list = [(ix, iy) for ix in range(nx) for iy in range(ny)]
           
        
            #logic to get grid points - 29th March 2021
            m = np.int(np.floor(nx*ny/num_images))
            #print('bbbbb')
            #print(m)
            
            row_needed = []
            def multiples(m, num_images):
                for i in range(num_images):
                    row_needed.append(i*m)
        
            multiples(m,num_images)
            
            #print('aaaa')
            #print(num_images)
            #print(np.array(row_needed).flatten())
            select_coords = np.array(coords_list)[np.array(row_needed).flatten()]
            
            print(type(select_coords))
            ################
            
            #tsne = np.asarray(coords_list)
            tsne = select_coords
            df_Xtsne = pd.DataFrame(tsne)
            df_Xtsne.to_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE_seeall.csv'), header=None, index=None)

            
        # 26th spril 2021
        # save the tSNE points to bve read later on irrespective of whoch option was chosen
        df_Xtsne_touse = pd.DataFrame(tsne)
        df_Xtsne_touse.to_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE_touse.csv'), header=None, index=None)
            
        tx, ty = tsne[:,0], tsne[:,1]
        tx = (tx-np.min(tx)) / (np.max(tx) - np.min(tx))
        ty = (ty-np.min(ty)) / (np.max(ty) - np.min(ty))
    

        
        ##identify the markers - 26th March 2021
        mm = np.asarray(LABELS_MARKERS)[chk_box_marker]
        
        tiles = []
        mc = []
        wc = []
        
        if LABELS_MARKERS[0] in mm:
            mc.append(0)
            wc.append(800)
        if LABELS_MARKERS[1] in mm:
            mc.append(1)
            wc.append(400)
        if LABELS_MARKERS[2] in mm:
            mc.append(2)
            wc.append(400)
        if LABELS_MARKERS[3] in mm:
            mc.append(3)
            wc.append(100)
        if LABELS_MARKERS[4] in mm:
            mc.append(4)
            wc.append(100)
        if LABELS_MARKERS[5] in mm:
            mc.append(5)
            wc.append(100)


        marker_choice = np.array(mc)
        weight_choice = np.array(wc)


        
        inds = range(len(marker_image_list))

        N = len(inds)

        
        for k in range(len(inds)):
            #print('KKKKK', k)
            image_all_2 = []
            image_all_1 = []
           

            
            im = io.imread(marker_image_list[inds[k]])#,plugin='matplotlib')
            #print (im.shape)    
            
            for m in range(len(marker_choice)):
                image = im[marker_choice[m]]

                median_filtered = scipy.ndimage.median_filter(image, size=1)
                image_all_2.append(median_filtered)


                # apply threshold

                thresh = threshold_otsu(image_all_2[m])





                print(thresh)
                #thresh = 12
                bw = closing(image > thresh, square(1))

                # remove artifacts connected to image border
                cleared_imtsne = clear_border(bw)


                image_all_1.append(cleared_imtsne*weight_choice[m]) #10

            tl = sum(image_all_1)



            tile_1 = Image.fromarray(np.uint8(cm.viridis(tl)*255))




            old_size = tile_1.size
            #print('resp')
            #print(pat_ind_list[inds[k]])

            #print(rb_val)

            if(rb_val==1): # for the border
                new_size = (old_size[0]+30, old_size[1]+30)


   
                if(resp_list[inds[k]] == 'Response 1'):
                    new_im = Image.new("RGB", new_size,color='yellow')
                elif(resp_list[inds[k]] == 'Response 2'):
                    new_im = Image.new("RGB", new_size,color='red')

            elif(rb_val==0): 
                new_size = (old_size[0]+1, old_size[1]+1)
                new_im = Image.new("RGB", new_size,color='black')

            new_im.paste(tile_1, (int((new_size[0]-old_size[0])/2),int((new_size[1]-old_size[1])/2)))



            file_name = os.path.join(path_wd+'/output_tiles/image_tsne_tils'+str(k)+'.png')



            new_im.save(file_name) #2nd Feb 2021

            # Read Images 


            tile_2 = Image.open(file_name)
            rs = max(1, tile_2.width/max_dim, tile_2.height/max_dim)
            tile_2 = tile_2.resize((int(tile_2.width/rs), int(tile_2.height/rs)), Image.ANTIALIAS) #commented antialias on 9th Feb 2021

            full_image.paste(tile_2, (int((width-max_dim)*ty[k]), int((height-max_dim)*tx[k])),mask=tile_2.convert('RGBA'))



        matplotlib.pyplot.figure(figsize = (25,20))
        plt.imshow(full_image)
        plt.axis('off')
        #return (full_image)



        ###5th March 2021
        k = random.randint(2,500)
        file_name = os.path.join(path_wd+'/image_tSNE_GUI/static/image_tsne_tils_all_'+str(k)+'.png')
        full_image.save(file_name)
        #file_name_1 = "/Users/4470526/Google Drive/Project/Projects_Moffitt/NSCLC_grant_aim2/image_tSNE_code_copy/bokeh_GUI/bokeh-branch-2.3/examples/app/image_tSNE_GUI/static/image_tsne_tils_all_"+str(k)+".png"
        #full_image.save(file_name_1)

 
        rotated_img     = full_image.rotate(90)
        file_name_rot = os.path.join(path_wd+'/image_tSNE_GUI/static/image_tsne_tils_all_'+str(k)+'_rot.png')#"/Users/4470526/Google Drive/Project/Projects_Moffitt/NSCLC_grant_aim2/image_tSNE_code_copy/bokeh_GUI/bokeh-branch-2.3/examples/app/image_tSNE_GUI/static/image_tsne_tils_all_"+str(k)+"_rot.png" 


        rotated_img.save(file_name_rot)
        
        
        
        
    ##code for reshuffling
    
    elif (rb_shf_val==1): # shuffle images 
        
        #pick the tsne file to start with
        if(rb_rs_val==1):
            df_Xtsne = pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE_user.csv'), index_col=None, header= None)
            print('usertsne')
        elif(rb_rs_val==0):
            df_Xtsne = pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE.csv'), index_col=None, header= None)
            print('origtsne')
            
        elif(rb_rs_val==2):
            df_Xtsne = pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE_seeall.csv'), index_col=None, header= None)
            print('seealltsne')
            
        tsne = np.array(df_Xtsne)   
        
        
        # 26th spril 2021
        # save the tSNE points to be read later on, irrespective of which option was chosen
        df_Xtsne_touse = pd.DataFrame(tsne)
        df_Xtsne_touse.to_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE_touse.csv'), header=None, index=None)
        
        
        tx, ty = tsne[:,0], tsne[:,1]
        tx = (tx-np.min(tx)) / (np.max(tx) - np.min(tx))
        ty = (ty-np.min(ty)) / (np.max(ty) - np.min(ty))
            
        full_image = Image.new('RGBA', (width, height))
        
        size = [256,256]

        
        ##identify the markers - 26th March 2021
        mm = np.asarray(LABELS_MARKERS)[chk_box_marker]
        
        tiles = []
        mc = []
        wc = []
        
        if LABELS_MARKERS[0] in mm:
            mc.append(0)
            wc.append(800)
        if LABELS_MARKERS[1] in mm:
            mc.append(1)
            wc.append(400)
        if LABELS_MARKERS[2] in mm:
            mc.append(2)
            wc.append(400)
        if LABELS_MARKERS[3] in mm:
            mc.append(3)
            wc.append(100)
        if LABELS_MARKERS[4] in mm:
            mc.append(4)
            wc.append(100)
        if LABELS_MARKERS[5] in mm:
            mc.append(5)
            wc.append(100)
    
                                
        marker_choice = np.array(mc)#[3,5]
        weight_choice = np.array(wc)#[800,100]

        inds = range(len(marker_image_list))

        N = len(inds)
        shf_N = np.array(range(N))     
        random.shuffle(shf_N)
        
        for k in range(N):
            #print('KKKKK', k)
            image_all_2 = []
            image_all_1 = []
           

            
            im = io.imread(marker_image_list[shf_N[k]])#,plugin='matplotlib')
            print (im.shape)    
         
            for m in range(len(marker_choice)):
                image = im[marker_choice[m]]

                median_filtered = scipy.ndimage.median_filter(image, size=1)
                image_all_2.append(median_filtered)



                # apply threshold
                '''
                if (marker_choice[m] == 4):
                    if (pat_ind_list[k]==43):
                        thresh = 11
                    elif (pat_ind_list[k]==41):
                        thresh = 6
                    elif (pat_ind_list[k]==54):
                        thresh = 2
                    elif (pat_ind_list[k]==36):
                        thresh = 2
                    elif (pat_ind_list[k]==37):
                        thresh = 3
                    elif (pat_ind_list[k]==34):
                        thresh = 4
                    elif (pat_ind_list[k]==47):
                        thresh = 14
                    elif (pat_ind_list[k]==44):
                        thresh = 1.44
                    elif (pat_ind_list[k]==33):
                        thresh = 12
                elif (marker_choice[m] == 5):
                    if (pat_ind_list[k]==43):
                        thresh = 3
                    elif (pat_ind_list[k]==41):
                        thresh = 7
                    elif (pat_ind_list[k]==54):
                        thresh = 3.4
                    elif (pat_ind_list[k]==36):
                        thresh = 3.77
                    elif (pat_ind_list[k]==37):
                        thresh = 5
                    elif (pat_ind_list[k]==34):
                        thresh = 3.2
                    elif (pat_ind_list[k]==47):
                        thresh = 8
                    elif (pat_ind_list[k]==44):
                        thresh = 4
                    elif (pat_ind_list[k]==33):
                        thresh = 11
                else :
                '''
                thresh = threshold_otsu(image_all_2[m])

                                
                print(thresh)
            
                bw = closing(image > thresh, square(1))

                # remove artifacts connected to image border
                cleared_imtsne = clear_border(bw)


                image_all_1.append(cleared_imtsne*weight_choice[m]) #10

            tl = sum(image_all_1)
            tile_1 = Image.fromarray(np.uint8(cm.viridis(tl)*255))
            old_size = tile_1.size

            print(pat_ind_list[shf_N[k]])


            

            if(rb_val==1): # for the border
                new_size = (old_size[0]+30, old_size[1]+30)

                #new_im = Image.new("RGB", new_size,color='yellow')   ## luckily, this is already black!

                if(resp_list[shf_N[k]] == 'Response 1'):
                    new_im = Image.new("RGB", new_size,color='yellow')
                elif(resp_list[shf_N[k]] == 'Response 2'):
                    new_im = Image.new("RGB", new_size,color='red')

            elif(rb_val==0): 
                new_size = (old_size[0]+1, old_size[1]+1)
                new_im = Image.new("RGB", new_size,color='black')

            new_im.paste(tile_1, (int((new_size[0]-old_size[0])/2),int((new_size[1]-old_size[1])/2)))                     
                                
            file_name = os.path.join(path_wd+'/output_tiles/image_tsne_tils'+str(k)+'.png')


            #plt.imsave(file_name,tile_1)
            new_im.save(file_name) #2nd Feb 2021


            # Read Images 


            tile_2 = Image.open(file_name)
            rs = max(1, tile_2.width/max_dim, tile_2.height/max_dim)
            tile_2 = tile_2.resize((int(tile_2.width/rs), int(tile_2.height/rs)), Image.ANTIALIAS) #commented antialias on 9th Feb 2021


            full_image.paste(tile_2, (int((width-max_dim)*ty[shf_N[k]]), int((height-max_dim)*tx[shf_N[k]])),mask=tile_2.convert('RGBA'))


        matplotlib.pyplot.figure(figsize = (25,20))
        plt.imshow(full_image)
        plt.axis('off')
        #return (full_image)



        ###5th March 2021
        k = random.randint(2,500)
        file_name = os.path.join(path_wd+'/image_tSNE_GUI/static/image_tsne_tils_all_'+str(k)+'.png')
        #file_name = "/Users/4470526/Downloads/Projects/Amer_images/output/output_figures/Boundary_runs/im_tsne_tiles/tils_CSBC/image_tsne_tils_all_"+str(k)+".png"
        full_image.save(file_name)
        #file_name_1 = "/Users/4470526/Google Drive/Project/Projects_Moffitt/NSCLC_grant_aim2/image_tSNE_code_copy/bokeh_GUI/bokeh-branch-2.3/examples/app/image_tSNE_GUI/static/image_tsne_tils_all_"+str(k)+".png"
        #full_image.save(file_name_1)

        rotated_img     = full_image.rotate(90)
        file_name_rot = file_name = os.path.join(path_wd+'/image_tSNE_GUI/static/image_tsne_tils_all_'+str(k)+'_rot.png')
        #"/Users/4470526/Google Drive/Project/Projects_Moffitt/NSCLC_grant_aim2/image_tSNE_code_copy/bokeh_GUI/bokeh-branch-2.3/examples/app/image_tSNE_GUI/static/image_tsne_tils_all_"+str(k)+"_rot.png" 


        rotated_img.save(file_name_rot)


                                
                                
                                
                                
                                
                                
    return([file_name_rot,tsne])





def button_callback():
    
    theme_t = theme_select.value
    if(theme_t == 'black'):
        curdoc().theme= theme_black
    elif(theme_t == 'gray'):
        curdoc().theme= theme_gray
    elif(theme_t == 'dark blue'):
        curdoc().theme= theme_blue
        
    jk = create_figure()
    layout.children[1] = jk[0]
    
    print('############')
    print('In Button callback')
    tsne3 = jk[1]
    print(type(tsne3))
    print(tsne3)

    p1_out = draw_tSNE_scatter(tsne3)
    p1 = p1_out[0]
    p2 = p1_out[1]
    p3 = p1_out[2]
    #layout.children[2] = p1
    source = p1_out[3]
    
    tab1 = Panel(child=p1, title="Response")
    tab2 = Panel(child=p2, title="Treatment")
    tab3 = Panel(child=p3, title="Cluster Annotations")

    tabs = Tabs(tabs=[ tab1, tab2, tab3 ])
    layout.children[2] = tabs
    
    return([p,p1,p2,p3,source,tabs]) 






def button2_callback():
    
    
    
    rb = radio_button_group.value
    
    rb_rs = radio_button_group_RS.value
    
    rb_shf = radio_button_group_Shf.value
    
    chk_box_marker = checkbox_group.active
    
    print(chk_box_marker)
    print(rb_shf)
    print(rb_rs)
    print(rb)
    
    #pk = get_cds_source()
   
    #source1 = pk


    #p1 = draw_tSNE_scatter(tsne_b2)
        
    print('######cds#####')
    print(source.data)
    #print('######cds#points####')
    #print(cds.selected_points)
    print('######cds#rows####')
    print(source.selected.indices)
    

    display_div.text = ""
    for selected_row in source.selected.indices:
        print(selected_row)
        display_div.text += f"Row {selected_row}: <br>"
    p1 = figure(tools="hover ,reset ,pan ,wheel_zoom ,box_select ,lasso_select")
    #layout.children[2] = create_figure2(source)
    return p1    
    


def create_figure():
   
    p = figure(tools=TOOLS,x_range=x_range, y_range=y_range,width=1000,height=1000)
    
    rb = radio_button_group.value
    
    rb_rs = radio_button_group_RS.value
    
    rb_shf = radio_button_group_Shf.value
    
    chk_box_marker = checkbox_group.active
    
    print(chk_box_marker)

    
    if(rb=='Yes'):
        rb_val = 1
    else:
        rb_val = 0
        
    if(rb_rs=='Generate new co-ordinates'):
        rb_rs_val = 1
    elif(rb_rs=='Use pre-defined co-ordinates'):
        rb_rs_val = 0
    elif(rb_rs == 'Arrange in rows'):
        rb_rs_val = 2
        
        
    if(rb_shf=='Yes'):
        rb_shf_val = 1
    else:
        rb_shf_val = 0
            

    
    if(len(chk_box_marker)==0):
        #file_name_1 = os.path.join(path_wd+'/user_inputs/image_tsne_tils_all_1_rot.png')#"image_tSNE_GUI/static/image_tsne_tils_all.png"
        file_name_1 = "image_tSNE_GUI/static/image_tsne_tils_all.png"
        p.image_url(url=[file_name_1], x=x_range[0],y=y_range[1],w=x_range[1]-x_range[0],h=y_range[1]-y_range[0])
        df_Xtsne = pd.read_csv(os.path.join(path_wd + '/user_inputs/metadata/X_imagetSNE.csv'), index_col=None, header= None)
        df_Xtsne.shape
        tsne_points = np.array(df_Xtsne )

    else:
    
 
        out1 = generate_image_tSNE(chk_box_marker, rb_val,rb_rs_val,rb_shf_val,LABELS_MARKERS)
        file_name_all = out1[0]
        tsne_points = out1[1]


        file_name_1 = file_name_all.split('app')[1]
        print('#########file_name_1########')  
        print(file_name_1)


        p.image_url(url=[file_name_1], x=x_range[0],y=y_range[1],w=x_range[1]-x_range[0],h=y_range[1]-y_range[0])


    return ([p,tsne_points])




desc = Div(text=open(os.path.join(path_wd + '/image_tSNE_GUI/desc.html')).read(), sizing_mode="stretch_width")

#desc = Div(text=open(os.path.join(path_wd + '/image_tSNE_GUI/desc.html').read(), sizing_mode="stretch_width")
radio_button_group = Select(value='Yes',
                          title='Image border',
                          width=200,
                          options=RB_1)

radio_button_group_RS = Select(value='Generate new co-ordinates',
                          title='tSNE co-ordinates',
                          width=220,
                          options=RB_2)

radio_button_group_Shf = Select(value='Yes',
                          title='Shuffle images',
                          width=200,
                          options=RB_3)

theme_select = Select(value = 'black',
                      title='Theme', 
                      width = 200, 
                      options = TS_1)


checkbox_group = CheckboxGroup(labels=LABELS_MARKERS)#, active=[0, 1])


button = Button(label='Run', width=100, button_type="success")

## added to activate Run button
button.on_click(button_callback)



button2 = Button(label='Run2', width=100, button_type="success")


## added to activate Run button
##
button2.on_click(button2_callback)


print('########### new p1#')



######

# set up layout
out2 = create_figure()
print(out2)
p = out2[0]
tsne2 = out2[1]
print('############')
print(type(tsne2))
print(tsne2)
p11_out = draw_tSNE_scatter(tsne2)
p1 = p11_out[0]
p2 = p11_out[1]
p3 = p11_out[2]

tab1 = Panel(child=p1, title="Response")
tab2 = Panel(child=p2, title="Treatment")
tab3 = Panel(child=p3, title="Cluster Annotations")

tabs = Tabs(tabs=[ tab1, tab2, tab3 ])



selects = column(desc, checkbox_group,  radio_button_group, radio_button_group_RS, radio_button_group_Shf,  theme_select, button, width=520) # was 420




layout=row(selects,p, tabs)#,selected_points)#create_figure())






#doc = curdoc()
curdoc().theme= theme_black


# add to document
curdoc().add_root(layout)

curdoc().title = "MIsTic: Image tSNE viewer"



# cd image_tSNE_code/bokeh_GUI/bokeh-branch-2.3/examples/app
# bokeh serve --port 5098 --show image_tSNE_GUI

# whole GUI in index.html
# indi canvas in theme.yaml
