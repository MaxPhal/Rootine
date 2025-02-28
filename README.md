# Code of Rootine v 2.0 

Author: Maxime Phalempin and Steffen Schlüter
Email: maxime.phalempin@ufz.de or steffen.schlueter@ufz.de
Release Date: 2020/20/01

The script was used in the manuscript:
```` 
Phalempin, M., Lippold, E., Vetterlein, D., & Schlüter, S. (2021). An improved method for the segmentation of roots from X-ray computed tomography 3D images: Rootine v. 2. Plant Methods, 17, 1-19.
````
If you used Rootine for your root segmentation, please make sure to cite the above reference. 

#1. General description:
This is a shell script that implements the Rootine v.2 workflow on a Linux OS (more information for Windows users at the end of this script).
Two macro files are needed to run this script ('Rootinev2_macro_part1.ijm' and 'Rootinev2_macro_part3.ijm'). 
The first part of the macro performs the stitching, the greyscale value drift correction and prepares the data for filtering with the 3D NLM. 
The filtering of the image is performed from the command line in this executable file.  
The second part of the macro performs the rest of the operations (i.e. from pot wall detection to the analysis of the segmented root system).

#2. Installation:
This script requires Fiji/ImageJ (https://imagej.net/Fiji/Downloads, the version used in this study was ImageJ 1.53c)  and the UNLM filter based on ITK (https://www.nitrc.org/projects/unlmeans) - consult website for installation
The ImageJ plugins 'MorpholibJ' and '3D Hysteresis Thresholding' can be installed by adding the update site IJPB-plugins and 3D ImageJ Suite in the Fiji Menu Help / Update ... / Manage Update Site.
The ImageJ plugin 'Attenuation correction' can be retrieved from this website (http://imagejdocu.tudor.lu/doku.php?id=plugin:stacks:attenuation_correction:start) and installed manually. 

#3. Directories description:
The Rootinev2_macro_part1.ijm and Rootinev2_macro_part3.ijm need to be moved into the PATH_TO_MACRO folder.
The folders that contain the executables for FIJI and UNLM need to be set in PATH_TO_FIJI and PATH_TO_UNLM respectively.

#4. Image file and input description:
This script is tailored for two .mhd files called imagename_top.mhd and imagename_bot.mhd that are concatenated into one image prior to the greyscale value drift correction.  
Make sure that the image size in the x-y dimensions has even numbers. Otherwise there are issues with rescaling of the image.
If you want to start with one image or a different file format (e.g. stack of .tiff files), adapt the content of Rootinev2_macro_part1.ijm. The macro parts requires three input text files containing values described down below. These text files should be named :

"imagename_vol_info.txt"
"imagename_ROI_param.txt"
"imagename_param.txt"

These input files need to be present in the "PATH_TO_IMAGE" folder. Note that all parameters and their values should always be separated with a tabulation.
IMPORTANT : the txt file reading function is sensitive to empty characters. Avoid the presence of empty lines!

"imagename_vol_info.txt" : contains some volume information required to open, concatenate the top and bottom images and analyse RSA.
#INPUT 
nr_layers            		- the number of layers into which the input image will divided. /!\ See important note below /!\ 
overlap_slice_bot_stack	- the number of the slice where the overlap region starts in the bottom image
overlap_slice_top_stack 	- the number of the slice where the overlap region starts in the top image
analysed_depth1_bot		- the number of the slice corresponding to the bottom of the "depth1" layer
analysed_depth1_top		- the number of the slice corresponding to the top of the "depth1" layer
analysed_depth2_bot		- the number of the slice corresponding to the bottom of the "depth2" layer
analysed_depth2_top		- the number of the slice corresponding to the top of the "depth2" layer

# /!\ IMPORTANT NOTE /!\ 
Currently, the memory limitation of the "Tubeness Plugin" used in the "Rootinev2_macro_part3" is set to > 2GB. An error is thrown if the image size exceeds 2GB. It falls down on the user to determine an approriate number of layers into which the full image will be divided. The appropriate value is found by making sure 
that ((imagename_top.raw + imagename_bottom.raw)/nr_layers) <  2 GB. For the bottom and top test images that we provide (size = approx. 7.6 GB, after removing the overlapping areas), the number of layers is set to 4. The "nr_layers" parameter will also be used to split the images prior to 3D NLM. 

"imagename_ROI_param.txt" : contains the coordinates of the bounded ROI
#INPUT 
x1						- the x-coordinate of the bounded ROI in the first slice of the stack
y1					 	- the y-coordinate of the bounded ROI in the first slice of the stack
w1						- the width of the bounded ROI in the first slice of the stack
h1						- the height of the bounded ROI in the first slice of the stack
x2						- the x-coordinate of the bounded ROI in the middle slice of the stack
x2					 	- the y-coordinate of the bounded ROI in the middle slice of the stack
w2						- the width of the bounded ROI in the first middle of the stack
h2						- the height of the bounded ROI in the first middle of the stack
x3						- the x-coordinate of the bounded ROI in the last slice of the stack
x3					 	- the y-coordinate of the bounded ROI in the last slice of the stack
w3						- the width of the bounded ROI in the first last of the stack
h3						- the height of the bounded ROI in the first last of the stack

"imagename_param.txt" : contains the tunable parameters required to perform the root segmentation
INPUT 
blur_radius      		-  the standard deviation of the blur radius of the Gaussian filter kernel for "Unsharp masking"
mask_weight 	     	-  the weight of the mask for "Unsharp masking"
f_r			        -  the root greyvalue factor determined on a test image to retrieve the root greyvalue
R_r				    -  the root greyvalue range which sets the root grey value window during the background removal
dr_min			-  the diameter of the finest root in the image (expressed in number of pixels)
dr_max			-  the diameter of the biggest root in the image (expressed in number of pixels)
median_kernel_size    	-  the kernel size in XYZ dimension of the 3D Median filter
t_v				-  the vesselness threshold to be used during the false negatives recovery step
t_s				-  the size threshold to be used during the false negatives recovery step

Note that the contrast threshold is not to be defined in the param.txt but in this bash file. Please enter the value suited to the noise level in the image in line 113.

#5. Root system architecture analysis description. 
Note that this macro is tailored to analyse two depth layers for root length density (RLD) and root diameter distribution (RDD), as is described in the manuscript. Here, depth 1 refers to the deepest layer in the image stack and depth2 to the most upper layer.  If you want to use these functionality to investigate RSA traits, please enter the corresponding slice numbers of the bottom and top layer for each depth in the input file "vol_info.txt".

# -- Directories --
#Paths to important apps and macros - please adapt:
PATH_TO_FIJI=/home/phalempi/Desktop/programs/Fiji.app
PATH_TO_UNLM=/home/phalempi/Desktop/programs/UnbiasedNonLocalMeans/bin/Linux
PATH_TO_MACRO=/media/phalempi/macro
PATH_TO_IMAGE=/media/phalempi/images
IMAGE=07

# -- Workflow (LINUX OS) --
### Part 1: Image stichting, normalization and preparation for filtering 
REQUIRES - vol_info.txt
```` shell 
$PATH_TO_FIJI/ImageJ-linux64 -macro $PATH_TO_MACRO/Rootinev2_macro_part1.ijm "$PATH_TO_IMAGE/$IMAGE"
````

### Part 2 : Noise removal with a 3D Non-local means filter
# Filtering the image layers as defined by the nr_layers (see Important note above).

```` shell 
t_c=60 
for i in 01 02 03 04 # please adapt this sequence according to the number of layers defined. 
do
$PATH_TO_UNLM/UnbiasedNonLocalMeans  --hp 1.0 --rc 1,1,1 --rs 3,3,3 --sigma $t_c $PATH_TO_IMAGE/$IMAGE"_norm_"$i"".mhd $PATH_TO_IMAGE/$IMAGE"_norm_"$i"_nlm".mhd
done
````
### Part 3 : Preprocessing (continued), root segmentation and analysis 
# REQUIRES - vol_info.txt, ROI_param.txt and param.txt

```` shell 
$PATH_TO_FIJI/ImageJ-linux64 -macro $PATH_TO_MACRO/Rootinev2_macro_part3.ijm "$PATH_TO_IMAGE/$IMAGE"
````

# -- Workflow (Windows) --

We haven't tested the 3D non-local means filter for Windows yet. But there is a version of the 3D non-local means filter for Windows OS on the website (https://www.nitrc.org/projects/unlmeans). The ImageJ macros in part 1 and part 3 can be carried out in a Windows Powershell with a few adaptions:
```` shell
SET PATH_TO_FIJI=C:\Users\phalempi\Desktop\Image_Analysis\Fiji.app
SET PATH_TO_MACRO=C:\Users\phalempi\Desktop\Image_Analysis\macro
SET PATH_TO_IMAGE=C:\Users\phalempi\Desktop\Image_Analysis\image
SET IMAGE=imagename
````
# 1. STEP: Image concatenation and attenuation correction
```` shell
%PATH_TO_FIJI%\ImageJ-win64.exe -batch %$PATH_TO_MACRO/Rootinev2_macro_part1.ijm "$PATH_TO_IMAGE/$IMAGE"
````
# 3. STEP: Main Segmentation workflow
```` shell
# %PATH_TO_FIJI%\ImageJ-win64.exe -batch %$PATH_TO_MACRO/Rootinev2_macro_part3.ijm "$PATH_TO_IMAGE/$IMAGE"
````
