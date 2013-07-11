setuvo - "see tumors in vivo" 
microCT subcutaneous tumor segmentation software

Written by Rehan Ali [1] and Cigdem Demir-Gunduz [2] 
[1] Postdoctoral Scholar, Graves Lab, Stanford University Dept of Radiation Oncology
[2] Professor in Computer Science, Bilkent University, Turkey

Assistance by Tuende Szilagyi, Wolfson Medical Vision Labs, University of Oxford, UK

INTRODUCTION
============
setuvo is a Matlab application for semi-automatically finding the boundaries in subcutaneous tumors using preclinical micro-Computed Tomography (microCT) image data. 

INSTRUCTIONS
============

* COMPILATION OF MEX CODE

setuvo uses a MEX-based level set to improve performance. Before you can begin, this needs to be compiled.

If you're using a Mac, the level set is already compiled for you ('levelset3DC.mexmaci64'). You can start using setuvo straight away.

To compile the MEX code, you'll need the Matlab Compiler. 
1. Open Matlab. 
2. Go to the MEX directory in the setuvo folder.
3. Run the "run_me.m" script. This will compile the code and test that it works with a sample dataset.


* PERFORMING A SEGMENTATION

The following instructions will show you how to use setuvo, using a sample dataset that we've included.
1. Open Matlab.

2. Go to the Main directory in the setuvo folder.

3. Run the "segment_ct_tumor.m" script to start the setuvo GUI (graphical user interface).

4. In the "Path to CT image data" field, enter the following
	data/c3_M2R_CT_mouse.dcm

5. Click the "Load Image Data" button. The CT data of a mouse should show up in the main window, in an axial view.

6. You can browse through the different slices using the slider at the bottom of the image viewer.

7. Locate the subcutaneous tumors. There's one on each shoulder of this mouse. Decide which you want to segment.

8. You need to crop the image in three dimensions. We'll start by cropping in the z-direction (in and out of the image data). Find the slice just anterior to the tumor, and click "First Frame". Then find the slice just posterior to the tumor, and click "Last Frame". 

For the left tumor (a small tumor), you can use the z-slices:
	140	205

For the right tumor (which is larger and easier to see), you can use the z-slices:
	122	190

9. Now pick a slice somewhere in the middle of the tumor, and click "Draw ROI". Draw a box that comfortably covers the whole tumor. This will crop the image.

10. Click on "Select Seed Point". RIGHT-click somewhere inside the tumor. Somewhere in the middle will work best. This will initialise the algorithm (i.e. tell it where the tumor is).

11. In the Output Filename Prefix, enter
	test1
This will save the output in the results folder, and give it a prefix 'test1'. You can run the level set several times, each time with a different prefix, if you want to test the effects of changing the algorithm parameters.

12. Click "Segment". This will start the algorithm. You'll see some text appear in the main Matlab window. 

The algorithm should take a few seconds, after which you'll see the main view window update and scroll through the image data. The segmentation result will be superimposed in red, and the initialisation region will be highlighted in yellow/green.


* COMPARING AGAINST MANUAL SEGMENTATIONS

You can compare the results of the segmentation against a manually defined contour. To do this, you'll use the "Analyse Results" box at the bottom.

1. In the "Saved Filename Prefix", enter the prefix that you used earlier (in this case, test1).

2. In the "Path to Manual Seg data", enter the following:
	data/c3_M2R_CT_mouse.dcm-labels.mat

3. Click the "Compare Result vs Manual Seg" button. 

You'll see the main view screen update and scroll through the image set. The manual segmentation will be shown in green, and the setuvo segmentation will be in red. Areas of overlap will be in yellow.


* USING THE LEVEL SET OUTPUTS

The level set saves several files into the results folder. These include the following:

[prefix]_3Dfluthru_vs_manual.avi
	Video of manual vs setuvo results.

[prefix]_3Dflythru.avi
	Video of setuvo results.

[prefix]_convergence_plot.png
	Plot of the total number of voxels segmentated at each iteration of the level set.

[prefix]_dice_vs_t.png
	Plot of how the dice score (degree of overlap between manual and setuvo results) varies with level set iteration number.

[prefix]_fullsize_seg.bin
	Binary 3D volume of the same resolution as the original uncropped CT dataset, containing the segmentation result.

[prefix]_tmap.bin
	Integer-valued 3D volume of the same resolution as the original uncropped CT dataset, containing the levle set results at each iteration number. (Each value i represents the level set output at that iteration number i.)

[prefix].mat
	Matlab-readable outputs. Contains a struct with all the parameters used, and the various input, pre-processing and output images. (More details coming soon about this.)


PROBLEMS?
=========
If you have any problems, you can post a bug report on the Git page at 
	https://github.com/bobbyali/setuvo/issues


CHANGELOG
=========

v2.0 
- Incorporated MEX

v2.1 
- Fixed bug when saving outputs - used wrong string
- Added gradient fix for outer boundary points