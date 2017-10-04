# light-field-super-resolution
This MATLAB code implements the Graph-Based Light Field Super-Resolution framework described in [1][2].

The author of the code is:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Mattia Rossi ([email](rossi.mattia@gmail.com))*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Signal Processing Laboratory 4 (LTS4)*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Ecole Polytechnique Federale de Lausanne (Switzerland)*

If you use or adapt this code in your work (either as a stand-alone tool or as a component of any algorithm), you need to cite the appropriate articles [1][2].

This code is for academic purpose only: not for commercial or industrial activities.

## How to run the code

No installation is required. Just download the repository, copy it inside your MATLAB workspace, and follow the three steps below.

1. The Heidelberg and Stanford light field dataset that are mentioned in our articles [1][2] are provided (by their respective creators) with different file formats and reference systems. Download the two light field dataset and then run the script `importdata.m` to convert all the light fields to the same format and reference system. Please note that the two dataset have to be downloaded from their respective websites (references to the two datasets are provided in [1][2]).

2. Run the script `testdata.m` to generate the low resolution light fields from the data generated at the previous point.

3. Run the script `GBtest.m` to apply our Graph-Based (GB) Super-Resolution algorithm to the low resolution light fields generated at the previous point. The script `GBtest.m` allows the user to play with the parameters of the super-resolution algorithm. Since each low resolution light field is decomposed into sub light fields that are super-resolved separately and merged at the very final step, in GBtest.m it is also possible to activate the parallel reconstruction, that takes advantage of the MATLAB parfor.

## Input light field conventions

In the case the user wants to input its own light field, this MUST be stored in a MATLAB 2D cell array. Each entry of the cell array must be a light field view, in particular, the view must be an RGB image stored as a `height * width * 3` MATLAB matrix of type uint8.
Regarding the view order, the following convention is adopted:
- moving along a row of the cell array from left to right must be equivalent to move the camera horizontally from left to right in the 3D scene,
- moving along a column of the cell array from the top to the bottom must be equivalent to move the camera vertically from the top to the bottom in 3D the scene.

## A short note

In our code, the `(s,t)` angular reference system has its origin in the top left corner of the cell array described in the *Input light field conventions* section, with the horizontal angular axis pointing to the right and the vertical angular axis pointing to the bottom. However, in our code, the horizontal angular axis is named `s` and the vertical angular axis is named `t`, while in our articles [1][2] their names are swapped. The same holds true for the spatial axes `x` and `y`. It is ONLY a different naming of the axis, and NOT a different arrangement of the light field views.

## Reference
[1] Mattia Rossi and Pascal Frossard, *Graph-Based Light Field Super-Resolution*, in 19th IEEE International Workshop on Multimedia Signal Processing, 2017.

[2] Mattia Rossi and Pascal Frossard, *Light Field Super-Resolution Via Graph-Based Regularization*, in preparation, CoRR, 2017, [Available online]: http://arxiv.org/abs/1701.02141.
