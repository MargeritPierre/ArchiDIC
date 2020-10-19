# ArchiDIC: Two-scale Digital Image Correlation
*A suite of MATLAB scripts dedicated to the kinematical characterization of architectured structures*

### What do these scripts do ?
- `Meshing` the structure from a reference image
- Mesh-based `Global DIC` for the identification of *microscopic* kinematics
- Feature-based `Local DIC` for the quantification of *macroscopic* kinematics

### How to test the proposed procedures ?
1. `copy` the repository content in a folder
2. make this folder the MATLAB's `Current Folder`
3. run `main.m`



## Meshing the structure using a reference image
This procedure is intended to mesh complex geometries using a reference picture that is used to extract the object's boundaries.
It makes use of the [DistMesh](http://persson.berkeley.edu/distmesh/) procedure by Per-Olof Persson to create high-quality triangulations of the geometry.

### Extraction of the `Region Of Interest`
Starting from an image `I` of the geometry, the Region Of Interest (or mesh interior domain) consists in a *binary mask* with 1's *inside* the geometry and 0's *outside*.
It is created in three steps:
1. Substraction of the image background gray level `b = I[refPt]` so that `I = abs(I-b)`
2. Initialize the mast `M` by gray level thresholding using `M = I>t`
3. Order-statistic filtering of the mask with `M = `[`ordfilt2`](https://fr.mathworks.com/help/images/ref/ordfilt2.html)`(M,ord,K)`, which is equivalent to `M = conv2(M,K)>=ord`. The binary filter kernel `K` is a disk of radius `R`
The parameters of the procedure are then the reference frame `I`, the reference background point `refPt`, the threshold value `t`, the statistical order `ord` as well as the kernel radius `R`. 

### Computation of the `Signed Distance Function`

### Meshing


## Digital Image Correlation: generalities


## Global DIC


## Local DIC
