# upU-net_Perlin
upU-Net architecture for Perlin noise removal in fluorescence microscopy imaging.

#### Main Program

The main.m file has two uses, depending on the *training* option.

 - *training*=0: once the size of the images of the dataset is chosen (128 or 256), it loads the suitable net (**upUNet128.mat** or **upUNet256.mat**) and restore the images in the datastore, removing the Perlin noise.
 - *training*=1: it allows to choose the architecture for the upU-Net and then train it on one of the two datasets contained in the repository. One can generate new dataset of desired dimension using the file datasetGen.m
 
In case of usage of this code (or part of it), please cite

Benfenati, A. upU-Net Approaches for Background Emission Removal in Fluorescence Microscopy. J. Imaging 2022, 8, 142. https://doi.org/10.3390/jimaging8050142 
 
#### List of files

- **datasetGen.m**: it allows to create the dataset for experiments.synthetic microscopy images affected by Perlin, Gaussian and Poisson noise. The user can choose the option to blur the images via a PSF.
- **main.m**: main file
- **particleSIM.m**: MatLab function that generate synthetic microscopy images affected by Perlin, Gaussian and Poisson noise. The user can choose the option to blur the images via a PSF.
- **partionadataIMAGES.m**: custom function to divide the full dataset in training, validation and test dataset.
- **perlinNoise3D.m**: function to generate Perlin noise.
- **setNetwork.m**: it sets the network, provided the image size, the number of blocks and of the filters of the first block 

