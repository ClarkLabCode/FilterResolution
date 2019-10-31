# Extracting fast receptive fields from slow imaging frame rates

This code demonstrates how fast neural receptive fields can be extracted from multi-photon microscopy with slow frame rates. Details can be found in our [paper](https://doi.org/10.1038/s41467-019-12974-0). If you find this useful, please consider citing: 
```
@Article{Mano2019,
author={Mano, Omer and Creamer, Matthew S. and Matulis, Catherine A. and Salazar-Gatzimas, Emilio and Chen, Juyue and Zavatone-Veth, Jacob A. and Clark, Damon A.},
title={Using slow frame rate imaging to extract fast receptive fields},
journal={Nature Communications},
year={2019},
volume={10},
number={1},
pages={4979},
issn={2041-1723},
doi={10.1038/s41467-019-12974-0},
url={https://doi.org/10.1038/s41467-019-12974-0}
}
```

To use this code, click the green "Clone or download" button at the top of this page and choose "Download ZIP". Exctract the downloaded ZIP file, saving it somewhere convenient on your computer. To reproduce figures from our paper, open ```figureX_panels.m``` in MATLAB (where X is the figure number you'd like to generate). Run the script and choose "Change Folder" when prompted.

To plot responses using different alignment methods in figure 4, change the value of ```choosePlot``` on line ```8```. Setting ```choosePlot``` to ```1``` will use data acquired through a linescan protocol, while all other values will use data acquired as 2d images. Besides the linescan protocol, all other methods use the same underlying data, but align stimulus and response in different ways. We show that even with slow frame rates, choosing the correct alignment allows us to recover the fast neural receptive field.

To exctract voxel timing filters from your own data, you can use the function ```calculateVoxelTimingFilter```. You can run ```help calculateVoxelTimingFilter``` to see documentation for this function (also available as a comment at the top of the file), and you can see a demo of how to use the function in ```demoVoxelTiming.m```
