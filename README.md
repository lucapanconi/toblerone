# Topological Boundary Line Estimation using Recurrence of Neighbouring Emissions
Topological Boundary Line Estimation using Recurrence of Neighbouring Emissions, or TOBLERONE, is a topological data analysis algorithm which uses persistent homology to identify objects of arbitrary shape in grayscale images. The programme is run entirely through [**RStudio**](https://www.rstudio.com/products/rstudio/) and requires only the [**tiff**](https://cran.r-project.org/web/packages/tiff/index.html) package to function.

## Installation
The package can be installed using pre-existing functions in the [**devtools**](https://cran.r-project.org/web/packages/devtools/index.html) package. It is highly recommended that you use **RStudio**. If you do not have **devtools** installed, enter the following line into the command prompt:
```{r}
install.packages("devtools")
```
Then, simply enter the line below into the command prompt:
```{r}
devtools::install_github("lucapanconi/toblerone")
```

## Usage
In the example code below, we implement the TOBLERONE library, upload an example image, define the image name and specify an output location for results. This package includes an example dataset called "sj_example", which features images of *Schizosaccharomyces japonicus* cells. This is given by a list with two entries, corresponding to the disordered and ordered images respectively.
```{r}
#Implement TOBLERONE library.
library(toblerone)

#Upload example image.
example_image <- sj_example[[1]]

#Define example image name.
example_image_name <- "example_image"

#Define output location for results (ideally create a new folder for this).
output_location <- "C:/Users/user_name/Desktop/output_folder/"

#Run TOBLERONE on the example image.
toblerone(example_image, example_image_name, output_location)
```

Alternatively, we can upload images using the **readTIFF** function, which is included with the **tiff** package and comes pre-installed with TOBLERONE.

```{r}
#Upload example images.
example_image <- tiff::readTIFF("C:/Users/user_name/Desktop/image_folder/example_image.tif")
```

If this is the first time analysing the image, the filtrations matrix (a pre-cursor to topological analysis) will be calculated - this may take a few minutes, but will be saved in the output location for future use. From here, follow the instructions given in the command prompt or read the additional information below. Note that the image you upload does not necessarily need to be a .tif file, but must be in matrix form, ideally with intensity values between 0 and 1.

The package also contains functions for automatically calculating and extracting Generalised Polarisation (GP) values, which represent the degree of membrane order.
```{r}
#Upload images corresponding to ordered and disordered channels.
disordered_channel <- sj_example[[1]]
ordered_channel <- sj_example[[2]]

#Get all GP values from line profiles (one for each cell) identified by TOBLERONE.
gp_values <- get_gp_values(ordered_channel, disordered_channel, example_image_name, output_location)

#Get mean GP values from line profiles identified by TOBLERONE.
mean_gp_values <- get_mean_gp(ordered_channel, disordered_channel, example_image_name, output_location)
```

## Algorithm
This section contains information on how TOBLERONE works. You do not need to know this to use the algorithm, but you may find it useful. The algorithm technically only requires one parameter per image, but offers the ability to segment images into sub-images for more specific analysis. At each stage, the important data will be saved so that, if the algorithm is suddenly cut off, it can be picked up from where it left off without needing to wait again. The algorithmic process is as follows:
1. **Filtrations:** Derives the filtrations and coordinates matrices which stores the filtration value of each pixel. First the image undergoes a simple background/foreground segmentation. Then, each pixel is assigned a value from 1 to the total number of pixels, from highest to lowest intensity value (this is known as the filtrations matrix). If the filtrations and coordinates matrices are already available at the output location, they will be uploaded and used instead.
2. **Primary Segmentation:** The first application of persistent homology. This is used to gain a rough idea of where each object is within the image. The user is asked to input a number between 0 and 1, representing the persistence threshold. Once given, the persistent homology representation will be calculated and objects derived. If the user is not satisfied with the segmentation, they can input a different threshold. If the segmentation is sufficient, the user can enter “Y” (in upper or lower case) to accept. For the *S. japonicus* example dataset, a persistence threshold of 0.68 was found to give a good segmentation.
3. **Secondary Segmentation:** For each distinct object found in primary segmentation, a window will be drawn around the object and persistent homology can be repeated with a different threshold. This can be used to split up objects in close proximity or elucidate those not fully identified by the original segmentation. For each object, the user is asked for a threshold and the persistent homology representation (for that window) will be calculated. This can be repeated as many times as necessary. If the secondary segmentation is sufficient, the user can enter “Y” to accept. If the user decides not to include that particular object, they can enter “N”. The window size will be automatically determined by the base window size parameter, but can be altered: the window can be increased in the horizontal direction by entering a number greater than 1 (the floor of the number will be taken to be the window size) and in the vertical direction by entering a number less than or equal to -1 (the floor of the absolute value will be taken to be the window size).
4. **Boundary Identification:** Once secondary segmentation is complete and all objects have been identified, a variant of the [**Swinging Arm method**](http://www.geosensor.net/papers/galton06.GISCIENCE.pdf) will be used to determine the boundary of each individual object. From then, each layer within the loop will be calculated – this is done by observing each of the four adjacent pixels around each pixel in the loop, starting from the pixel which is closest to 90° anticlockwise (since the loop itself is oriented anticlockwise), and counting the number of other loop pixels within that row and column (to determine if it is inside or outside the loop – only interior pixels are added). Ultimately this gives a list of the layered loops within each object. This process may take a few minutes, so ensure you are satisfied with your segmentation before proceeding.
5. **Layer Selection:** At this stage, the user can enter a number to display the interior loops (1 being the exterior loop, 2 the next loop in, and so on). If the loop is "pinched off" at any point and splits into separate loops, these loops will be displayed together but saved separately. Once an appropriate loop has been found, the user can enter “Y” to accept the loop or “N” to discard. Once the loops of all objects have been chosen, all loops will be saved alongside images in the output folder.

There are several parameters which control the user interface. These include:
* **failsafe** – Boolean value which dictates whether failsafe measures will be used. If TRUE and data is available in the output folder, the user will be given the choice of which stage to start at. If FALSE, then entire process (except deriving filtrations and coordinates if they are available) will be completed.
* **use_image_decomposition** – Boolean value which dictates whether the image will undergo initial decomposition into foreground and background. Recommended to be kept true unless exceptionally dark image is used.
* **check_objects** – Boolean value which dictates whether secondary segmentation will be used. Recommended to be kept true.
* **base_window_size** – An integer value which describes the number of pixels that will be included in both the horizontal and vertical directions on either side of the image. The smaller this is, the more compact the window will be around each segmented image and the faster the computation time.
* **identify_boundary** – Boolean value which dictates whether Swinging Arm method will be performed to identify object boundaries and the layers beneath.
* **smooth_image** – Boolean value which dictates whether identified objects will be smoothed by removing one pixel-wide protrusions during boundary identification. Recommended to be kept true.
* **colourblind_palette** – Boolean value which controls whether the colourblind palette will be used throughout.
* **brightness_scale** – Amplifies the brightness by a given non-negative value. Can be used to brighten especially dark images. Defaults to 1.

## Outputs
The output folder will be populated with a series of different data files. Do not remove or delete these, as TOBLERONE will automatically use them if they are available (to save unnecessary recalculations). The outputs of the algorithm will be saved under:
* **loop_image** files - These are tif images displaying the whole image with a single loop highlighted in each.
* **loop** files - These are csv files in which each row corresponds to the coordinates of a pixel in the loop. All pixels are in order of anti-clockwise orientation.

Additionally, an **input_parameters** file will be given, which records each input made by the user during the algorithm. The results of the segmentation and boundary estimation can be displayed using any of the functions below:
```{r}
#Display primary segmentation.
display_segmentation(output_location, example_image_name, stage = 1)

#Display secondary segmentation.
display_segmentation(output_location, example_image_name, stage = 2)

#Display loop boundaries with a green border.
display_boundaries(output_location, example_image_name, colour = c(0, 255, 0))
```

## 3D TOBLERONE
A version of TOBLERONE for 3D stacks and videos is now available. Stacks must be input as arrays with no more than three dimensions. Only primary segmentation is currently possible. The example below highlights an implementation:
```{r}
#Upload individual images from a given location to create a stack.
location <- "C:/Users/user_name/Desktop/stack_folder/"
setwd(location)
files <- lapply(list.files(pattern = ".tif"), function(i){
  return(tiff::readTIFF(i))
})
m <- nrow(files[[1]])
n <- ncol(files[[1]])
frames <- length(files)
stack <- array(0, dim = c(m, n, frames))
for(i in 1:frames){
  stack[,,i] <- files[[i]]
}

#Give stack a name.
stack_name <- "example_stack"

#Define output location.
output <- "C:/Users/user_name/Desktop/output_folder/"

#Run 3D TOBLERONE.
toblerone3D(stack, stack_name, output)
```
From here, the algorithm will function in the same way as 2D TOBLERONE.

## License
Licensed under GNU General Public License v3.0.