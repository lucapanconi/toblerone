# Topological Boundary Line Estimation using Recurrence of Neighbouring Emissions
Topological Boundary Line Estimation using Recurrence of Neighbouring Emissions, or TOBLERONE, is a topological data analysis algorithm which uses persistent homology to identify objects of arbitrary shape in grayscale images. The programme is run entirely through RStudio and requires only the [**tiff**](https://cran.r-project.org/web/packages/tiff/index.html) package to function.

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
In the example code below, we implement the TOBLERONE library, upload an example image, define the image name and specify an output location for results. The **readTIFF** function is included with the **tiff** package and comes pre-installed with TOBLERONE.
```{r}
#Implement TOBLERONE library.
library(toblerone)

#Upload example image.
example_image <- readTIFF("C:/Users/user_name/Desktop/image_folder/example_image.tif")

#Define example image name.
example_image_name <- "example_image"

#Define output location for results.
output_location <- "C:/Users/user_name/Desktop/output_folder/"

#Run TOBLERONE on the example image.
toblerone(example_image, example_image_name, output_location)
```
From here, follow the instructions given in the command prompt or read the additional information below. Note that the image you upload does not necessarily need to be a .tif file, but must be in matrix form, ideally with intensity values between 0 and 1.

## Algorithm
This section contains information on how TOBLERONE works. You do not need to know this to use the algorithm, but you may find it useful. The algorithm technically only requires one parameter per image, but offers the ability to segment images into sub-images for more specific analysis. At each stage, the important data will be saved so that, if the algorithm is suddenly cut off, it can be picked up from where it left off without needing to wait again. The algorithmic process is as follows:
1. **Filtrations:** Derives the filtrations and coordinates matrices which stores the filtration value of each pixel. First the image undergoes a simple background/foreground segmentation. Then, each pixel is assigned a value from 1 to the total number of pixels, from highest to lowest intensity value (this is known as the filtrations matrix). If the filtrations and coordinates matrices are already available at the output location, they will be uploaded and used instead.
2. **Primary Segmentation:** The first application of persistent homology. This is used to gain a rough idea of where each object is within the image. The user is asked to input a number between 0 and 1, representing the persistence threshold. Once given, the persistent homology representation will be calculated and objects derived. If the user is not satisfied with the segmentation, they can input a different threshold. If the segmentation is sufficient, the user can enter “Y” (in upper or lower case) to accept.
3. **Secondary Segmentation:** For each distinct object found in primary segmentation, a window will be drawn around the object and persistent homology can be repeated with a different threshold. This can be used to split up objects in close proximity or elucidate those not fully identified by the original segmentation. For each object, the user is asked for a threshold and the persistent homology representation (for that window) will be calculated. This can be repeated as many times as necessary. If the secondary segmentation is sufficient, the user can enter “Y” to accept. If the user decides not to include that particular object, they can enter “N”. The window size will be automatically determined by the base window size parameter, but can be altered: the window can be increased in the horizontal direction by entering a number greater than 1 (the floor of the number will be taken to be the window size) and in the vertical direction by entering a number less than or equal to -1 (the floor of the absolute value will be taken to be the window size).
4. **Boundary Identification:** Once secondary segmentation is complete and all objects have been identified, a variant of the [**Swinging Arm method**](http://www.geosensor.net/papers/galton06.GISCIENCE.pdf) will be used to determine the boundary of each individual object. From then, each layer within the loop will be calculated – this is done by observing each of the four adjacent pixels around each pixel in the loop, starting from the pixel which is closest to 90° anticlockwise (since the loop itself is oriented anticlockwise), and counting the number of other loop pixels within that row and column (to determine if it is inside or outside the loop – only interior pixels are added). Ultimately this gives a list of the layered loops within each object.
5. **Layer Selection:** At this stage, the user can enter a number to display the interior loops (1 being the exterior loop, 2 the next loop in, and so on). If the loop is "pinched off" at any point and splits into separate loops, these loops will be displayed together but saved separately. Once an appropriate loop has been found, the user can enter “Y” to accept the loop or “N” to discard. Once the loops of all objects have been chosen, all loops will be saved alongside images in the output folder.

There are several parameters which control the user interface. These include:
* **failsafe** – Boolean value which dictates whether failsafe measures will be used. If TRUE and data is available in the output folder, the user will be given the choice of which stage to start at. If FALSE, then entire process (except deriving filtrations and coordinates if they are available) will be completed.
* **use_image_decomposition** – Boolean value which dictates whether the image will undergo initial decomposition into foreground and background. Recommended to be kept true unless exceptionally dark image is used.
* **check_objects** – Boolean value which dictates whether secondary segmentation will be used. Recommended to be kept true.
* **base_window_size** – An integer value which describes the number of pixels that will be included in both the horizontal and vertical directions on either side of the image. The smaller this is, the more compact the window will be around each segmented image and the faster the computation time.
* **identify_boundary** – Boolean value which dictates whether Swinging Arm method will be performed to identify object boundaries and the layers beneath.
* **smooth_image** – Boolean value which dictates whether identified objects will be smoothed by removing one pixel-wide protrusions during boundary identification. Recommended to be kept true.
* **colourblind_palette** – Boolean value which controls whether the colourblind palette will be used throughout.

## License
Licensed under GNU General Public License v3.0.