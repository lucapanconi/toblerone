# Topological Boundary Line Estimation using Recurrence of Neighbouring Emissions
---
An R package for the implementation of TOpological Boundary Line Estimation using Recurrence Of Neighbouring Emissions (or TOBLERONE). This package contains a function for performing topological image analysis on grayscale images to better identify regions corresponding to underlying objects. Although it was initially designed for cell segmentation, TOBLERONE can deconstruct images of any biological structure.
## Installation
---
The package can be installed using pre-existing functions in the **devtools** (add link here) package. It is highly recommended that you use **RStudio** (add link here). Simply enter the line below into the command prompt.
```{r}
install_github("lucapanconi/toblerone")
```

## Usage
---
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

From here, follow the instructions given in the command prompt.

## License
---
Licensed under GNU General Public License v3.0.