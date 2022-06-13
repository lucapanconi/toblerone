#' Calculates the Generalised Polarisation (GP) image of returns the GP values along the line profiles identified by TOBLERONE.
#'
#' @param image_1 A matrix representing the first image (preferably ordered channel).
#' @param image_2 A matrix representing the second image (preferably disordered channel).
#' @param image_name A string representing the name of the image originally used in TOBLERONE.
#' @param output A string denoting the output file location as used in TOBLERONE originally.
#' @return A list of the GP values along the line profiles identified by TOBLERONE.
#' @export
get_gp_values <- function(image_1, image_2, image_name, output){
  #Calculate GP image.
  gp_image <- create_gp_image(image_1, image_2)
  
  #Find all loops.
  setwd(output)
  #Get file names.
  file_names <- list.files(pattern = paste(image_name, "_loop_", sep = ""))
  image_names <- list.files(pattern = paste(image_name, "_loop_image_", sep = ""))
  file_names <- file_names[!file_names %in% image_names]
  
  #Iterate over all loop labels.
  all_gp_values <- list()
  found <- 0
  for(loop_label in file_names){
    #Upload the loop.
    loop <- read.csv(loop_label)
    found <- found + 1
    
    #Create list of GP values.
    gp_values <- c()
    for(pixel in 1:nrow(loop)){
      gp_values <- append(gp_values, gp_image[loop[pixel,1], loop[pixel,2]])
    }
    
    #Add mean GP value to list.
    all_gp_values[[found]] <- gp_values
  }
  #Return all mean GP values.
  return(all_gp_values)
}