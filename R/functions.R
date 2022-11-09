#' Perform TOBLERONE on a given image.
#'
#' This function performs TOBLERONE on a given grayscale image. It is recommended
#' that the image is represented in matrix form with all values between 0 and 1.
#'
#' @param image A matrix representing the image to be analysed.
#' @param image_name A string representing the name of the image, of which all files will be saved under.
#' @param output A string denoting the output file location. All results and calculations will be output here.
#' @return A set of images and "loop" csv files corresponding to the pixel coordinates of the boundaries found by TOBLERONE.
#' @export
#Perform TOBLERONE on a given image.
toblerone <- function(image, image_name, output, failsafe = TRUE, use_image_decomposition = TRUE, check_objects = TRUE,
                       base_window_size = 1, identify_boundary = TRUE, smooth_image = TRUE, colourblind_palette = FALSE,
                       brightness_scale = 1){
  #Define empty inputs list.
  inputs <- data.frame()
  
  #Define pre-set RGB colours.
  if(colourblind_palette){
    colour_list <- (rbind(
      c(51, 34, 136),
      c(17, 19, 51),
      c(68, 170, 153),
      c(136, 204, 238),
      c(221, 204, 119),
      c(204, 102, 119),
      c(170, 68, 153),
      c(136, 34, 85)
    ) / 255)
  } else {
    colour_list <- rbind(
      c(255, 153, 153),
      c(255, 204, 153),
      c(255, 255, 153),
      c(153, 255, 153),
      c(153, 255, 255),
      c(153, 153, 255),
      c(204, 153, 255),
      c(255, 153, 204)
    ) / 255
  }
  adjacent_colour <- c(2,3,4,5,6,7,8,1)
  
  #Define offsets used throughout.
  offsets <- cbind(c(-1, -1, 0, 1, 1, 1, 0, -1), c(0, -1, -1, -1, 0, 1, 1, 1))
  
  #Generate filtrations matrix and coordinates data frame.
  #Upload image, which represents the density function, f. Adjust brightness and normalise.
  m <- nrow(image)
  n <- ncol(image)
  image <- image * brightness_scale
  image <- image + (image > 1) * 1 * (matrix(1, nrow = m, ncol = n) - image)
  f <- image/max(image)
  grid::grid.newpage()
  grid::grid.raster(f)
  
  #Set working directory.
  setwd(output)
  
  #Determine whether the filtrations matrix and coordinates data frame have already been saved at the given location.
  #If so, use them, otherwise calculate them.
  #Extract all files associated with the image name.
  file_names <- list.files(pattern = image_name)
  if(paste(image_name, "_filtrations.csv", sep = "") %in% file_names && paste(image_name, "_coordinates.csv", sep = "") %in% file_names){
    filtrations <- as.matrix(read.csv(paste(output, image_name, "_filtrations.csv", sep = "")))
    coordinates <- as.data.frame(read.csv(paste(output, image_name, "_coordinates.csv", sep = "")))
    old_filtrations <- filtrations
    old_coordinates <- coordinates
  } else {
    #Give output message.
    cat("Calculating filtrations...")
    #Determine count of each intensity value.
    scaled_f <- floor(f * 255)
    count <- do.call(c, lapply(0:255, function(i){
      return(sum(scaled_f == i))
    }))
    #Perform discretised image decomposition.
    if(use_image_decomposition){
      #Calculate break-off intensity.
      mean_count <- mean(count)
      for(i in 1:(length(count) - 1)){
        if(count[i] > mean_count && count[i + 1] < mean_count){
          break_off <- i
          break
        }
      }
    } else {
      #Set break-off to 0.
      break_off <- 0
    }
    #Record coordinates.
    active <- which(scaled_f + 1 > break_off)
    coordinates <- data.frame(cbind(active %% m + (active %% m == 0) * m, (active - 1) %/% m + 1, scaled_f[active]))
    colnames(coordinates) <- c("i", "j", "threshold_number")
    coordinates <- coordinates[order(coordinates$threshold_number, decreasing = TRUE),]
    coordinates <- coordinates[,1:2]
    
    #Create filtrations matrix.
    filtrations <- matrix(Inf, nrow = m, ncol = n)
    filtrations[as.matrix(coordinates)] <- 1:nrow(coordinates)
    
    #The entry (i,j) of filtrations now corresponds to the filtration value of pixel (i,j). Save both the filtration matrix and
    #coordinates data frame for future use.
    write.csv(filtrations, paste(output, image_name, "_filtrations.csv", sep = ""), row.names = FALSE)
    write.csv(coordinates, paste(output, image_name, "_coordinates.csv", sep = ""), row.names = FALSE)
    old_filtrations <- filtrations
    old_coordinates <- coordinates
  }
  
  #If failsafe mode is active, check for saved data to skip unnecessary steps.
  checks <- rep(0, 3)
  if(failsafe){
    #Checks correspond to:
    #1 - Initial segmentation. First images list.
    #2 - Secondary segmentation. Second images list.
    #3 - Loop identification. Layer loops list.
    #4 - Layer selection.
    allowed_stages <- c(1)
    if(paste(image_name, "_primary_segmentation.RData", sep = "") %in% file_names){
      primary_segmentation <- readRDS(paste(output, image_name, "_primary_segmentation.RData", sep = ""))
      checks[1] <- 1
      allowed_stages <- append(allowed_stages, 2)
    }
    if(paste(image_name, "_secondary_segmentation.RData", sep = "") %in% file_names){
      secondary_segmentation <- readRDS(paste(output, image_name, "_secondary_segmentation.RData", sep = ""))
      windows <- readRDS(paste(output, image_name, "_analysis_windows.RData", sep = ""))
      checks[2] <- 1
      allowed_stages <- append(allowed_stages, 3)
    }
    if(paste(image_name, "_layer_identification.RData", sep = "") %in% file_names){
      layer_identification <- readRDS(paste(output, image_name, "_layer_identification.RData", sep = ""))
      checks[3] <- 1
      allowed_stages <- append(allowed_stages, 4)
    }
    if(paste(image_name, "_input_parameters.csv", sep = "") %in% file_names){
      input_parameters <- read.csv(paste(output, image_name, "_input_parameters.csv", sep = ""))
    }
    if(max(checks) == 1){
      cat("Existing data found. Select a number to initiate TOBLERONE from the following stages:\n")
      cat("1: Primary Segmentation.\n")
      if(checks[1] == 1){
        cat("2: Secondary segmentation.\n")
        inputs <- input_parameters[input_parameters[,1] == 1,]
      }
      if(checks[2] == 1){
        cat("3: Boundary identification.\n")
      }
      if(checks[3] == 1){
        cat("4: Layer selection.\n")
        inputs <- rbind(input_parameters[input_parameters[,1] == 1,], input_parameters[input_parameters[,1] == 2,])
      }
      input <- as.numeric(readline("Select stage: "))
      while(TRUE){
        if(input %in% allowed_stages){
          stage <- input
          break
        } else {
          input <- as.numeric(readline("Enter a number given above: "))
        }
      }
    } else {
      stage <- 1
    }
  } else {
    stage <- 1
  }
  
  #Primary Segmentation
  #Repeat until suitable threshold is chosen.
  number_of_pixels <- nrow(coordinates)
  if(stage == 1){
    filtrations <- cbind(Inf, rbind(Inf, filtrations, Inf), Inf)
    cat("\nPrimary Segmentation: Enter a number between 0 and 1 as a persistence threshold.")
    first_time <- TRUE
    while(TRUE){
      #Get persistence threshold.
      if(first_time){
        input <- readline(prompt = "Enter number: ")
        first_time <- FALSE
      } else {
        input <- readline(prompt = "Enter number, or type Y if finished: ")
      }
      inputs <- rbind(inputs, c(1, input))
      if(tolower(input) == "y"){
        break
      } else {
        #Set persistence threshold.
        persistence_threshold <- as.numeric(input)
        
        #Undertake image segmentation via persistent homology.
        #Create arrays: r contains the roots of each pixel, e contains the set number to which each pixel currently belongs and
        #U contains the sets of pixels.
        r <- rep(0, number_of_pixels)
        e <- rep(0, number_of_pixels)
        U <- list()
        
        #Iterate over all pixels starting from the highest density. Perform persistent homology to generate the partition.
        set_number <- 0
        for(i in 1:number_of_pixels){
          #Extract coordinates of i.
          coordinate_1 <- coordinates[i,1]
          coordinate_2 <- coordinates[i,2]
          
          #Let N be the set of neighbours of i with index lower than i.
          neighbour_filtration <- filtrations[cbind(coordinate_1 + offsets[,1] + 1, coordinate_2 + offsets[,2] + 1)]
          N <- neighbour_filtration[neighbour_filtration < i]
          
          #Check length of N.
          if(length(N) == 0){
            #If N is empty, then i is a peak. Create a new subset in U containing just i. Update r and e.
            r[i] <- i
            set_number <- set_number + 1
            e[i] <- set_number
            U[[set_number]] <- c(i)
          }else{
            #Organise the neighbours in N in descending order of the density of their roots.
            R <- r[N]
            neighbours <- as.data.frame(cbind(R, N))
            neighbours <- neighbours[order(neighbours$R),]
            
            #If N is non-empty, there must be at least one neighbour of higher density. Cycle through all neighbours.
            potentials <- unique(neighbours[f[as.matrix(coordinates[neighbours[,1],1:2])] <= f[coordinate_1,coordinate_2] + persistence_threshold,1])
            
            if(length(potentials) == 0){
              #Node does not connect to any neighbours. Create new set for node.
              r[i] <- i
              set_number <- set_number + 1
              e[i] <- set_number
              U[[set_number]] <- c(i)
            }else{
              j <- potentials[1]
              U[[e[j]]] <- append(U[[e[j]]], i)
              r[i] <- j
              e[i] <- e[j]
              if(length(potentials) > 1){
                #Merge the sets of the potentials.
                for(k in potentials[-1]){
                  e_k <- U[[e[k]]]
                  U[[e[j]]] <- append(U[[e[j]]], U[[e[k]]])
                  e[e_k] <- e[j]
                  r[e_k] <- r[j]
                }
              }
            }
          }
        }
        
        #Output only the sets with root density greater than the persistence threshold.
        #If using whole image, add all objects to the same image.
        images <- list()
        
        #Set up colour image.
        whole_colour <- array(f, dim = c(m, n, 3))
        
        #Determine objects.
        found <- 0
        current_colour <- 1
        for(i in unique(r)){
          if(f[coordinates[i,1], coordinates[i,2]] >= persistence_threshold){
            #Initialise new image.
            found <- found + 1
            image <- matrix(0L, nrow = m, ncol = n)
            #Update colours.
            cols <- colour_list[current_colour,]
            current_colour <- adjacent_colour[current_colour]
            #Fill images.
            image[as.matrix(coordinates[U[[e[i]]],])] <- 1
            whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),1)] <- cols[1]
            whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),2)] <- cols[2]
            whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),3)] <- cols[3]
            #Store image.
            images[[found]] <- image
          }
        }
        #Display whole colour image.
        grid::grid.newpage()
        grid::grid.raster(whole_colour)
        
        #Output number of objects found and time taken.
        if(found == 1){
          cat("Found 1 object.")
        } else{
          cat("Found ", found, " objects.")
        }
      }
    }
    primary_segmentation <- images
    primary_segmentation[[length(primary_segmentation) + 1]] <- f
    saveRDS(primary_segmentation, paste(output, image_name, "_primary_segmentation.RData", sep = ""))
    stage <- stage + 1
    #Save inputs.
    colnames(inputs) <- c("stage", "input")
    write.csv(inputs, paste(output, image_name, "_input_parameters.csv", sep = ""), row.names = FALSE)
  } else if(checks[1] == 1){
    images <- primary_segmentation[-length(primary_segmentation)]
  }
  
  #Secondary Segmentation
  #If checking objects, repeat for all found objects.
  if(stage == 2){
    filtrations <- old_filtrations
    coordinates <- old_coordinates
    filtrations <- cbind(Inf, rbind(Inf, filtrations, Inf), Inf)
    windows <- list()
    display_image <- array(f, dim = c(m, n, 3))
    display_colour <- 1
    if(check_objects){
      new_images <- list()
      new_found <- 0
      cat("\nSecondary Segmentation: Enter a number between 0 and 1 as a persistence threshold for the highlighted segment.")
      for(image in images){
        #Create display image.
        cols <- colour_list[display_colour,]
        display_colour <- adjacent_colour[display_colour]
        original_display <- display_image
        display_image_outline <- display_image
        #Create red outline.
        for(i in 2:(m-1)){
          for(j in 2:(n-1)){
            if(image[i,j] == 1){
              # display_image_outline[cbind(i + offsets[,1],j + offsets[,2],1)] <- 1
              # display_image_outline[cbind(i + offsets[,1],j + offsets[,2],1)] <- 0
              # display_image_outline[cbind(i + offsets[,1],j + offsets[,2],1)] <- 0
              for(k in 1:8){
                display_image_outline[i + offsets[k,1],j + offsets[k,2],] <- c(1,0,0)
              }
            }
          }
        }
        #Fill in colour.
        indicator <- image == 1
        display_image[,,1] <- cols[1] * indicator + display_image[,,1] * !indicator
        display_image[,,2] <- cols[2] * indicator + display_image[,,2] * !indicator
        display_image[,,3] <- cols[3] * indicator + display_image[,,3] * !indicator
        display_image_outline[,,1] <- cols[1] * indicator + display_image_outline[,,1] * !indicator
        display_image_outline[,,2] <- cols[2] * indicator + display_image_outline[,,2] * !indicator
        display_image_outline[,,3] <- cols[3] * indicator + display_image_outline[,,3] * !indicator
        
        #Display the image.
        grid::grid.newpage()
        grid::grid.raster(display_image_outline)
        
        #Reset window size and create initial window.
        window_size_i <- base_window_size
        window_size_j <- base_window_size
        
        #Determine i and j boundaries.
        active <- which(image == 1)
        i_values <- active %% m + (active %% m == 0) * m
        j_values <- (active - 1) %/% m + 1
        
        #Determine image boundaries.
        i_min <- max(min(i_values) - window_size_i, 1)
        i_max <- min(max(i_values) + window_size_i, m)
        j_min <- max(min(j_values) - window_size_j, 1)
        j_max <- min(max(j_values) + window_size_j, n)
        
        #Ask whether to keep image.
        input <- readline(prompt = "Enter number, type Y if finished or type N to discard: ")
        inputs <- rbind(inputs, c(2, input))
        if(tolower(input) == "y"){
          #Add image to new images.
          new_found <- new_found + 1
          new_images[[new_found]] <- image
          #Store window.
          windows[[new_found]] <- rbind(c(i_min, i_max), c(j_min, j_max))
        } else if(tolower(input) == "n"){
          #Reset display
          display_image <- original_display
        } else {
          #Reset display
          display_image <- original_display
          
          #Extract the part of the image that corresponds to that object plus a given window.
          new_filtrations <- old_filtrations[i_min:i_max, j_min:j_max]
          new_coordinates <- old_coordinates[old_coordinates[,1] >= i_min & old_coordinates[,1] <= i_max &
                                               old_coordinates[,2] >= j_min & old_coordinates[,2] <= j_max,]
          new_coordinates[,1] <- new_coordinates[,1] - i_min + 1
          new_coordinates[,2] <- new_coordinates[,2] - j_min + 1
          
          #Display current image.
          new_f <- f[i_min:i_max, j_min:j_max]
          new_image <- image[i_min:i_max, j_min:j_max]
          new_display_image <- display_image[i_min:i_max, j_min:j_max,]
          grid::grid.newpage()
          grid::grid.raster(new_display_image)
          
          #Calculate new number of pixels.
          number_of_pixels <- nrow(new_coordinates)
          
          #Update filtrations matrix.
          new_filtrations[as.matrix(new_coordinates[1:number_of_pixels,])] <- 1:number_of_pixels
          
          #Repeat until a suitable threshold is found.
          coordinates <- new_coordinates[,1:2]
          filtrations <- cbind(Inf, rbind(Inf, new_filtrations, Inf), Inf)
          
          #Set persistence threshold.
          persistence_threshold <- as.numeric(input)
          
          #Undertake image segmentation via persistent homology.
          #Create arrays: r contains the roots of each pixel, e contains the set number to which each pixel currently belongs and
          #U contains the sets of pixels.
          r <- rep(0, number_of_pixels)
          e <- rep(0, number_of_pixels)
          U <- list()
          
          #Iterate over all pixels starting from the highest density. Perform persistent homology to generate the partition.
          set_number <- 0
          for(i in 1:number_of_pixels){
            #Extract coordinates of i.
            coordinate_1 <- coordinates[i,1]
            coordinate_2 <- coordinates[i,2]
            
            #Let N be the set of neighbours of i with index lower than i.
            neighbour_filtration <- filtrations[cbind(coordinate_1 + offsets[,1] + 1, coordinate_2 + offsets[,2] + 1)]
            N <- neighbour_filtration[neighbour_filtration < i]
            
            #Check length of N.
            if(length(N) == 0){
              #If N is empty, then i is a peak. Create a new subset in U containing just i. Update r and e.
              r[i] <- i
              set_number <- set_number + 1
              e[i] <- set_number
              U[[set_number]] <- c(i)
            }else{
              #Organise the neighbours in N in descending order of the density of their roots.
              R <- r[N]
              neighbours <- as.data.frame(cbind(R, N))
              neighbours <- neighbours[order(neighbours$R),]
              
              #If N is non-empty, there must be at least one neighbour of higher density. Cycle through all neighbours.
              potentials <- unique(neighbours[new_f[as.matrix(coordinates[neighbours[,1],1:2])] <= new_f[coordinate_1,coordinate_2] + persistence_threshold,1])
              
              if(length(potentials) == 0){
                #Node does not connect to any neighbours. Create new set for node.
                r[i] <- i
                set_number <- set_number + 1
                e[i] <- set_number
                U[[set_number]] <- c(i)
              }else{
                j <- potentials[1]
                U[[e[j]]] <- append(U[[e[j]]], i)
                r[i] <- j
                e[i] <- e[j]
                if(length(potentials) > 1){
                  #Merge the sets of the potentials.
                  for(k in potentials[-1]){
                    e_k <- U[[e[k]]]
                    U[[e[j]]] <- append(U[[e[j]]], U[[e[k]]])
                    e[e_k] <- e[j]
                    r[e_k] <- r[j]
                  }
                }
              }
            }
          }
          
          #Output only the sets with root density greater than the persistence threshold.
          saved_new_images <- list()
          whole_new_image <- matrix(0L, nrow = i_max - i_min + 1, ncol = j_max - j_min + 1)
          
          #Set up colour image.
          whole_colour <- array(f[i_min:i_max, j_min:j_max], dim = c(i_max - i_min + 1, j_max - j_min + 1, 3))
          
          #Set up shifted coordinates.
          shifted_coordinates <- coordinates
          shifted_coordinates[,1] <- shifted_coordinates[,1] + i_min - 1
          shifted_coordinates[,2] <- shifted_coordinates[,2] + j_min - 1
          
          #Determine objects.
          found <- 0
          current_colour <- 1
          for(i in unique(r)){
            if(new_f[coordinates[i,1], coordinates[i,2]] >= persistence_threshold){
              #Initialise new image.
              found <- found + 1
              new_image <- matrix(0L, nrow = i_max - i_min + 1, ncol = j_max - j_min + 1)
              #Update colours.
              cols <- colour_list[current_colour,]
              current_colour <- adjacent_colour[current_colour]
              #Fill images.
              new_image[as.matrix(coordinates[U[[e[i]]],])] <- 1
              whole_new_image[as.matrix(coordinates[U[[e[i]]],])] <- 1
              whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),1)] <- cols[1]
              whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),2)] <- cols[2]
              whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),3)] <- cols[3]
              display_image[cbind(as.matrix(shifted_coordinates[U[[e[i]]],]),1)] <- cols[1]
              display_image[cbind(as.matrix(shifted_coordinates[U[[e[i]]],]),2)] <- cols[2]
              display_image[cbind(as.matrix(shifted_coordinates[U[[e[i]]],]),3)] <- cols[3]
              #Add to saved images list.
              saved_new_images[[found]] <- new_image
            }
          }
          #Display whole colour image.
          grid::grid.newpage()
          grid::grid.raster(whole_colour)
          
          #Repeat until segmentation is accepted.
          cat("\nEnter a number between 0 and 1 as a persistence threshold for the highlighted segment.\nEnter a positive whole number to alter the width of the given window.\nEnter a negative whole number to alter the height of the given window.")
          while(TRUE){
            #Get input.
            input <- readline(prompt = "Enter number, type Y to take individual segments, type W to take the whole image, or type N to discard: ")
            inputs <- rbind(inputs, c(2, input))
            if(tolower(input) == "w"){
              #Bring the image back to its original size.
              clean <- matrix(0L, nrow = m, ncol = n)
              for(i in 1:(i_max - i_min + 1)){
                for(j in 1:(j_max - j_min + 1)){
                  if(whole_new_image[i,j] == 1){
                    clean[i + i_min - 1, j + j_min - 1] <- 1
                  }
                }
              }
              new_found <- new_found + 1
              new_images[[new_found]] <- clean
              
              #Store window.
              windows[[new_found]] <- rbind(c(i_min, i_max), c(j_min, j_max))
              break
            } else if(tolower(input) == "y"){
              #Bring all images back to original size and add to list.
              for(saved_new_image in saved_new_images){
                clean <- matrix(0L, nrow = m, ncol = n)
                for(i in 1:(i_max - i_min + 1)){
                  for(j in 1:(j_max - j_min + 1)){
                    if(saved_new_image[i,j] == 1){
                      clean[i + i_min - 1, j + j_min - 1] <- 1
                    }
                  }
                }
                new_found <- new_found + 1
                new_images[[new_found]] <- clean
                
                #Store window.
                windows[[new_found]] <- rbind(c(i_min, i_max), c(j_min, j_max))
              }
              break
            } else if(tolower(input) == "n"){
              #Reset display
              display_image <- original_display
              break
            } else {
              #Reset display
              display_image <- original_display
              
              #Interpret numeric input.
              input <- as.numeric(input)
              
              #If the input is a whole number greater than equal to 1, change the window size to that value. Otherwise update the
              #persistence threshold and perform TOBLERONE again.
              if(input >= 1){
                #Set window size.
                window_size_j <- floor(input)
                
                #Determine new image boundaries.
                i_min <- max(min(i_values) - window_size_i, 1)
                i_max <- min(max(i_values) + window_size_i, m)
                j_min <- max(min(j_values) - window_size_j, 1)
                j_max <- min(max(j_values) + window_size_j, n)
                
                #Extract the part of the image that corresponds to that object plus a given window.
                new_filtrations <- old_filtrations[i_min:i_max, j_min:j_max]
                new_coordinates <- old_coordinates[old_coordinates[,1] >= i_min & old_coordinates[,1] <= i_max &
                                                     old_coordinates[,2] >= j_min & old_coordinates[,2] <= j_max,]
                new_coordinates[,1] <- new_coordinates[,1] - i_min + 1
                new_coordinates[,2] <- new_coordinates[,2] - j_min + 1
                
                #Display current image.
                new_f <- f[i_min:i_max, j_min:j_max]
                new_image <- image[i_min:i_max, j_min:j_max]
                new_display_image <- display_image[i_min:i_max, j_min:j_max,]
                grid::grid.newpage()
                grid::grid.raster(new_display_image)
                
                #Calculate new number of pixels.
                number_of_pixels <- nrow(new_coordinates)
                
                #Update filtrations matrix.
                new_filtrations[as.matrix(new_coordinates[1:number_of_pixels,])] <- 1:number_of_pixels
                
                #Finalise.
                coordinates <- new_coordinates[,1:2]
                filtrations <- cbind(Inf, rbind(Inf, new_filtrations, Inf), Inf)
              } else if(input <= -1){
                #Set window size.
                window_size_i <- floor(-input)
                
                #Determine new image boundaries.
                i_min <- max(min(i_values) - window_size_i, 1)
                i_max <- min(max(i_values) + window_size_i, m)
                j_min <- max(min(j_values) - window_size_j, 1)
                j_max <- min(max(j_values) + window_size_j, n)
                
                #Determine new image boundaries.
                i_min <- max(min(i_values) - window_size_i, 1)
                i_max <- min(max(i_values) + window_size_i, m)
                j_min <- max(min(j_values) - window_size_j, 1)
                j_max <- min(max(j_values) + window_size_j, n)
                
                #Extract the part of the image that corresponds to that object plus a given window.
                new_filtrations <- old_filtrations[i_min:i_max, j_min:j_max]
                new_coordinates <- old_coordinates[old_coordinates[,1] >= i_min & old_coordinates[,1] <= i_max &
                                                     old_coordinates[,2] >= j_min & old_coordinates[,2] <= j_max,]
                new_coordinates[,1] <- new_coordinates[,1] - i_min + 1
                new_coordinates[,2] <- new_coordinates[,2] - j_min + 1
                
                #Display current image.
                new_f <- f[i_min:i_max, j_min:j_max]
                new_image <- image[i_min:i_max, j_min:j_max]
                new_display_image <- display_image[i_min:i_max, j_min:j_max,]
                grid::grid.newpage()
                grid::grid.raster(new_display_image)
                
                #Calculate new number of pixels.
                number_of_pixels <- nrow(new_coordinates)
                
                #Update filtrations matrix.
                new_filtrations[as.matrix(new_coordinates[1:number_of_pixels,])] <- 1:number_of_pixels
                
                #Finalise.
                coordinates <- new_coordinates[,1:2]
                filtrations <- cbind(Inf, rbind(Inf, new_filtrations, Inf), Inf)
              } else {
                #Set persistence threshold.
                persistence_threshold <- input
                
                #Undertake image segmentation via persistent homology.
                #Create arrays: r contains the roots of each pixel, e contains the set number to which each pixel currently belongs and
                #U contains the sets of pixels.
                r <- rep(0, number_of_pixels)
                e <- rep(0, number_of_pixels)
                U <- list()
                
                #Iterate over all pixels starting from the highest density. Perform persistent homology to generate the partition.
                set_number <- 0
                for(i in 1:number_of_pixels){
                  #Extract coordinates of i.
                  coordinate_1 <- coordinates[i,1]
                  coordinate_2 <- coordinates[i,2]
                  
                  #Let N be the set of neighbours of i with index lower than i.
                  neighbour_filtration <- filtrations[cbind(coordinate_1 + offsets[,1] + 1, coordinate_2 + offsets[,2] + 1)]
                  N <- neighbour_filtration[neighbour_filtration < i]
                  
                  #Check length of N.
                  if(length(N) == 0){
                    #If N is empty, then i is a peak. Create a new subset in U containing just i. Update r and e.
                    r[i] <- i
                    set_number <- set_number + 1
                    e[i] <- set_number
                    U[[set_number]] <- c(i)
                  }else{
                    #Organise the neighbours in N in descending order of the density of their roots.
                    R <- r[N]
                    neighbours <- as.data.frame(cbind(R, N))
                    neighbours <- neighbours[order(neighbours$R),]
                    
                    #If N is non-empty, there must be at least one neighbour of higher density. Cycle through all neighbours.
                    potentials <- unique(neighbours[new_f[as.matrix(coordinates[neighbours[,1],1:2])] <= new_f[coordinate_1,coordinate_2] + persistence_threshold,1])
                    
                    if(length(potentials) == 0){
                      #Node does not connect to any neighbours. Create new set for node.
                      r[i] <- i
                      set_number <- set_number + 1
                      e[i] <- set_number
                      U[[set_number]] <- c(i)
                    }else{
                      j <- potentials[1]
                      U[[e[j]]] <- append(U[[e[j]]], i)
                      r[i] <- j
                      e[i] <- e[j]
                      if(length(potentials) > 1){
                        #Merge the sets of the potentials.
                        for(k in potentials[-1]){
                          e_k <- U[[e[k]]]
                          U[[e[j]]] <- append(U[[e[j]]], U[[e[k]]])
                          e[e_k] <- e[j]
                          r[e_k] <- r[j]
                        }
                      }
                    }
                  }
                }
                
                #Output only the sets with root density greater than the persistence threshold.
                saved_new_images <- list()
                whole_new_image <- matrix(0L, nrow = i_max - i_min + 1, ncol = j_max - j_min + 1)
                
                #Set up colour image.
                whole_colour <- array(f[i_min:i_max, j_min:j_max], dim = c(i_max - i_min + 1, j_max - j_min + 1, 3))
                
                #Set up shifted coordinates.
                shifted_coordinates <- coordinates
                shifted_coordinates[,1] <- shifted_coordinates[,1] + i_min - 1
                shifted_coordinates[,2] <- shifted_coordinates[,2] + j_min - 1
                
                #Determine objects.
                found <- 0
                current_colour <- 1
                for(i in unique(r)){
                  if(new_f[coordinates[i,1], coordinates[i,2]] >= persistence_threshold){
                    #Initialise new image.
                    found <- found + 1
                    new_image <- matrix(0L, nrow = i_max - i_min + 1, ncol = j_max - j_min + 1)
                    #Update colours.
                    cols <- colour_list[current_colour,]
                    current_colour <- adjacent_colour[current_colour]
                    #Fill images.
                    new_image[as.matrix(coordinates[U[[e[i]]],])] <- 1
                    whole_new_image[as.matrix(coordinates[U[[e[i]]],])] <- 1
                    whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),1)] <- cols[1]
                    whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),2)] <- cols[2]
                    whole_colour[cbind(as.matrix(coordinates[U[[e[i]]],]),3)] <- cols[3]
                    display_image[cbind(as.matrix(shifted_coordinates[U[[e[i]]],]),1)] <- cols[1]
                    display_image[cbind(as.matrix(shifted_coordinates[U[[e[i]]],]),2)] <- cols[2]
                    display_image[cbind(as.matrix(shifted_coordinates[U[[e[i]]],]),3)] <- cols[3]
                    #Add to saved images list.
                    saved_new_images[[found]] <- new_image
                  }
                }
                #Display whole colour image.
                grid::grid.newpage()
                grid::grid.raster(whole_colour)
              }
            }
          }
        }
      }
    }
    
    #Smooth images.
    if(smooth_image){
      adjacent_pixels <- cbind(c(0,-1,0,1), c(1,0,-1,0))
      updates <- lapply(1:length(new_images), function(i){
        #Initialise empty lists.
        updated_images <- list()
        updated_windows <- list()
        #Extract image and window.
        image_to_change <- new_images[[i]]
        current_window <- windows[[i]]
        current_i_min <- current_window[1,1] + 1
        current_i_max <- current_window[1,2] + 1
        current_j_min <- current_window[2,1] + 1
        current_j_max <- current_window[2,2] + 1
        #Create filler rows and columns.
        image_to_change <- cbind(0, rbind(0, image_to_change, 0), 0)
        changed <- TRUE
        while(changed){
          changed <- FALSE
          cropped <- image_to_change[current_i_min:current_i_max, current_j_min:current_j_max]
          active <- which(cropped == 1)
          mcrop <- current_i_max - current_i_min + 1
          coords <- cbind(active %% mcrop + (active %% mcrop == 0) * mcrop + current_i_min - 1, (active - 1) %/% mcrop + 1 + current_j_min - 1)
          for(j in 1:nrow(coords)){
            if(sum(image_to_change[as.matrix(cbind(adjacent_pixels[,1] + coords[j, 1], adjacent_pixels[,2] + coords[j, 2]))]) <= 1){
              image_to_change[coords[j, 1], coords[j, 2]] <- 0
              changed <- TRUE
            }
          }
        }
        
        #Check changed image. If it contains more than one connected component, split image.
        #Check through all pixels.
        cropped <- image_to_change[(current_i_min + 1):(current_i_max + 1), (current_j_min + 1):(current_j_max + 1)]
        active <- which(cropped == 1)
        mcrop <- current_i_max - current_i_min + 1
        coords <- cbind(active %% mcrop + (active %% mcrop == 0) * mcrop + current_i_min - 1, (active - 1) %/% mcrop + 1 + current_j_min - 1)
        for(j in 1:nrow(coords)){
          if(image_to_change[coords[j,1], coords[j,2]] == 1){
            #If an active pixel is found, initialise a new image and copy over all connected pixels.
            # updated_images_found <- updated_images_found + 1
            image_to_change[coords[j,1], coords[j,2]] <- 0
            updated_image <- matrix(0L, nrow = m, ncol = n)
            pixels_to_check <- data.frame(coords[j,1], coords[j,2])
            while(nrow(pixels_to_check) > 0){
              #Add current pixel to updated image.
              updated_image[pixels_to_check[1,1] - 1, pixels_to_check[1,2] - 1] <- 1
              for(o in 1:4){
                if(image_to_change[pixels_to_check[1,1] + adjacent_pixels[o,1], pixels_to_check[1,2] + adjacent_pixels[o,2]] == 1){
                  pixels_to_check <- rbind(pixels_to_check, c(pixels_to_check[1,1] + adjacent_pixels[o,1], pixels_to_check[1,2] + adjacent_pixels[o,2]))
                  image_to_change[pixels_to_check[1,1] + adjacent_pixels[o,1], pixels_to_check[1,2] + adjacent_pixels[o,2]] <- 0
                }
              }
              pixels_to_check <- pixels_to_check[-1,]
            }
            #Add updated image to list.
            updated_images[[length(updated_images) + 1]] <- updated_image
            updated_windows[[length(updated_windows) + 1]] <- windows[[i]]
          }
        }
        
        #Return image.
        return(list(updated_images, updated_windows))
      })
      #Extract updated images and windows.
      updated_images <- lapply(updates, function(i){
        return(i[[1]][[1]])
      })
      #Extract updated images and windows.
      updated_windows <- lapply(updates, function(i){
        return(i[[2]][[1]])
      })
    }
    
    #Save secondary segmentation.
    new_images <- updated_images
    windows <- updated_windows
    secondary_segmentation <- new_images
    secondary_segmentation[[length(secondary_segmentation) + 1]] <- f
    saveRDS(secondary_segmentation, paste(output, image_name, "_secondary_segmentation.RData", sep = ""))
    saveRDS(windows, paste(output, image_name, "_analysis_windows.RData", sep = ""))
    stage <- stage + 1
    #Save inputs.
    colnames(inputs) <- c("stage", "input")
    write.csv(inputs, paste(output, image_name, "_input_parameters.csv", sep = ""), row.names = FALSE)
    #Display the image.
    grid::grid.newpage()
    grid::grid.raster(display_image)
  } else if(checks[2] == 1){
    new_images <- secondary_segmentation[-length(secondary_segmentation)]
  }
  
  #Layer Identification
  #Once all objects have been appropriately re-analysed, construct the new images list.
  images <- new_images
  #Identify the boundary of each segmented image found by persistent homology.
  #Determines the boundary of each object, if required.
  if(identify_boundary){
    #Give output message.
    cat("Identifying boundaries...")
    #Determine exterior loops.
    if(stage == 3){
      #Define next squares and adjacent squares.
      next_squares <- c(7,8,1,2,3,4,5,6)
      adjacent_squares <- c(2,3,4,5,6,7,8,1)
      #If there are no images left (all have been deleted), stop.
      images_available <- TRUE
      number_of_images <- length(images)
      if(number_of_images == 0){
        images_available <- FALSE
      }else{
        #Iterate over each image. Scan each pixel left to right, then top to bottom, upon hitting a bright pixel, perform grid
        #sweep and add to loops list.
        loops <- lapply(1:number_of_images, function(index){
          #Extract image and window.
          image <- images[[index]]
          window <- windows[[index]]
          #Construct window for analysis.
          i_min <- window[1,1]
          i_max <- window[1,2]
          j_min <- window[2,1]
          j_max <- window[2,2]
          new_m <- i_max - i_min + 1
          new_n <- j_max - j_min + 1
          image <- image[i_min:i_max, j_min:j_max]
          
          #Find loops.
          active <- which(image == 0)
          clear <- matrix(FALSE, nrow = new_m + 2, ncol = new_n + 2)
          coords <- data.frame(cbind(active %% new_m + (active %% new_m == 0) * new_m, (active - 1) %/% new_m + 1))
          new <- do.call(rbind, lapply(1:nrow(coords), function(i){
            return(cbind(offsets[,1] + coords[i,1], offsets[,2] + coords[i,2]))
          }))
          clear[new + 1] <- TRUE
          active <- which(image == 1 & clear[2:(new_m + 1), 2:(new_n + 1)])
          coords <- cbind(active %% new_m + (active %% new_m == 0) * new_m, (active - 1) %/% new_m + 1)
          clear <- matrix(0, nrow = new_m, ncol = new_n)
          clear[coords] <- 1:nrow(coords)
          
          #Initialise loop node list.
          i <- coords[1,1]
          j <- coords[1,2]
          initial_i <- i
          initial_j <- j
          loop_nodes_x <- c(i)
          loop_nodes_y <- c(j)
          
          #Initialise first start square and current square.
          start_square <- 3
          current_square <- 3
          
          #Perform grid sweep.
          finished <- FALSE
          while(!finished){
            #Get coordinates of current square.
            square_i <- i + offsets[current_square, 1]
            square_j <- j + offsets[current_square, 2]
            #Check if current square is bright.
            if(clear[square_i, square_j] > 0){
              #If already in loop nodes, stop process.
              if(initial_i == square_i && initial_j == square_j){
                finished <- TRUE
              } else{
                #If it isn't, add the current square to the loop nodes and break the loop.
                loop_nodes_x <- append(loop_nodes_x, square_i)
                loop_nodes_y <- append(loop_nodes_y, square_j)
                #Set next start square accordingly.
                current_square <- next_squares[current_square]
                start_square <- next_squares[current_square]
                #Set new current.
                i <- square_i
                j <- square_j
              }
            }else{
              #Increase current square to check next square.
              current_square <- adjacent_squares[current_square]
            }
            #Check if no more points have been found.
            if(current_square == start_square){
              finished <- TRUE
            }
          }
          
          #Return loop.
          return(cbind(loop_nodes_x + i_min - 1, loop_nodes_y + j_min - 1))
        })
      }
    }
    
    #Determine all interior loops.
    #Get all good loops.
    if(stage == 3){
      good_loops <- list()
      found <- 0
      for(i in 1:length(loops)){
        loop <- loops[[i]]
        loop_length <- nrow(loop)
        if((loop[1,1] - loop[loop_length, 1]) ** 2 + (loop[1,2] - loop[loop_length, 2]) ** 2 <= 2){
          found <- found + 1
          good_loops[[found]] <- loop
        }
      }
      
      #Set up layer loops list.
      if(length(good_loops) > 0){
        #For each good loop, create a list in layer loops to store that loop and all loops inside it. For each loop,
        #keep adding all adjacent pixels as their own loop until no more are found.
        layer_loops <- lapply(1:length(images), function(i){
          image <- images[[i]]
          loop <- good_loops[[i]]
          layer_loops <- list(loop)
          check_image <- image == 1
          check_image[loop] <- FALSE
          while(sum(check_image) > 0){
            #Find loops.
            new <- do.call(rbind, lapply(1:nrow(loop), function(i){
              return(cbind(offsets[,1] + loop[i,1], offsets[,2] + loop[i,2]))
            }))
            clear <- matrix(FALSE, nrow = m + 2, ncol = n + 2)
            clear[new + 1] <- TRUE
            active <- which(check_image & clear[2:(m + 1), 2:(n + 1)])
            loop <- cbind(active %% m + (active %% m == 0) * m, (active - 1) %/% m + 1)
            layer_loops[[length(layer_loops) + 1]] <- loop
            check_image[loop] <- FALSE
          }
          return(layer_loops)
        })
      }
      saveRDS(layer_loops, paste(output, image_name, "_layer_identification.RData", sep = ""))
      stage <- stage + 1
    } else if(checks[3] == 1){
      layer_loops <- layer_identification
    }
    
    #Select the loops you want to keep.
    if(stage == 4){
      cat("\nLayer Identification: Enter a positive whole number to view the corresponding layer.")
      found <- length(layer_loops)
      if(found > 0){
        #For each good loop, display the first loop layer and allow the user to input a loop number.
        final_loops <- list()
        loops_found <- 0
        for(i in 1:found){
          current_loop <- layer_loops[[i]][[1]]
          image_colour <- array(f, dim = c(m, n, 3))
          image <- f
          for(j in 1:nrow(current_loop)){
            image_colour[current_loop[j,1], current_loop[j,2],] <- 0
            image_colour[current_loop[j,1], current_loop[j,2], 1] <- 1
            image[current_loop[j,1], current_loop[j,2]] <- 1
          }
          grid::grid.newpage()
          grid::grid.raster(image_colour)
          #Repeat until suitable layer has been found.
          while(TRUE){
            #Take input.
            input <- readline(prompt = "Enter number, type Y if finished or type N to discard: ")
            inputs <- rbind(inputs, c(4, input))
            current_loops <- list(current_loop)
            if(tolower(input) == "y"){
              #Add loop to final loops list.
              for(current_loop in current_loops){
                loops_found <- loops_found + 1
                final_loops[[loops_found]] <- current_loop
              }
              break
            } else if(tolower(input) == "n"){
              break
            } else {
              #Interpret numeric input.
              input <- ceiling(as.numeric(input))
              if(input < 1 || input > length(layer_loops[[i]])){
                cat(paste("Not a valid layer number, enter a number between 1 and ", length(layer_loops[[i]]), ".", sep = ""))
              } else {
                if(is.list(layer_loops[[i]][[input]])){
                  current_loops <- layer_loops[[i]][[input]]
                  current_loop <- data.frame()
                  for(loop in current_loops){
                    current_loop <- rbind(current_loop, loop)
                  }
                  current_loop <- as.matrix(current_loop)
                } else {
                  current_loop <- layer_loops[[i]][[input]]
                  current_loops <- list(current_loop)
                }
                image_colour <- array(f, dim = c(m, n, 3))
                #image <- f#
                for(j in 1:nrow(current_loop)){
                  image_colour[current_loop[j,1], current_loop[j,2],] <- 0
                  image_colour[current_loop[j,1], current_loop[j,2], 1] <- 1
                  #image[current_loop[j,1], current_loop[j,2]] <- 1#
                }
                grid::grid.newpage()
                grid::grid.raster(image_colour)
              }
            }
          }
        }
        
        #Save all loops.
        if(length(final_loops) > 0){
          for(i in 1:length(final_loops)){
            loop <- final_loops[[i]]
            loop_length <- nrow(loop)
            image_colour <- array(f, dim = c(m, n, 3))
            for(k in 1:loop_length){
              image_colour[loop[k,1], loop[k,2],] <- 0
              image_colour[loop[k,1], loop[k,2], 1] <- 1
            }
            tiff::writeTIFF(image_colour, paste(output, image_name, "_loop_image_", i, ".tif", sep = ""))
            write.csv(loop, paste(output, image_name, "_loop_", i, ".csv", sep = ""), row.names = FALSE)
          }
        }
      }
    }
  }
  
  #Save inputs.
  colnames(inputs) <- c("stage", "input")
  write.csv(inputs, paste(output, image_name, "_input_parameters.csv", sep = ""), row.names = FALSE)
}

#' Display a segmentation result from TOBLERONE.
#'
#' This function plots the objects found from either the primary or
#' secondary segmentation.
#'
#' @param image_name A string representing the name of the image.
#' @param output A string denoting the output file location. As defined in TOBLERONE, this is the location of the segmentation outputs.
#' @param stage Set equal to 1 for primary segmentation and 2 for secondary segmentation. Defaults to 1.
#' @return Displays an image of the whole segmentation.
#' @export
display_segmentation <- function(output, image_name, stage = 1, colourblind_palette = FALSE){
  #Update working directory.
  setwd(output)
  #Select segmentation to analyse.
  if(stage == 1){
    segmentation_name <- "_primary_segmentation.RData"
  } else if(stage == 2){
    segmentation_name <- "_secondary_segmentation.RData"
  } else {
    stop("Please enter a valid stage number (1 for primary, 2 for secondary).")
  }
  #Get file names.
  file_names <- list.files(pattern = image_name)
  #Check if primary segmentation exists.
  if(paste(image_name, segmentation_name, sep = "") %in% file_names){
    segmentation <- readRDS(paste(output, image_name, segmentation_name, sep = ""))
    f <- segmentation[[length(segmentation)]]
    m <- nrow(f)
    n <- ncol(f)
    images <- segmentation[-length(segmentation)]

    #Define pre-set RGB colours.
    if(colourblind_palette){
      colour_list <- rbind(
        c(51, 34, 136),
        c(17, 19, 51),
        c(68, 170, 153),
        c(136, 204, 238),
        c(221, 204, 119),
        c(204, 102, 119),
        c(170, 68, 153),
        c(136, 34, 85)
      ) / 255
    } else {
      colour_list <- rbind(
        c(255, 153, 153),
        c(255, 204, 153),
        c(255, 255, 153),
        c(153, 255, 153),
        c(153, 255, 255),
        c(153, 153, 255),
        c(204, 153, 255),
        c(255, 153, 204)
      ) / 255
    }
    adjacent_colour <- c(2,3,4,5,6,7,8,1)
    current_colour <- 1

    #Set up colour image.
    whole_colour <- array(f, dim = c(m, n, 3))

    #Add objects to colour image.
    for(image in images){
      #Update colours.
      cols <- colour_list[current_colour,]
      current_colour <- adjacent_colour[current_colour]
      #Fill images.
      for(i in 1:m){
        for(j in 1:n){
          if(image[i,j] == 1){
            whole_colour[i,j,] <- cols
          }
        }
      }
    }
    #Display whole colour image.
    grid::grid.newpage()
    grid::grid.raster(whole_colour)
  } else {
    #Output error message.
    stop("Segmentation not found.")
  }
}

#' Display all boundaries found on a given image.
#'
#' This function plots the boundaries of all objects found from TOBLERONE.
#'
#' @param image_name A string representing the name of the image.
#' @param output A string denoting the output file location. As defined in TOBLERONE, this is the location of the segmentation outputs.
#' @param colour A vector of length 3 representing the RGB colour code that defines the colour of the displayed loops. Defaults to red.
#' @return Displays an image of the boundaries found from TOBLERONE.
#' @export
display_boundaries <- function(output, image_name, colour = c(255, 0, 0)){
  #Update working directory.
  setwd(output)
  #Get file names.
  file_names <- list.files(pattern = paste(image_name, "_loop_image_", sep = ""))
  #Check if at least one loop exists.
  if(paste(image_name, "_loop_image_1.tif", sep = "") %in% file_names){
    #Set up colour image.
    whole_colour <- tiff::readTIFF(paste(image_name, "_loop_image_1.tif", sep = ""))

    #Determine number of loops.
    number_of_loops <- length(file_names)

    #Add objects to colour image.
    if(max(colour) > 255 || min(colour) < 0){
      stop("Colour must be RGB with all values between 0 and 255.")
    } else {
      #Update colour code.
      colour <- colour / 255

      #Add each loop to image if it exists.
      for(i in 1:number_of_loops){
        #Try to find loop.
        loop_name <- list.files(pattern = paste(image_name, "_loop_", i, sep = ""))
        #If loop is found, read it and add to image.
        if(length(loop_name) > 0){
          loop <- read.csv(paste(image_name, "_loop_", i, ".csv", sep = ""))
          for(j in 1:nrow(loop)){
            whole_colour[loop[j,1], loop[j,2], ] <- colour
          }
        }
      }

      #Display whole colour image.
      grid::grid.newpage()
      grid::grid.raster(whole_colour)
    }
  } else {
    #Output error message.
    stop("No boundaries found.")
  }
}

#' Calculates the Generalised Polarisation (GP) image of two given images.
#'
#' @param image_1 A matrix representing the first image (preferably ordered channel).
#' @param image_2 A matrix representing the second image (preferably disordered channel).
#' @param normalise Boolean value dictating whether GP values should be normalised between 0 and 1. Defaults to FALSE.
#' @param display Boolean value dictating whether GP image will be displayed. Defaults to FALSE.
#' @return The GP image in matrix form.
#' @export
create_gp_image <- function(image_1, image_2, normalise = FALSE, display = FALSE){
  #Get rows and columns.
  m <- nrow(image_1)
  n <- ncol(image_1)
  if(m == nrow(image_2) && n == nrow(image_2)){
    gp_image <- matrix(0L, nrow = m, ncol = n)
    for(i in 1:m){
      for(j in 1:n){
        sum <- image_1[i,j] + image_2[i,j]
        if(sum == 0){
          gp_image[i,j] <- 0
        } else {
          gp_image[i,j] <- (image_1[i,j] - image_2[i,j]) / sum
        }
      }
    }
    if(normalise){
      if(min(gp_image) == max(gp_image)){
        gp_image <- matrix(0L, nrow = m, ncol = n)
      } else {
        gp_image <- (gp_image - min(gp_image)) / (max(gp_image) - min(gp_image))
      }
    }
    if(display){
      if(min(gp_image) == max(gp_image)){
        normalised_gp_image <- matrix(0L, nrow = m, ncol = n)
      } else {
        normalised_gp_image <- (gp_image - min(gp_image)) / (max(gp_image) - min(gp_image))
      }
      grid::grid.newpage()
      grid::grid.raster(normalised_gp_image)
    }
    return(gp_image)
  } else {
    stop("Image sizes do not match.")
  }
}

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

#' Calculates the Generalised Polarisation (GP) image of two given images and takes the mean GP value of each loop.
#'
#' @param image_1 A matrix representing the first image (preferably ordered channel).
#' @param image_2 A matrix representing the second image (preferably disordered channel).
#' @param image_name A string representing the name of the image originally used in TOBLERONE.
#' @param output A string denoting the output file location as used in TOBLERONE originally.
#' @return A vector of the average GP values from each loop.
#' @export
get_mean_gp <- function(image_1, image_2, image_name, output){
  #Calculate GP image.
  gp_image <- create_gp_image(image_1, image_2)

  #Find all loops.
  setwd(output)
  #Get file names.
  file_names <- list.files(pattern = paste(image_name, "_loop_", sep = ""))
  image_names <- list.files(pattern = paste(image_name, "_loop_image_", sep = ""))
  file_names <- file_names[!file_names %in% image_names]

  #Iterate over all loop labels.
  mean_gp_values <- c()
  for(loop_label in file_names){
    #Upload the loop.
    loop <- read.csv(loop_label)

    #Create list of GP values.
    gp_values <- c()
    for(pixel in 1:nrow(loop)){
      gp_values <- append(gp_values, gp_image[loop[pixel,1], loop[pixel,2]])
    }

    #Add mean GP value to list.
    mean_gp_values <- append(mean_gp_values, mean(gp_values))
  }
  #Return all mean GP values.
  return(mean_gp_values)
}

#' Perform 3D TOBLERONE on a given stack.
#'
#' This function performs TOBLERONE on a given stack of grayscale images. It is recommended
#' that the image is represented in matrix form with all values between 0 and 1.
#'
#' @param stack An array representing the stack of images to be analysed.
#' @param stack_name A string representing the name of the stack, of which all files will be saved under.
#' @param output A string denoting the output file location. All results and calculations will be output here.
#' @return A set of images corresponding to the pixel coordinates of the segmentation found by TOBLERONE.
#' @export
toblerone3D <- function(stack, stack_name, output, use_image_decomposition = TRUE, brightness_scale = 1, colourblind_palette = FALSE){
  #Define empty inputs list.
  inputs <- data.frame()
  
  #Define pre-set RGB colours.
  if(colourblind_palette){
    colour_list <- (rbind(
      c(51, 34, 136),
      c(17, 19, 51),
      c(68, 170, 153),
      c(136, 204, 238),
      c(221, 204, 119),
      c(204, 102, 119),
      c(170, 68, 153),
      c(136, 34, 85)
    ) / 255)
  } else {
    colour_list <- rbind(
      c(255, 153, 153),
      c(255, 204, 153),
      c(255, 255, 153),
      c(153, 255, 153),
      c(153, 255, 255),
      c(153, 153, 255),
      c(204, 153, 255),
      c(255, 153, 204)
    ) / 255
  }
  adjacent_colour <- c(2,3,4,5,6,7,8,1)
  
  #Define offsets used throughout (3D).
  offsets <- cbind(c(-1, -1, 0, 1, 1, 1, 0, -1, 0, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 0),
                   c(0, -1, -1, -1, 0, 1, 1, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0),
                   c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1))
  
  #Generate filtrations matrix and coordinates data frame.
  #Upload stack, which represents the density function, f. Adjust brightness and normalise.
  dimensions <- dim(stack)
  m <- dimensions[1]
  n <- dimensions[2]
  frames <- dimensions[3]
  stack <- stack * brightness_scale
  stack <- stack + (stack > 1) * 1 * (array(1, dim = c(m, n, frames)) - stack)
  f <- stack / max(stack)
  
  #Set working directory.
  setwd(output)
  
  #Determine whether the filtrations matrix and coordinates data frame have already been saved at the given location.
  #If so, use them, otherwise calculate them.
  #Extract all files associated with the image name.
  file_names <- list.files(pattern = stack_name)
  if(paste(stack_name, "_filtrations.RData", sep = "") %in% file_names && paste(stack_name, "_coordinates.csv", sep = "") %in% file_names){
    filtrations <- readRDS(paste(output, stack_name, "_filtrations.RData", sep = ""))
    coordinates <- as.data.frame(read.csv(paste(output, stack_name, "_coordinates.csv", sep = "")))
    old_filtrations <- filtrations
    old_coordinates <- coordinates
  } else {
    #Give output message.
    cat("Calculating filtrations...")
    #Determine count of each intensity value.
    scaled_f <- floor(f * 255)
    count <- do.call(c, lapply(0:255, function(i){
      return(sum(scaled_f == i))
    }))
    #Perform discretised image decomposition.
    if(use_image_decomposition){
      #Calculate break-off intensity.
      mean_count <- mean(count)
      for(i in 1:(length(count) - 1)){
        if(count[i] > mean_count && count[i + 1] < mean_count){
          break_off <- i
          break
        }
      }
    } else {
      #Set break-off to 0.
      break_off <- 0
    }
    
    #Record coordinates.
    active <- which(scaled_f + 1 > break_off)
    vals <- scaled_f[active]
    framenumbers <- (active - 1) %/% (m * n) + 1
    active <- (active - 1) %% (m * n) + 1
    coordinates <- data.frame(cbind(active %% m + (active %% m == 0) * m, (active - 1) %/% m + 1, framenumbers, vals))
    colnames(coordinates) <- c("i", "j", "k", "threshold_number")
    coordinates <- coordinates[order(coordinates$threshold_number, decreasing = TRUE),]
    coordinates <- coordinates[,1:3]
    
    #Create filtrations array.
    filtrations <- array(Inf, dim = c(m, n, frames))
    filtrations[as.matrix(coordinates)] <- 1:nrow(coordinates)
    
    #The entry (i,j, k) of filtrations now corresponds to the filtration value of voxel (i,j,k).
    #Save both the filtration matrix and coordinates data frame for future use.
    saveRDS(filtrations, paste(output, stack_name, "_filtrations.RData", sep = ""))
    write.csv(coordinates, paste(output, stack_name, "_coordinates.csv", sep = ""), row.names = FALSE)
    old_filtrations <- filtrations
    old_coordinates <- coordinates
    
  }
  
  #Primary Segmentation
  #Repeat until suitable threshold is chosen.
  old_filtrations <- filtrations
  old_coordinates <- coordinates
  filtrations <- array(Inf, dim = c(m + 2, n + 2, frames + 2))
  filtrations[2:(m + 1), 2:(n + 1), 2:(frames + 1)] <- old_filtrations
  cat("\nPrimary Segmentation: Enter a number between 0 and 1 as a persistence threshold.")
  first_time <- TRUE
  number_of_voxels <- nrow(coordinates)
  while(TRUE){
    #Get persistence threshold.
    if(first_time){
      input <- readline(prompt = "Enter number: ")
      first_time <- FALSE
    } else {
      input <- readline(prompt = "Enter number, S to show video, or type Y if finished: ")
    }
    inputs <- rbind(inputs, c(1, input))
    if(tolower(input) == "y"){
      break
    } else if(tolower(input) == "s"){
      #Show video.
      current_frame <- c(1, frames + 1, 2 * frames + 1)
      sleeptime <- 5 / frames
      for(i in 1:frames){
        grid::grid.newpage()
        grid::grid.raster(video_colour[,,current_frame])
        current_frame <- current_frame + 1
        Sys.sleep(sleeptime)
      }
    } else {
      #Set persistence threshold.
      persistence_threshold <- as.numeric(input)
      
      #Undertake image segmentation via persistent homology.
      #Create arrays: r contains the roots of each voxel, e contains the set number to which each voxel currently belongs and
      #U contains the sets of voxels.
      r <- rep(0, number_of_voxels)
      e <- rep(0, number_of_voxels)
      U <- list()
      
      #Iterate over all pixels starting from the highest density. Perform persistent homology to generate the partition.
      set_number <- 0
      for(i in 1:number_of_voxels){
        #Extract coordinates of i.
        coordinate_1 <- coordinates[i,1]
        coordinate_2 <- coordinates[i,2]
        coordinate_3 <- coordinates[i,3]
        
        #Let N be the set of neighbours of i with index lower than i.
        neighbour_filtration <- filtrations[cbind(coordinate_1 + offsets[,1] + 1, coordinate_2 + offsets[,2] + 1, coordinate_3 + offsets[,3] + 1)]
        N <- neighbour_filtration[neighbour_filtration < i]
        
        #Check length of N.
        if(length(N) == 0){
          #If N is empty, then i is a peak. Create a new subset in U containing just i. Update r and e.
          r[i] <- i
          set_number <- set_number + 1
          e[i] <- set_number
          U[[set_number]] <- c(i)
        }else{
          #Organise the neighbours in N in descending order of the density of their roots.
          R <- r[N]
          neighbours <- as.data.frame(cbind(R, N))
          neighbours <- neighbours[order(neighbours$R),]
          
          #If N is non-empty, there must be at least one neighbour of higher density. Cycle through all neighbours.
          potentials <- unique(neighbours[f[as.matrix(coordinates[neighbours[,1],1:3])] <= f[coordinate_1,coordinate_2,coordinate_3] + persistence_threshold,1])
          
          if(length(potentials) == 0){
            #Node does not connect to any neighbours. Create new set for node.
            r[i] <- i
            set_number <- set_number + 1
            e[i] <- set_number
            U[[set_number]] <- c(i)
          }else{
            j <- potentials[1]
            U[[e[j]]] <- append(U[[e[j]]], i)
            r[i] <- j
            e[i] <- e[j]
            if(length(potentials) > 1){
              #Merge the sets of the potentials.
              for(k in potentials[-1]){
                e_k <- U[[e[k]]]
                U[[e[j]]] <- append(U[[e[j]]], U[[e[k]]])
                e[e_k] <- e[j]
                r[e_k] <- r[j]
              }
            }
          }
        }
      }
      
      #Output only the sets with root density greater than the persistence threshold.
      #If using whole image, add all objects to the same image.
      videos <- list()
      
      #Set up colour image.
      video_colour <- array(stack, dim = c(m, n, 3 * frames))
      
      #Determine objects.
      found <- 0
      current_colour <- 1
      for(i in unique(r)){
        if(f[coordinates[i,1], coordinates[i,2], coordinates[i,3]] >= persistence_threshold){
          #Initialise new image.
          found <- found + 1
          video <- array(0, dim = c(m, n, frames))
          #Update colours.
          cols <- colour_list[current_colour,]
          current_colour <- adjacent_colour[current_colour]
          #Fill images.
          video[as.matrix(coordinates[U[[e[i]]],])] <- 1
          video_colour[as.matrix(coordinates[U[[e[i]]],])] <- cols[1]
          video_colour[as.matrix(cbind(coordinates[U[[e[i]]],1:2], coordinates[U[[e[i]]],3] + frames))] <- cols[2]
          video_colour[as.matrix(cbind(coordinates[U[[e[i]]],1:2], coordinates[U[[e[i]]],3] + 2 * frames))] <- cols[3]
          #Store image.
          videos[[found]] <- video
        }
      }
      #Output number of objects found and time taken.
      if(found == 1){
        cat("Found 1 object.")
      } else{
        cat("Found ", found, " objects.")
      }
      #Save video.
      saveRDS(videos, paste(output, "Segmentation.RData"))
      current_frame <- c(1, frames + 1, 2 * frames + 1)
      for(i in 1:(frames)){
        writeTIFF(video_colour[,,current_frame], paste(output, stack_name, "_frame_", i, ".tif", sep = ""))
        current_frame <- current_frame + 1
      } 
    }
  }
}