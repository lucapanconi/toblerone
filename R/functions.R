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
toblerone <- function(image, image_name, output, failsafe = TRUE, use_image_decomposition = TRUE, check_objects = TRUE,
                      base_window_size = 1, identify_boundary = TRUE, smooth_image = TRUE, colourblind_palette = FALSE){
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
  #Upload image, which represents the density function, f.
  f <- image/max(image)
  m <- nrow(f)
  n <- ncol(f)
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
    #Perform discretised image decomposition.
    if(use_image_decomposition){
      #Determine count of each intensity value.
      scaled_f <- floor(f * 255)
      values <- 0:255
      count <- rep(0, 256)
      for(i in 1:m){
        for(j in 1:n){
          count[scaled_f[i,j] + 1] <- count[scaled_f[i,j] + 1] + 1
        }
      }
      #Calculate break-off intensity.
      mean_count <- mean(count)
      for(i in 1:(length(count) - 1)){
        if(count[i] > mean_count && count[i + 1] < mean_count){
          break_off <- i
          break
        }
      }
      #Record coordinates.
      coordinates <- data.frame()
      for(i in 1:m){
        for(j in 1:n){
          if(scaled_f[i,j] + 1 > break_off){
            coordinates <- rbind(coordinates, c(i,j,scaled_f[i,j]))
          }
        }
      }
      colnames(coordinates) <- c("i", "j", "threshold_number")
      coordinates <- coordinates[order(coordinates$threshold_number, decreasing = TRUE),]
      coordinates <- coordinates[,1:2]

      #Create filtrations matrix.
      filtrations <- matrix(Inf, nrow = m, ncol = n)
      number_of_pixels <- nrow(coordinates)
      for(i in 1:number_of_pixels){
        filtrations[coordinates[i,1], coordinates[i,2]] <- i
      }
    } else {
      #Determine each unique threshold to test and organise pixels into vectors corresponding to their threshold values. Use this to
      #determine a filtration value for each pixel.
      number_of_pixels <- m * n
      entries <- c()
      for(i in 1:m){
        entries <- append(entries, f[i,])
      }
      thresholds <- sort(unique(entries), decreasing = TRUE)
      number_of_thresholds <- length(thresholds)
      threshold_numbers <- 1:number_of_thresholds
      coordinates <- data.frame()
      for(i in 1:m){
        for(j in 1:n){
          for(k in threshold_numbers){
            if(f[i,j] == thresholds[k]){
              coordinates <- rbind(coordinates, c(i, j, k))
              break
            }
          }
        }
      }
      colnames(coordinates) <- c("i", "j", "threshold_number")
      coordinates <- coordinates[order(coordinates$threshold_number),]
      #The row of coordinates now corresponds to the filtration value of the pixel given by coordinates (i,j). We no longer need the
      #third column, so we remove it.
      coordinates <- coordinates[,1:2]
      number_of_pixels <- m * n
      filtrations <- matrix(0L, nrow = m, ncol = n)
      for(i in 1:number_of_pixels){
        filtrations[coordinates[i,1], coordinates[i,2]] <- i
      }
    }
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
          N <- c()
          for(k in 1:8){
            neighbour_filtration <- filtrations[coordinate_1 + offsets[k,1] + 1, coordinate_2 + offsets[k,2] + 1]
            if(neighbour_filtration < i){
              N <- append(N, neighbour_filtration)
            }
          }

          #Check length of N.
          if(length(N) == 0){
            #If N is empty, then i is a peak. Create a new subset in U containing just i. Update r and e.
            r[i] <- i
            set_number <- set_number + 1
            e[i] <- set_number
            U[[set_number]] <- c(i)
          }else{
            #Organise the neighbours in N in descending order of the density of their roots.
            R <- c()
            for(k in N){
              R <- append(R, r[k])
            }
            neighbours <- as.data.frame(cbind(R, N))
            neighbours <- neighbours[order(neighbours$R),]

            #If N is non-empty, there must be at least one neighbour of higher density. Cycle through all neighbours.
            check <- TRUE
            potentials <- c()
            for(j in 1:nrow(neighbours)){
              if(f[coordinate_1,coordinate_2] + persistence_threshold >= f[coordinates[neighbours[j,1],1],coordinates[neighbours[j,1],2]]){
                potentials <- append(potentials, neighbours[j,1])
                check <- FALSE
              }
            }
            potentials <- unique(potentials)
            if(check){
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
                  for(l in e_k){
                    e[l] <- e[j]
                    r[l] <- r[j]
                  }
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
            for(k in U[[e[i]]]){
              image[coordinates[k,1], coordinates[k,2]] <- 1
              whole_colour[coordinates[k,1], coordinates[k,2],] <- cols
            }
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
              for(k in 1:8){
                display_image_outline[i + offsets[k,1],j + offsets[k,2],] <- c(1,0,0)
              }
            }
          }
        }
        #Fill in colour.
        for(i in 1:m){
          for(j in 1:n){
            if(image[i,j] == 1){
              display_image[i,j,] <- cols
              display_image_outline[i,j,] <- cols
            }
          }
        }

        #Display the image.
        grid::grid.newpage()
        grid::grid.raster(display_image_outline)

        #Reset window size and create initial window.
        window_size_i <- base_window_size
        window_size_j <- base_window_size

        #Determine i and j boundaries.
        i_values <- c()
        j_values <- c()
        for(i in 1:m){
          for(j in 1:n){
            if(image[i,j] == 1){
              i_values <- append(i_values, i)
              j_values <- append(j_values, j)
            }
          }
        }
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
          new_coordinates <- data.frame()
          for(i in 1:(i_max - i_min + 1)){
            for(j in 1:(j_max - j_min + 1)){
              if(new_filtrations[i,j] != Inf){
                new_coordinates <- rbind(new_coordinates, c(old_coordinates[new_filtrations[i,j],1] - i_min + 1,
                                                            old_coordinates[new_filtrations[i,j],2] - j_min + 1,
                                                            new_filtrations[i,j]))
              }
            }
          }
          #Re-order coordinates in ascending filtration value.
          colnames(new_coordinates) <- c("i", "j", "threshold_number")
          new_coordinates <- new_coordinates[order(new_coordinates$threshold_number),]
          #Display current image.
          new_f <- f[i_min:i_max, j_min:j_max]
          new_image <- image[i_min:i_max, j_min:j_max]
          new_display_image <- display_image[i_min:i_max, j_min:j_max,]
          grid::grid.newpage()
          grid::grid.raster(new_display_image)

          #Calculate new number of pixels.
          number_of_pixels <- nrow(new_coordinates)

          #Update filtrations matrix.
          for(pixel in 1:number_of_pixels){
            new_filtrations[new_coordinates[pixel,1], new_coordinates[pixel,2]] <- pixel
          }

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
            N <- c()
            for(k in 1:8){
              neighbour_filtration <- filtrations[coordinate_1 + offsets[k,1] + 1, coordinate_2 + offsets[k,2] + 1]
              if(neighbour_filtration < i){
                N <- append(N, neighbour_filtration)
              }
            }

            #Check length of N.
            if(length(N) == 0){
              #If N is empty, then i is a peak. Create a new subset in U containing just i. Update r and e.
              r[i] <- i
              set_number <- set_number + 1
              e[i] <- set_number
              U[[set_number]] <- c(i)
            }else{
              #Organise the neighbours in N in descending order of the density of their roots.
              R <- c()
              for(k in N){
                R <- append(R, r[k])
              }
              neighbours <- as.data.frame(cbind(R, N))
              neighbours <- neighbours[order(neighbours$R),]

              #If N is non-empty, there must be at least one neighbour of higher density. Cycle through all neighbours.
              check <- TRUE
              potentials <- c()
              for(j in 1:nrow(neighbours)){
                if(new_f[coordinate_1, coordinate_2] + persistence_threshold >= new_f[coordinates[neighbours[j,1],1],coordinates[neighbours[j,1],2]]){
                  potentials <- append(potentials, neighbours[j,1])
                  check <- FALSE
                }
              }
              potentials <- unique(potentials)
              if(check){
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
                    for(l in e_k){
                      e[l] <- e[j]
                      r[l] <- r[j]
                    }
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
              for(k in U[[e[i]]]){
                new_image[coordinates[k,1], coordinates[k,2]] <- 1
                whole_new_image[coordinates[k,1], coordinates[k,2]] <- 1
                whole_colour[coordinates[k,1], coordinates[k,2],] <- cols
                display_image[coordinates[k,1] + i_min - 1, coordinates[k,2] + j_min - 1,] <- cols
              }
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
            input <- readline(prompt = "Enter number, type Y to take whole image, type S to take individual segments, or type N to discard: ")
            inputs <- rbind(inputs, c(2, input))
            if(tolower(input) == "y"){
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
            } else if(tolower(input) == "s"){
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
              #persistenc threshold and perform TOBLERONE again.
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
                new_coordinates <- data.frame()
                for(i in 1:(i_max - i_min + 1)){
                  for(j in 1:(j_max - j_min + 1)){
                    if(new_filtrations[i,j] != Inf){
                      new_coordinates <- rbind(new_coordinates, c(old_coordinates[new_filtrations[i,j],1] - i_min + 1,
                                                                  old_coordinates[new_filtrations[i,j],2] - j_min + 1,
                                                                  new_filtrations[i,j]))
                    }
                  }
                }
                #Re-order coordinates in ascending filtration value.
                colnames(new_coordinates) <- c("i", "j", "threshold_number")
                new_coordinates <- new_coordinates[order(new_coordinates$threshold_number),]
                #Display current image.
                new_f <- f[i_min:i_max, j_min:j_max]
                new_image <- image[i_min:i_max, j_min:j_max]
                new_display_image <- display_image[i_min:i_max, j_min:j_max,]
                grid::grid.newpage()
                grid::grid.raster(new_display_image)

                #Calculate new number of pixels.
                number_of_pixels <- nrow(new_coordinates)

                #Update filtrations matrix.
                for(pixel in 1:number_of_pixels){
                  new_filtrations[new_coordinates[pixel,1], new_coordinates[pixel,2]] <- pixel
                }

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

                #Extract the part of the image that corresponds to that object plus a given window.
                new_filtrations <- old_filtrations[i_min:i_max, j_min:j_max]
                new_coordinates <- data.frame()
                for(i in 1:(i_max - i_min + 1)){
                  for(j in 1:(j_max - j_min + 1)){
                    if(new_filtrations[i,j] != Inf){
                      new_coordinates <- rbind(new_coordinates, c(old_coordinates[new_filtrations[i,j],1] - i_min + 1,
                                                                  old_coordinates[new_filtrations[i,j],2] - j_min + 1,
                                                                  new_filtrations[i,j]))
                    }
                  }
                }
                #Re-order coordinates in ascending filtration value.
                colnames(new_coordinates) <- c("i", "j", "threshold_number")
                new_coordinates <- new_coordinates[order(new_coordinates$threshold_number),]
                #Display current image.
                new_f <- f[i_min:i_max, j_min:j_max]
                new_image <- image[i_min:i_max, j_min:j_max]
                new_display_image <- display_image[i_min:i_max, j_min:j_max,]
                grid::grid.newpage()
                grid::grid.raster(new_display_image)

                #Calculate new number of pixels.
                number_of_pixels <- nrow(new_coordinates)

                #Update filtrations matrix.
                for(pixel in 1:number_of_pixels){
                  new_filtrations[new_coordinates[pixel,1], new_coordinates[pixel,2]] <- pixel
                }

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
                  N <- c()
                  for(k in 1:8){
                    neighbour_filtration <- filtrations[coordinate_1 + offsets[k,1] + 1, coordinate_2 + offsets[k,2] + 1]
                    if(neighbour_filtration < i){
                      N <- append(N, neighbour_filtration)
                    }
                  }

                  #Check length of N.
                  if(length(N) == 0){
                    #If N is empty, then i is a peak. Create a new subset in U containing just i. Update r and e.
                    r[i] <- i
                    set_number <- set_number + 1
                    e[i] <- set_number
                    U[[set_number]] <- c(i)
                  }else{
                    #Organise the neighbours in N in descending order of the density of their roots.
                    R <- c()
                    for(k in N){
                      R <- append(R, r[k])
                    }
                    neighbours <- as.data.frame(cbind(R, N))
                    neighbours <- neighbours[order(neighbours$R),]

                    #If N is non-empty, there must be at least one neighbour of higher density. Cycle through all neighbours.
                    check <- TRUE
                    potentials <- c()
                    for(j in 1:nrow(neighbours)){
                      if(new_f[coordinate_1, coordinate_2] + persistence_threshold >= new_f[coordinates[neighbours[j,1],1],coordinates[neighbours[j,1],2]]){
                        potentials <- append(potentials, neighbours[j,1])
                        check <- FALSE
                      }
                    }
                    potentials <- unique(potentials)
                    if(check){
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
                          for(l in e_k){
                            e[l] <- e[j]
                            r[l] <- r[j]
                          }
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
                    for(k in U[[e[i]]]){
                      new_image[coordinates[k,1], coordinates[k,2]] <- 1
                      whole_new_image[coordinates[k,1], coordinates[k,2]] <- 1
                      whole_colour[coordinates[k,1], coordinates[k,2],] <- cols
                      display_image[coordinates[k,1] + i_min - 1, coordinates[k,2] + j_min - 1,] <- cols
                    }
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

    secondary_segmentation <- new_images
    secondary_segmentation[[length(secondary_segmentation) + 1]] <- f
    saveRDS(secondary_segmentation, paste(output, image_name, "_secondary_segmentation.RData", sep = ""))
    saveRDS(windows, paste(output, image_name, "_analysis_windows.RData", sep = ""))
    stage <- stage + 1
  } else if(checks[2] == 1){
    new_images <- secondary_segmentation[-length(secondary_segmentation)]
  }
  #Display the image.
  grid::grid.newpage()
  grid::grid.raster(display_image)

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
      adjacent_pixels <- cbind(c(0,-1,0,1), c(1,0,-1,0))
      #If there are no images left (all have been deleted), stop.
      images_available <- TRUE
      number_of_images <- length(images)
      if(number_of_images == 0){
        images_available <- FALSE
      }else{
        #Iterate over each image. Scan each pixel left to right, then top to bottom, upon hitting a bright pixel, perform grid
        #sweep and add to loops list.
        loops <- list()
        found <- 0
        for(index in 1:number_of_images){
          empty <- c()
          image <- images[[index]]
          window <- windows[[index]]
          #Construct window for analysis.
          i_min <- window[1,1]
          i_max <- window[1,2]
          j_min <- window[2,1]
          j_max <- window[2,2]
          image <- rbind(0, cbind(0, image[i_min:i_max, j_min:j_max], 0), 0)
          new_m <- i_max - i_min + 2
          new_n <- j_max - j_min + 2

          #Smooth image.
          if(smooth_image){
            changed <- TRUE
            while(changed){
              changed <- FALSE
              for(k in 2:new_m){
                for(l in 2:new_n){
                  if(image[k,l] == 1){
                    total <- 0
                    for(p in 1:4){
                      if(image[k + adjacent_pixels[p,1],l + adjacent_pixels[p,2]] == 0){
                        total <- total + 1
                      }
                    }
                    if(total >= 3){
                      image[k,l] <- 0
                      changed <- TRUE
                    }
                  }
                }
              }
            }
          }

          #Find loops.
          for(k in 2:new_m){
            for(l in 2:new_n){
              if(image[k,l] == 1){
                #New loop has been found.
                found <- found + 1

                #Initialise loop node list.
                i <- k
                j <- l
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
                  if(image[square_i, square_j] == 1){
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
                #Upon finishing, we remove everything inside the loop from the image.
                #Starting from the first pixel in the loop, check all adjacent pixels - if they are bright, record them, delete
                #them and then check their adjacent pixels. Repeat until no more unchecked pixels remain.
                unchecked <- data.frame(loop_nodes_x, loop_nodes_y)
                record <- rbind(FALSE, cbind(FALSE, matrix(TRUE, nrow = new_m, ncol = new_n), FALSE), FALSE)
                while(nrow(unchecked) > 0){
                  coordinate <- unchecked[1,]
                  image[coordinate[1,1], coordinate[1,2]] <- 0
                  coordinate <- coordinate + 1
                  for(s in 1:8){
                    new_coordinate <- coordinate + offsets[s,]
                    if(record[new_coordinate[1,1], new_coordinate[1,2]]){
                      record[new_coordinate[1,1], new_coordinate[1,2]] <- FALSE
                      new_coordinate <- new_coordinate - 1
                      if(image[new_coordinate[1,1], new_coordinate[1,2]] == 1){
                        unchecked <- rbind(unchecked, new_coordinate)
                      }
                    }
                  }
                  record[coordinate[1,1], coordinate[1,2]] <- FALSE
                  unchecked <- unchecked[-1,]
                }

                #Save the corresponding loop and remove from layered image.
                loops[[found]] <- cbind(loop_nodes_x + i_min - 2, loop_nodes_y + j_min - 2)
              }
            }
          }
        }
      }
      #Remove empty images from images list.
      if(length(empty) > 0){
        images[-empty]
        windows[-empty]
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
        layer_loops <- list()
        #For each good loop, create a list in layer loops to store that loop and all loops inside it. For each loop,
        #keep adding all adjacent pixels as there own loop until no more are found.
        for(i in 1:found){
          #Create list containing original loop.
          layer_loops[[i]] <- list(good_loops[[i]])
          counter <- 1
          check_image <- rbind(FALSE, cbind(FALSE, matrix(TRUE, nrow = m, ncol = n), FALSE), FALSE)
          while(counter == length(layer_loops[[i]])){
            #Extract latest loop.
            loop <- layer_loops[[i]][[counter]]
            #Create lists of length equal to number of rows and columns. Each entry in the row list contains the columns
            #in which that row appears and each entry in the columns list contains the rows in which that column appears.
            #Create image to store where loop actually is.
            appears_in_row <- vector(mode = "list", length = m)
            appears_in_col <- vector(mode = "list", length = n)
            for(j in 1:nrow(loop)){
              pixel_row <- as.numeric(loop[j,1])
              pixel_col <- as.numeric(loop[j,2])
              appears_in_row[[pixel_row]] <- append(appears_in_row[[pixel_row]], pixel_col)
              appears_in_col[[pixel_col]] <- append(appears_in_col[[pixel_col]], pixel_row)
              check_image[pixel_row + 1, pixel_col + 1] <- FALSE
            }
            #Clean up lists so intersections are not counted twice for consecutive pixels.
            for(j in 1:m){
              pixels_recorded <- sort(appears_in_row[[j]])
              number_of_pixels_in_row <- length(pixels_recorded)
              if(number_of_pixels_in_row > 1){
                appears_in_row[[j]] <- sort(appears_in_row[[j]])
                for(k in number_of_pixels_in_row:2){
                  if(pixels_recorded[k] == pixels_recorded[k-1] + 1){
                    appears_in_row[[j]] <- appears_in_row[[j]][-k]
                  }
                }
              }
            }
            for(j in 1:n){
              pixels_recorded <- sort(appears_in_col[[j]])
              number_of_pixels_in_col <- length(pixels_recorded)
              if(number_of_pixels_in_col > 1){
                appears_in_col[[j]] <- sort(appears_in_col[[j]])
                for(k in number_of_pixels_in_col:2){
                  if(pixels_recorded[k] == pixels_recorded[k-1] + 1){
                    appears_in_col[[j]] <- appears_in_col[[j]][-k]
                  }
                }
              }
            }
            #Initialise new loop.
            new_loop <- data.frame()
            #For each pixel, check immediately adjacent neighbours.
            previous_pixel_row <- as.numeric(loop[nrow(loop),1])
            previous_pixel_col <- as.numeric(loop[nrow(loop),2])
            for(j in 1:nrow(loop)){
              current_pixel_row <- as.numeric(loop[j,1])
              current_pixel_col <- as.numeric(loop[j,2])
              difference_row <- current_pixel_row - previous_pixel_row
              difference_col <- current_pixel_col - previous_pixel_col
              if(difference_row == 1 && difference_col >= 0){
                adjacent_pixels <- cbind(c(0,-1,0,1), c(1,0,-1,0))
              } else if(difference_row == -1 && difference_col <= 0){
                adjacent_pixels <- cbind(c(0,1,0,-1), c(-1,0,1,0))
              } else if(difference_col == 1 && difference_row <= 0){
                adjacent_pixels <- cbind(c(-1,0,1,0), c(0,-1,0,1))
              } else {
                adjacent_pixels <- cbind(c(1,0,-1,0), c(0,1,0,-1))
              }

              #This only works for cells with geometries which are not too complex - any which happen to have random infoldings
              #will not always work.
              for(k in 1:4){
                pixel_row <- as.numeric(current_pixel_row + adjacent_pixels[k,1])
                pixel_col <- as.numeric(current_pixel_col + adjacent_pixels[k,2])
                #Check if pixel is not already in loop.
                if(check_image[pixel_row + 1, pixel_col + 1]){
                  #Check if pixel has odd number of intersections above, below, left and right of it.
                  pix_rec_row <- appears_in_row[[pixel_row]]
                  pix_rec_col <- appears_in_col[[pixel_col]]
                  if(length(pix_rec_row[pix_rec_row > pixel_col]) > 0 && length(pix_rec_row[pix_rec_row < pixel_col]) > 0 &&
                     length(pix_rec_col[pix_rec_col > pixel_row]) > 0 && length(pix_rec_col[pix_rec_col < pixel_row]) > 0){
                    #Update image to include newly found pixels (so they are not counted twice)
                    check_image[pixel_row + 1, pixel_col + 1] <- FALSE
                    new_loop <- rbind(new_loop, c(pixel_row, pixel_col))
                  }
                }
              }
              previous_pixel_row <- current_pixel_row
              previous_pixel_col <- current_pixel_col
            }

            #Add 1 to counter. If the number of loops does not equal counter, then no new loop was found, so we stop.
            counter <- counter + 1

            #If at least one pixel was found, add loop to list of layer loops.
            new_loop_length <- nrow(new_loop)
            if(new_loop_length > 1){
              #Check the loop for any breaks.
              for(point in 1:(new_loop_length - 1)){
                if((new_loop[point,1] - new_loop[point + 1,1]) ** 2 + (new_loop[point,2] - new_loop[point + 1,2]) ** 2 > 2){
                  break
                }
              }
              #If a break is found, repeat the legacy code to find the boundaries.
              if(point != new_loop_length - 1){
                #Recreate image but remove all previously found loops.
                image <- images[[i]]
                window <- windows[[i]]
                for(layer_to_delete in layer_loops[[i]]){
                  for(p in 1:nrow(layer_to_delete)){
                    image[layer_to_delete[p,1], layer_to_delete[p,2]] <- 0
                  }
                }
                original_image <- image
                #Repeat legacy code, adding all newly found loops, until no new loops are found. If multiple loops are
                #found, store all of them in a new list in layer loops. Display all but save individually.
                while(TRUE){
                  #Construct window for analysis.
                  i_min <- window[1,1]
                  i_max <- window[1,2]
                  j_min <- window[2,1]
                  j_max <- window[2,2]
                  image <- rbind(0, cbind(0, image[i_min:i_max, j_min:j_max], 0), 0)
                  new_m <- i_max - i_min + 2
                  new_n <- j_max - j_min + 2

                  #Smooth image.
                  if(smooth_image){
                    changed <- TRUE
                    while(changed){
                      changed <- FALSE
                      for(k in 2:new_m){
                        for(l in 2:new_n){
                          if(image[k,l] == 1){
                            total <- 0
                            for(p in 1:4){
                              if(image[k + adjacent_pixels[p,1],l + adjacent_pixels[p,2]] == 0){
                                total <- total + 1
                              }
                            }
                            if(total >= 3){
                              image[k,l] <- 0
                              changed <- TRUE
                            }
                          }
                        }
                      }
                    }
                  }

                  #Find loops.
                  newly_found_loops <- list()
                  newly_found <- 0
                  for(k in 2:new_m){
                    for(l in 2:new_n){
                      if(image[k,l] == 1){
                        #New loop has been found.
                        newly_found <- newly_found + 1
                        #Initialise loop node list.
                        q <- k
                        j <- l
                        initial_i <- q
                        initial_j <- j
                        loop_nodes_x <- c(q)
                        loop_nodes_y <- c(j)

                        #Initialise first start square and current square.
                        start_square <- 3
                        current_square <- 3

                        #Perform grid sweep.
                        finished <- FALSE
                        while(!finished){
                          #Get coordinates of current square.
                          square_i <- q + offsets[current_square, 1]
                          square_j <- j + offsets[current_square, 2]
                          #Check if current square is bright.
                          if(image[square_i, square_j] == 1){
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
                              q <- square_i
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
                        #Upon finishing, we remove everything inside the loop from the image.
                        #Starting from the first pixel in the loop, check all adjacent pixels - if they are bright, record them, delete
                        #them and then check their adjacent pixels. Repeat until no more unchecked pixels remain.
                        unchecked <- data.frame(loop_nodes_x, loop_nodes_y)
                        record <- rbind(FALSE, cbind(FALSE, matrix(TRUE, nrow = new_m, ncol = new_n), FALSE), FALSE)
                        while(nrow(unchecked) > 0){
                          coordinate <- unchecked[1,]
                          image[coordinate[1,1], coordinate[1,2]] <- 0
                          coordinate <- coordinate + 1
                          for(s in 1:8){
                            new_coordinate <- coordinate + offsets[s,]
                            if(record[new_coordinate[1,1], new_coordinate[1,2]]){
                              record[new_coordinate[1,1], new_coordinate[1,2]] <- FALSE
                              new_coordinate <- new_coordinate - 1
                              if(image[new_coordinate[1,1], new_coordinate[1,2]] == 1){
                                unchecked <- rbind(unchecked, new_coordinate)
                              }
                            }
                          }
                          record[coordinate[1,1], coordinate[1,2]] <- FALSE
                          unchecked <- unchecked[-1,]
                        }

                        #Save the corresponding loop and remove from layered image.
                        newly_found_loops[[newly_found]] <- cbind(loop_nodes_x + i_min - 2, loop_nodes_y + j_min - 2)
                      }
                    }
                  }
                  #If no loops were found, stop, otherwise add all newly found loops to layer loops.
                  if(newly_found == 0){
                    counter <- 0
                    break
                  } else {
                    #Add to layer loops.
                    layer_loops[[i]][[counter]] <- newly_found_loops
                    counter <- counter + 1
                    #Remove all found loops from original image and repeat.
                    image <- original_image
                    for(newly_found_loop in newly_found_loops){
                      for(p in 1:nrow(newly_found_loop)){
                        image[newly_found_loop[p,1], newly_found_loop[p,2]] <- 0
                      }
                    }
                    original_image <- image
                  }
                }
              } else {
                #If no breaks were found, proceed as before.
                layer_loops[[i]][[counter]] <- as.matrix(new_loop)
              }
            }
          }
        }
        saveRDS(layer_loops, paste(output, image_name, "_layer_identification.RData", sep = ""))
        stage <- stage + 1
      }
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
