## code to prepare `marsupials` dataset goes here
marsupials <- read.csv('data-raw/marsupials.csv')
use_data(marsupials, overwrite = TRUE)


marsupials_tree <- read.nexus('data-raw/marsupials_TTOL.trees')
use_data(marsupials_tree, overwrite = TRUE)
