## code to prepare `marsupial_tree` dataset goes here
marsupial_tree <- ape::read.nexus("data-raw/marsupials.trees")
usethis::use_data(marsupial_tree, overwrite = TRUE)
