## code to prepare `marsupials` dataset goes here

marsupials = read.table("data-raw/marsupials.txt", sep = "\t", header = T, stringsAsFactors = F)
usethis::use_data(marsupials, overwrite = TRUE)
