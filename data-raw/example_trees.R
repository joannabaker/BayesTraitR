Artiodactyl_trees <- ape::read.nexus("data-raw/Artiodactyl.trees")
Bird_tree <- ape::read.nexus("data-raw/Bird.trees")
Mammal_trees <- ape::read.nexus("data-raw/Mammal.trees")
Marsupial_tree <- ape::read.nexus("data-raw/Marsupials.trees")
NortheastBantu_tree <- ape::read.nexus("data-raw/NortheastBantu.trees")
Primates_trees <- ape::read.nexus("data-raw/Primates.trees")

usethis::use_data(Artiodactyl_trees, Bird_tree, Mammal_trees, Marsupial_tree, NortheastBantu_tree, Primates_trees, overwrite = TRUE)
