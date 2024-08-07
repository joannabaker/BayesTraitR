Artiodactyl = read.table("data-raw/Artiodactyl.txt", sep = "\t", header = F)
  colnames(Artiodactyl) = c("tip_label", "multistate")

ArtiodactylMLIn = readLines("data-raw/ArtiodactylMLIn.txt")

BirdHetCom = readLines("data-raw/BirdHetCom.txt")

BirdTerritory = read.table("data-raw/BirdTerritory.txt", sep = "\t", header = F)
colnames(BirdTerritory) = c("tip_label", "Territory")

HarmonicMeanLh = readLines("data-raw/HarmonicMeanLh.txt")

MammalBody = read.table("data-raw/MammalBody.txt", sep = "\t", header = F)
colnames(MammalBody) = c("tip_label", "Body")

MammalBrainBody = read.table("data-raw/MammalBrainBody.txt", sep = "\t", header = F)
colnames(MammalBrainBody) = c("tip_label", "Brain", "Body")

MammalBrainBodyGt = read.table("data-raw/MammalBrainBodyGt.txt", sep = "\t", header = F)
colnames(MammalBrainBodyGt) = c("tip_label", "Brain", "Body", "Gt")

MammalBrainBodySampleData = read.table("data-raw/MammalBrainBodySampleData.txt", sep = "\t", header = F)

MammalModelB = read.table("data-raw/MammalModelB.txt", sep = "\t", header = F)
colnames(MammalModelB) = c("tip_label", "trend")

Marsupials = read.table("data-raw/Marsupials.txt", sep = "\t", header = F)
colnames(MammalModelB) = c("tip_label", "Body")


NortheastBantu = read.table("data-raw/NortheastBantu.txt", sep = "\t", header = F)
colnames(NortheastBantu) = c("tip_label", "long", "lat")

Primates = read.table("data-raw/Primates.txt", sep = "\t", header = F)
colnames(Primates) = c("tip_label", "trait1", "trait2")




usethis::use_data(Artiodactyl, ArtiodactylMLIn, BirdHetCom, BirdTerritory, HarmonicMeanLh, MammalBody, MammalBrainBody, MammalBrainBodyGt, MammalBrainBodySampleData, MammalModelB, Marsupials, NortheastBantu, Primates, overwrite = TRUE)
