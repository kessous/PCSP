library(countrycode)
library(sp)
library(rgdal)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(elevatr)
library(readr)
library(magrittr)
#library(speciesgeocodeR)
library(raster)
#library(scrubr)
library(maps)
library(abjutils)
#install.packages("rnaturalearthdata")
library(rnaturalearthdata)
#install.packages("flora")
library(flora)
library(vegan)


##################################### GBIF #######################################

raw_occ_gbif<-read_tsv("C:/Users/igor.kessous/0006303-230810091245214.csv",quote="")
#x<-unique(raw_occ_gbif$taxonRank)
sp_vr_occ_gbif<-raw_occ_gbif[!(raw_occ_gbif$taxonRank %in% c('FAMILY','ORDER', 'CLASS','GENUS')), ]

#write.csv(sp_vr_occ_gbif, file = 'gbif_filtered_infraspecies.csv')


# CLEAN THE DATASET to use coordinate
# mind that data often contain errors, so careful inspection and cleaning are necessary! 
# here we'll first remove records of absence or zero-abundance (if any):
names(sp_vr_occ_gbif)
sort(unique(sp_vr_occ_gbif$individualCount))  # notice if some points correspond to zero abundance
sort(unique(sp_vr_occ_gbif$occurrenceStatus))  # check for different indications of "absent", which could be in different languages! and remember that R is case-sensitive
absence_rows <- which(sp_vr_occ_gbif$individualCount == 0 |sp_vr_occ_gbif$occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"))
length(absence_rows)
if (length(absence_rows) > 0) {
  sp_vr_occ_gbif <- sp_vr_occ_gbif[-absence_rows, ]
}

############### CoordinateCleaner ##################
sp_vr_occ_gbif$countryCode <-  countrycode(sp_vr_occ_gbif$countryCode, origin =  'iso2c', destination = 'iso3c')

sp_vr_occ_gbif2<-cc_val(
  sp_vr_occ_gbif,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "clean",
  verbose = TRUE
)


flags <- clean_coordinates(x = sp_vr_occ_gbif2,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros", "countries")) 

summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")
sp_vr_occ_gbif_cl <- sp_vr_occ_gbif2[flags$.summary,]

#setwd('/Users/xbaroc/Documents/R/campos_de_altitude')
#write.csv(occ_splink2_cl, file= 'data_splink_clean.csv')
  

sp_vr_occ_gbif_cl2<-data.frame(filter(sp_vr_occ_gbif_cl, sp_vr_occ_gbif_cl$decimalLatitude <= -21.5))
sp_vr_occ_gbif_cl2<-data.frame(filter(sp_vr_occ_gbif_cl2, sp_vr_occ_gbif_cl2$decimalLatitude >= -22.3))
sp_vr_occ_gbif_cl2<-data.frame(filter(sp_vr_occ_gbif_cl2, sp_vr_occ_gbif_cl2$decimalLongitude <= -46))
sp_vr_occ_gbif_cl2<-data.frame(filter(sp_vr_occ_gbif_cl2, sp_vr_occ_gbif_cl2$decimalLongitude >= -47))

e<-extent(c(-47,-46,-22.31,-21.49))
elev <- get_elev_raster(raster(e, vals = 0), z = 11, prj = "+proj=longlat +datum=WGS84 +no_defs")
extent(elev)
#elev_crop=crop(elev,e)
#plot(elev_crop)
plot_data <- raster::extract(elev_crop, raster::extent(elev_crop), buffer = 10)
plot_df <- data.frame(raster::xyFromCell(elev_crop, 1:ncell(elev_crop)), value = plot_data)
elev2<-fortify(plot_df, maxpixels = 1000000000)

gg <- ggplot() +
  geom_raster(data = elev2, aes(x = x, y = y, fill = value)) + 
  geom_point(data = sp_vr_occ_gbif_cl2, aes(x = decimalLongitude, y = decimalLatitude), colour = "red", size = 0.3)+
  scale_fill_gradient(name = "Elevation (m)\n", low = "white", high = "darkblue", limits=c(500,2000)) +
  coord_equal() +labs(x = "Longitude", y = "Latitude")
gg


#removing occurences with elevation lower than 1000m
sp_vr_occ_gbif_cl2$elevation<-raster::extract(elev_crop,cbind.data.frame(sp_vr_occ_gbif_cl2$decimalLongitude,sp_vr_occ_gbif_cl2$decimalLatitude))
nrow(df)
hist(df$elevation)
df = sp_vr_occ_gbif_cl2 %>% filter(elevation>1000)
df2 = sp_vr_occ_gbif_cl2 %>% filter(elevation>1300)

gg2 <- ggplot() +
  geom_raster(data = elev2, aes(x = x, y = y, fill = value)) + 
  geom_point(data = df, aes(x = decimalLongitude, y = decimalLatitude), colour = "red", size = 0.3)+
  scale_fill_gradient(name = "Elevation (m)\n", low = "white", high = "darkblue", limits=c(500,2000)) +
  coord_equal() +labs(x = "Longitude", y = "Latitude")
gg2


gg3 <- ggplot() +
  geom_raster(data = elev2, aes(x = x, y = y, fill = value)) + 
  geom_point(data = df, aes(x = decimalLongitude, y = decimalLatitude), colour = "black", size = 0.2) + 
  scale_fill_viridis_c(name = "Elevation (m)\n", limits=c(500,2000), direction = -1, option = "plasma") +
  coord_equal() +labs(x = "Longitude", y = "Latitude")
gg3  
  

pdf('map_pcs_plants.pdf',width = 10,height = 10)
print(gg3)
dev.off()

pdf('map_pcs_plants2.pdf',width = 10,height = 10)
print(gg2)
dev.off()

pdf('map_pcs_plants_all.pdf',width = 10,height = 10)
print(gg)
dev.off()

#Including data without coordinate
non_coordinate<-data.frame(subset(sp_vr_occ_gbif, is.na(sp_vr_occ_gbif$decimalLatitude)))
non_coordinate_mg<- subset(non_coordinate, grepl("Minas Gerais|MG", non_coordinate$stateProvince, ignore.case = T))
non_coordinate_sp<- subset(non_coordinate, grepl("PAULO|SP", non_coordinate$stateProvince, ignore.case = T))


non_coordinate_mg2<- subset(non_coordinate_mg, grepl("Caldas|óleo|oleo|andradas|ibitiura de minas|ibitiúra de minas|serrote|pocinhos do rio verde|fazenda ponte alta|campinas|morada dos passaros|morada dos pássaros|vale das antas|bortolan|moinhos", non_coordinate_mg$locality, ignore.case = T))
non_coordinate_sp2<- subset(non_coordinate_sp, grepl("são joão da boa vista|sao joão da boa vista|sao joao da boa vista|são joao da boa vista|fazenda monte alegre|são roque da fartura|sao roque da fartura|ventania|bortolan|fazenda recreio|cachoeira da cascatinha|cachoeira do coqueiro torto|pico do gavião|pico do gaviao|aguas da prata|águas da prata|divinolandia|divinolândia", non_coordinate_sp$locality, ignore.case = T))

non_coordinate_total<-rbind(non_coordinate_mg2, non_coordinate_sp2)


total_vouchers<-rbind(df,non_coordinate_total)

########Including data from JABOT AFR ############
list_jabot_afr<-read.delim(file = "/Users/xbaroc/Desktop/dwca-afr-v1.28/occurrence.txt", header = TRUE, sep = "\t", dec = ".")
list_jabot_afr<-list_jabot_afr[!(list_jabot_afr$taxonRank %in% c('family','genus')), ]
list_jabot_afr$countryCode <- "BR"
list_jabot_afr$countryCode <-  countrycode(list_jabot_afr$countryCode, origin =  'iso2c', destination = 'iso3c')

list_jabot_afr2<-cc_val(
  list_jabot_afr,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  value = "clean",
  verbose = TRUE
)

list_jabot_afr2<-data.frame(filter(list_jabot_afr2, list_jabot_afr2$decimalLatitude <= -21.5))
list_jabot_afr2<-data.frame(filter(list_jabot_afr2, list_jabot_afr2$decimalLatitude >= -22.3))
list_jabot_afr2<-data.frame(filter(list_jabot_afr2, list_jabot_afr2$decimalLongitude <= -46))
list_jabot_afr2<-data.frame(filter(list_jabot_afr2,list_jabot_afr2$decimalLongitude >= -47))

non_coordinate_afr<-data.frame(subset(list_jabot_afr, is.na(list_jabot_afr$decimalLatitude)))
non_coordinate_mg_afr<- subset(non_coordinate_afr, grepl("Minas Gerais|MG", non_coordinate_afr$stateProvince, ignore.case = T))
non_coordinate_sp_afr<- subset(non_coordinate_afr, grepl("PAULO|SP", non_coordinate_afr$stateProvince, ignore.case = T))
non_coordinate_mg2_afr<- subset(non_coordinate_mg_afr, grepl("Caldas|óleo|oleo|andradas|ibitiura de minas|ibitiúra de minas|serrote|pocinhos do rio verde|fazenda ponte alta|campinas|morada dos passaros|morada dos pássaros|vale das antas|bortolan|moinhos", non_coordinate_mg_afr$locality, ignore.case = T))
non_coordinate_mg3_afr<- subset(non_coordinate_mg_afr, grepl("Caldas|óleo|oleo|andradas|ibitiura de minas|ibitiúra de minas|serrote|pocinhos do rio verde|fazenda ponte alta|campinas|morada dos passaros|morada dos pássaros|vale das antas|bortolan|moinhos", non_coordinate_mg_afr$municipality, ignore.case = T))
non_coordinate_sp2_afr<- subset(non_coordinate_sp_afr, grepl("são joão da boa vista|sao joão da boa vista|sao joao da boa vista|são joao da boa vista|fazenda monte alegre|são roque da fartura|sao roque da fartura|ventania|bortolan|fazenda recreio|cachoeira da cascatinha|cachoeira do coqueiro torto|pico do gavião|pico do gaviao|aguas da prata|águas da prata|divinolandia|divinolândia", non_coordinate_sp_afr$locality, ignore.case = T))
non_coordinate_sp3_afr<- subset(non_coordinate_sp_afr, grepl("são joão da boa vista|sao joão da boa vista|sao joao da boa vista|são joao da boa vista|fazenda monte alegre|são roque da fartura|sao roque da fartura|ventania|bortolan|fazenda recreio|cachoeira da cascatinha|cachoeira do coqueiro torto|pico do gavião|pico do gaviao|aguas da prata|águas da prata|divinolandia|divinolândia", non_coordinate_sp_afr$municipality, ignore.case = T))
non_coordinate_total_afr<-unique(rbind(non_coordinate_sp2_afr,non_coordinate_sp3_afr, non_coordinate_mg2_afr, non_coordinate_mg3_afr))
total_vouchers_afr<-rbind(list_jabot_afr2,non_coordinate_total_afr)
colnames(total_vouchers_afr)
colnames(total_vouchers)
common_columns <- intersect(colnames(total_vouchers_afr), colnames(total_vouchers))
total_vouchers_afr2<-total_vouchers_afr %>% select(common_columns)
total_vouchers2<-total_vouchers %>% select(common_columns)
total_vouchers_jabot_gbif <- rbind(total_vouchers_afr2, total_vouchers2)



#Including data from FJBPC
list_pcs<-read.csv("/Users/xbaroc/Desktop/Radiation_project/list_pcs.csv", sep = ";")

list_pcs$species<-NA
for (i in 1:length(list_pcs$species)) {
  list_pcs$scientificName[i] <- paste(list_pcs$Genero[i], list_pcs$Especie[i], sep = " ")
}

total_vouchers_jabot_gbif2<-merge(total_vouchers_jabot_gbif, list_pcs, all.x=T, all.y=T)
#species_list<-unique(df[10])
#species_list<-unique(total_vouchers[10])
species_list<-unique(total_vouchers_jabot_gbif2[1])
species_list$scientificName<-sapply(species_list$scientificName, remove.authors)
species_list<-unique(species_list)
species_list_flora<- get.taxa(species_list$scientificName,life.form = T, habitat = T, vegetation.type = T,
                          vernacular = T, states = T, establishment = T,domain = T, endemism = T)




#Final species list
#campos<-subset(species_list_flora, grepl("Campo|Cerrado", species_list_flora$vegetation.type))
campos<-subset(species_list_flora, grepl("Campo", species_list_flora$vegetation.type))

campos_final<- unique(campos[-c(11,12)])
########################################################################################################################



#Final vouchers list

campos<-rename(campos, species = original.search)
total_vouchers_jabot_gbif<-rename(total_vouchers_jabot_gbif, species = scientificName)
total_vouchers_jabot_gbif$species_noauthor<-sapply(total_vouchers_jabot_gbif$species, remove.authors)
list_campos<-subset(total_vouchers_jabot_gbif, species_noauthor %in% campos$species)
campos_species<-campos[,c(2,12)]
campos_species<- campos_species %>%
  rename(species_noauthor = species)
list_campos2<-merge(list_campos, campos_species)
unique(list_campos2$scientific.name)


#Final map

gg4 <- ggplot() +
  geom_raster(data = elev2, aes(x = x, y = y, fill = value)) + 
  geom_point(data = list_campos2, aes(x = decimalLongitude, y = decimalLatitude), colour = "red", size = 0.3)+
  scale_fill_gradient(name = "Elevation (m)\n", low = "white", high = "darkblue", limits=c(500,2000)) +
  coord_equal() +labs(x = "Longitude", y = "Latitude")
gg4

################################################################



#Estimating richness
#Total

count_df <- list_campos2 %>%
  group_by(year, scientific.name) %>%
  summarize(count = n()) %>%
  ungroup()

# Reshape the data frame to create a presence/absence and species count matrix
library(dplyr)
library(tidyr)

pa_df <- count_df %>%
  pivot_wider(names_from = scientific.name, values_from = count, values_fill = 0) %>%
  mutate(Total = rowSums(.[, -which(names(.) == "year")]))
#pa_df <- count_df %>%
 # pivot_wider(names_from = scientific.name, values_from = count, values_fill = 0) %>%
  #mutate(Total = rowSums(select(., -year)))

rownames(pa_df) <- NULL
pa_df <- pa_df[,-1]

#get richness estimators (for each sample, cumulative)
poolaccum_tot <- poolaccum(pa_df)
specpool_tot <- specpool(pa_df)
specpool_tot

#plot all: obs richness and  estimators
plot(poolaccum_tot)

#build the species accumulation curve & rarefaction curve (expected)
specaccum_tot <- specaccum(pa_df,method = "rarefaction")

#plot the curve with some predefined settings
plot(specaccum_tot,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#build a expected curve (randomization for boxplot comparison)
specaccum_tot2 <- specaccum(pa_df, "random")
#plot both curves ("observed" vs "randomized")
plot(specaccum_tot,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(specaccum_tot2, col="yellow", add=TRUE, pch="+")

#Final plot
pdf('rarefaction.pdf',width = 12,height = 10)
plot(specaccum_tot,ci.type="poly", col="darkolivegreen", lwd=3, ci.lty=0, ci.col="darkolivegreen",  xlab="Years", ylab="Number of species")
boxplot(specaccum_tot2, add = TRUE, pch = "+", col = rgb(0, 0, 0, alpha = 0.2))
dev.off()


#### rever os NA do pacote flora reduzindo os is.na da porcentagem
#### Adicionar lista do edmilson
#### rever mais nomes de localidades que possam ser incluidos na cloud


#### Similarity angiosperms
KL<-read.csv("/Users/xbaroc/Documents/Artigos/Publicados/Campos de altitude_papers/CA_January_2023/AppS1_R1.csv")
KL <-subset(KL, grepl("Campo de Altitude", KL$vegetation.type))
KL2<- KL[c(2,16:30)]
str(KL)

CA<-campos_final
CA$pcs=1

values_to_exclude <- c(
  'Amblystegiaceae', 'Anemiaceae', 'Araucariaceae', 'Aspleniaceae', 'Aytoniaceae',
  'Bartramiaceae', 'Blechnaceae', 'Bryaceae', 'Cyatheaceae', 'Dennstaedtiaceae',
  'Dicranaceae', 'Dryopteridaceae', 'Dumortieraceae', 'Fissidentaceae',
  'Frullaniaceae', 'Funariaceae', 'Gleicheniaceae', 'Helicophyllaceae',
  'Heliotropiaceae', 'Hypnaceae', 'Isoetaceae', 'Jamesoniellaceae',
  'Lejeuneaceae', 'Lembophyllaceae', 'Lindsaeaceae', 'Lycopodiaceae',
  'Lygodiaceae', 'Marchantiaceae', 'Meteoriaceae', 'Metzgeriaceae',
  'Nephrolepidaceae', 'Ophioglossaceae', 'Osmundaceae', 'Plagiochilaceae',
  'Polypodiaceae', 'Polytrichaceae', 'Pteridaceae', 'Radulaceae',
  'Sphagnaceae', 'Thelypteridaceae', 'Thuidiaceae', 'Woodsiaceae'
)

# Use the filter() function from dplyr to exclude rows with specified values
library(dplyr)

CA2 <- CA %>%
  filter(!family %in% values_to_exclude)
CA3<-CA2[c(2,19)]
KLCA<-merge(KL2, CA3,all.x = T, all.y = T)
KLCA[is.na(KLCA)] <- 0

library(BinMat)
PUB_binary_matrix<-t(KLCA[,c(2:17)])

pdf('PUB_upgma_LSL2.pdf',width = 12,height = 10)
my_upgma<-BinMat::upgma(PUB_binary_matrix, method = 'binary', hclust ='mcquitty',  bts = 10000)
dev.off()

pdf('PUB_upgma_SSL2.pdf',width = 12,height = 10)
my_upgma<-BinMat::upgma(PUB_binary_matrix2, method = 'binary', hclust ='mcquitty',  bts = 10000)
dev.off()

cerrado <- read.csv("/Users/xbaroc/Desktop/Radiation_project/cerrado.csv", sep = ";")
values_to_exclude_cerrado <- c(
  'AMBLYSTEGIACEAE', 'ANEMIACEAE', 'ARAUCARIACEAE', 'ASPLENIACEAE', 'AYTONIACEAE',
  'BARTRAMIACEAE', 'BLECHNACEAE', 'BRYACEAE', 'CYATHEACEAE', 'DENNSTAEDTIACEAE',
  'DICRANACEAE', 'DRYOPTERIDACEAE', 'DUMORTIERACEAE', 'FISSIDENTACEAE',
  'FRULLANIACEAE', 'FUNARIACEAE', 'GLEICHENIACEAE', 'HELICOPHYLLACEAE',
  'HELIOTROPIACEAE', 'HYPNACEAE', 'ISOETACEAE', 'JAMESONIELLACEAE',
  'LEJEUNEACEAE', 'LEMBOPHYLLACEAE', 'LINDSAEACEAE', 'LYCOPODIACEAE',
  'LYGODIACEAE', 'MARCHANTIACEAE', 'METEORIACEAE', 'METZGERIACEAE',
  'NEPHROLEPIDACEAE', 'OPHIOGLOSSACEAE', 'OSMUNDACEAE', 'PLAGIOCHILACEAE',
  'POLYPODIACEAE', 'POLYTRICHACEAE', 'PTERIDACEAE', 'RADULACEAE',
  'SPHAGNACEAE', 'THELYPTERIDACEAE', 'THUIDIACEAE', 'WOODSIACEAE'
)
# Use the filter() function from dplyr to exclude rows with specified values
library(dplyr)
cerrado2 <- cerrado %>%
  filter(!FAMÍLIA %in% values_to_exclude_cerrado)
cerrado2<-cerrado2[c(1:24)]
columns_to_merge <- c(2, 3, 4, 6, 7)
cerrado2$species <- apply(cerrado2[, columns_to_merge], 1, function(row) paste(row, collapse = " "))
cerrado2$ITATIAIA<- NULL
cerrado2$ITATITOPO <- NULL
cerrado2$SERRA.NEGRA <- NULL
cerrado2$CAMPOS.SULINOS <- NULL
cerrado2$MACAEDECIMA <- NULL
cerrado2$CUCA <- NULL
cerrado2$BOCAINA <- NULL
cerrado2$ANTAS <- NULL
cerrado2$FRADE<- NULL
cerrado2$DESENGANO <- NULL
cerrado2$CERRADO.sl. <- NULL
your_data <- cerrado2[complete.cases(cerrado2[, c("CATOL", "GRMOG", "CIPNEW", "SSJ", "ALMAS")]), ]
cerrado2 <- cerrado2[!(is.na(cerrado2$CATOL) & is.na(cerrado2$GRMOG) & is.na(cerrado2$CIPNEW) & is.na(cerrado2$SSJ) & is.na(cerrado2$ALMAS)), ]
cerrado2<-unique(cerrado2)

library(flora)

cerrado_list_flora<- get.taxa(cerrado2$species,life.form = T, habitat = T, vegetation.type = T,
                              vernacular = T, states = T, establishment = T,domain = T, endemism = T)

cerrado2$original.search<-NULL
cerrado2$original.search=cerrado_list_flora$original.search
cerrado3<-cbind(cerrado2,cerrado_list_flora)
cerrado3<-cerrado3[c(9:13,17:35)]
cerrado3<-subset(cerrado3, grepl("Campo|Cerrado", cerrado3$vegetation.type))


KLCA2<-KLCA[,c(1, 15:29, 35)]
KLCA2[,c(2:17)] <- as.data.frame(lapply(KLCA2[,c(2:17)], function(x) as.numeric(as.character(x))))
str(KLCA2)
KLCA2[KLCA2 == 0] <- NA
cerrado4<-cerrado3[,c(6, 1:5)]
cerrado4[,c(2:6)]<- as.data.frame(lapply(cerrado4[,c(2:6)], function(x) as.numeric(as.character(x))))

KLCA_cerrado<-merge(KLCA2, cerrado4, all.x=T, all.y=T)
KLCA_cerrado[is.na(KLCA_cerrado)] <- 0
result <- KLCA_cerrado %>%
  group_by(scientific.name) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE))

result[, sapply(result, is.numeric)] <- lapply(result[, sapply(result, is.numeric)], function(x) ifelse(x > 1, 1, x))


library(BinMat)
PUB_binary_matrix<-t(result[,c(2:22)])

my_upgma<-BinMat::upgma(PUB_binary_matrix, method = 'binary', hclust ='mcquitty',  bts = 1000)

PUB_binary_matrix2 <- PUB_binary_matrix[!(row.names(PUB_binary_matrix) == "brigadeiro"), , drop = FALSE]

pdf('upgma.pdf',width = 12,height = 10)
my_upgma2<-BinMat::upgma(PUB_binary_matrix2, method = 'binary', hclust ='mcquitty',  bts = 10000)
dev.off()

#estimating numbers

library(flora)
library(abjutils)
library(dplyr)
library(ggplot2)
library(stringr)

CA$life.form<-rm_accent(CA$life.form)
CA$endemism<-rm_accent(CA$endemism)
CA$habitat<-rm_accent(CA$habitat)


cons_status<- CA %>%
  group_by(threat.status.mma2022) %>%
  count()
cons_status<-cons_status[order(cons_status$n, decreasing = T), ]

cons_status2<- CA %>%
  group_by(threat.status.cnc) %>%
  count()
cons_status2<-cons_status2[order(cons_status2$n, decreasing = T), ]

life_form<- CA %>%
  group_by(life.form) %>%
  count()
life_form<-life_form[order(life_form$n, decreasing = T), ]

habitat<- CA %>%
  group_by(habitat) %>%
  count()
habitat<-habitat[order(habitat$n, decreasing = T), ]

family<- CA %>%
  group_by(family) %>%
  count()
family<-family[order(family$n, decreasing = T), ]


data <- read.csv('/Users/xbaroc/Desktop/Radiation_project/PUB_merged_DATA_PCS.csv', sep= ';')


empty_bar <- 3
to_add <- data.frame(matrix(NA, empty_bar*4, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(unique(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))


label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)


base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))


grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-500,2500) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(0,6), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+40, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  
  geom_segment(data=base_data, aes(x = start, y = -30, xend = end, yend = -30), colour = "black", alpha=0.8, size=0.5 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -120, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p

pdf('PUB_merged_graphs.pdf',width = 10,height = 10)
p
dev.off()


### New map

install.packages("maptools", repos="http://R-Forge.R-project.org")
library(maptools)
library(raster)
library(elevatr)
e1 <- extent(c(-47,-46,-22.31,-21.49))
elev <- get_elev_raster(raster(e1, vals = 0), z = 12, prj = "+proj=longlat +datum=WGS84 +no_defs")

#srtm <- getData("SRTM", lon = -15.59972, lat = 27.965)
# crop to Gran Canaria & Tenerife

srtm_ct <- crop(elev, e1)
# plot slope and aspect with hill shades
slope <- terrain(srtm_ct, opt = "slope")
aspect <- terrain(srtm_ct, opt = "aspect")
hill <- hillShade(slope, aspect, angle = 45, direction = 45, normalize = TRUE)
pdf('map_new.pdf',width = 20,height = 20)
plot(hill, col = grey(0:100/100), legend = FALSE)
plot(srtm_ct, col = rainbow(25, alpha = 0.35), add = TRUE)
points(list_campos2$decimalLongitude, list_campos2$decimalLatitude, col = "black", pch = 16)
dev.off()


## Dissimilarity by geographic distance
library(vegan)
library(ggplot2)
library(jaccard)
library(tidyr)
geo<-read.csv("/Users/xbaroc/Desktop/Radiation_project/location.csv", sep = ";")
geo <- separate(geo, col = "Reference.coordinates..S..W.", into = c("decimalLatitude", "decimalLongitude"), sep = ",")
library(stringr)
geo$Mountaintop.area <- str_extract(geo$Mountaintop.area, "\\((.*?)\\)")
geo$Mountaintop.area <- gsub("[()]", "", geo$Mountaintop.area)

PUB_binary_matrix3<-PUB_binary_matrix2
rownames(PUB_binary_matrix3) <-c("PNIT","PNSO","PNSJ","PNCP","PECJ","PNSB","PEPP","ARAC","PETP","APAP","PEDS","PIMA","PESP","PDMI","PCSP","CATO","GRMG","SECI","SESJ","PIAM")
PUB_JAC<-vegdist(PUB_binary_matrix3, method = "jaccard")
PUB_mydf<-as.data.frame(as.matrix(PUB_JAC))
write.csv2(PUB_mydf,file = "dissimilarity.csv")
write.csv(geo, file = "geo.csv")

library(fossil)
frame=geo[,c(3,2,1)] #Ordem: Longitude, Latitude, Especie
rownames(frame)<-frame$Mountaintop.area
frame$decimalLongitude<- gsub("[^0-9.]", "", frame$decimalLongitude)
frame$decimalLongitude<-as.numeric(frame$decimalLongitude)
frame$decimalLongitude<- -abs(frame$decimalLongitude)
frame$decimalLatitude<- gsub("[^0-9.]", "", frame$decimalLatitude)
frame$decimalLatitude<-as.numeric(frame$decimalLatitude)
frame$decimalLatitude<- -abs(frame$decimalLatitude)
shp=earth.dist(frame, dist = T)
class(shp)
shp<-as.matrix(shp)
row.names(shp) <- rownames(frame)
colnames(shp) <- rownames(frame)
shp2<-as.data.frame(shp)

combinations1 <- combn(rownames(PUB_mydf), 2, simplify = TRUE)
combinations1 <- apply(combinations1, 2, function(x) paste(sort(x), collapse = " vs "))

df <- data.frame(locations = combinations1,
                 dissimilarity = PUB_mydf[lower.tri(PUB_mydf)])

# For shp2
combinations2 <- combn(rownames(shp2), 2, simplify = TRUE)
combinations2 <- apply(combinations2, 2, function(x) paste(sort(x), collapse = " vs "))

df2 <- data.frame(locations = combinations2,
                  geo = shp2[lower.tri(shp2)])

# Merge the dataframes
df3 <- merge(df, df2, by = "locations")
df3$log_geo<-log10(df3$geo)
df3$log_dis<-log10(df3$dissimilarity)

library(caper)
library(ggplot2)
library(ggpubr)
library(lsmeans)
model1<-glm(dissimilarity ~ geo, data = df3)
summary(model1)
qqnorm(residuals(model1, type="deviance"));qqline(residuals(model1, type="deviance"))
plot(predict(model1), residuals(model1))
shapiro.test(model1$residuals)
hist(residuals(model1))

model2<-glm(log_dis~ log_geo, data = df3)
summary(model2)
qqnorm(residuals(model2, type="deviance"));qqline(residuals(model2, type="deviance"))
plot(predict(model2), residuals(model2))
shapiro.test(model2$residuals)


ggscatter(df3, x = "geo", y = "dissimilarity", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = 'geo', ylab = "dissimilarity", alpha=0.1)

merged_matrix <- rbind(shp2, PUB_mydf)

#Removing accents

if (!requireNamespace("stringi", quietly = TRUE)) {
  install.packages("stringi")
}
library(stringi)

# Assuming campos_final is your dataframe
# Loop through each column and remove accents from the values
campos_final2<-campos_final
for (col in names(campos_final2)) {
  campos_final2[[col]] <- stri_trans_general(campos_final2[[col]], "Latin-ASCII")
}

write.csv2(campos_final2, file = "list_species_pcs.csv")
