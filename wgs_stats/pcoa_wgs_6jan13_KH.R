#Kirsten Hofmockel 1-6-13
#make sure you have these packages loaded
library(matR)
library(RJSONIO)
library(reshape)
library(vegan)
library(ggplot2)
library(RCurl)

#This is all Kevin Keegan's code for functions that he wrote#

source_https <- function(url, ...) {
# load package
require(RCurl)

# parse and evaluate each .R script
sapply(c(url, ...), function(u) {
eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
})
}

source_https("https://raw.github.com/DrOppenheimer/matR-apps/master/remove_singletons.r",
"https://raw.github.com/DrOppenheimer/matR-apps/master/collection-merge.R",
"https://raw.github.com/DrOppenheimer/matR-apps/master/norm_center_scale.r",
"https://raw.github.com/DrOppenheimer/matR-apps/master/plot_pcoa.from_object.r",
"https://raw.github.com/DrOppenheimer/matR-apps/master/heatmap_dendrogram.from_file.r")
#Now you can use Kevin's functions from the VM on your computer in R#


####do the first to using a new key generated from MG-RAST

mconfig$setAuth()

Micro_July.list <- c(
"CC12-Micro-July2012"="4506606.3",
"CC21-Micro-July2012"="4506616.3",
"CC35-Micro-July2012"="4506626.3",
"CC43-Micro-July2012"="4506636.3",
"P13-Micro-July2012"="4506646.3",
"P24-Micro-July2012"="4506656.3",
"P31-Micro-July2012"="4506666.3",
"P46-Micro-July2012"="4506676.3",
"PF15-Micro-July2012"="4506686.3",
"PF23-Micro-July2012"="4506696.3",
"PF32-Micro-July2012"="4506706.3",
"PF41-Micro-July2012"="4506716.3"
)

#to make hard copy
write.csv(Micro_July.list,"Micro_July.list.txt")

Micro_July.colors <- c(
"CC12-Micro-July2012"="red",
"CC21-Micro-July2012"="red",
"CC35-Micro-July2012"="red",
"CC43-Micro-July2012"="red",
"P13-Micro-July2012"="blue",
"P24-Micro-July2012"="blue",
"P31-Micro-July2012"="blue",
"P46-Micro-July2012"="blue",
"PF15-Micro-July2012"="green",
"PF23-Micro-July2012"="green",
"PF32-Micro-July2012"="green",
"PF41-Micro-July2012"="green"
)

# get the raw abundance data

Micro_July.class.count_data <- collection(Micro_July.list, count = c (entry = "count", annot = "organism", source = "M5RNA", level = "class"))

# remove singletons
Micro_July.class.count_data.singletons_removed <- remove_singletons(as.matrix(Micro_July.class.count_data$count))

# normalize, standardize, and scale
Micro_July.class.count_data.singletons_removed.normed <- norm_center_scale(Micro_July.class.count_data.singletons_removed)

# save processed CLASS data to file - with the sample names as the row.names(class)
Micro_July.class.count_data.singletons_removed.normed.NAMES <- Micro_July.class.count_data.singletons_removed.normed
dimnames(Micro_July.class.count_data.singletons_removed.normed.NAMES)[[2]] <- names(Micro_July.list) # here mgids are replaced with sample names
write.table(Micro_July.class.count_data.singletons_removed.normed.NAMES, file = "Micro_July.class.count_data.singletons_removed.normed.NAMES.12-13-12.txt", row.names=TRUE, col.names=NA, sep="\t", quote=FALSE, append=TRUE)

# create boxplots of raw and processed (singletons removed, norm-center-scale) data
split.screen(c(2,1))
screen(1)
boxplot(Micro_July.data$count, main = "raw_data", las=2, col=Micro_July.colors, names = names(Micro_July.data))
screen(2)
boxplot(Micro_July.data.singletons_removed.normed, main = "singletons_removed.norm-cent-scaled", las=2, col=Micro_July.colors, names = names(Micro_July.data))

# create pcoa
plot_pcoa.from_object(matrix=Micro_July.data.singletons_removed.normed, colors=Micro_July.colors)

# move to new terminal & run the following code:
# scp hofmockel.vm:/mnt/16s.by_fraction_and_month.12-10-12/
# matrix_object.euclidean.PCoA.png ./
# open matrix_object.euclidean.PCoA.png

#plot_pco_with_stats.9-18-12.pl --data_file #class.count_data.singletons_removed.normed.NAMES.12-13-12.txt --groups_list #july_oct.groups --num_perm 1000 --dist_method euclidean --cleanup 2>error.log &

#Linnux:
#plot_pco_with_stats.9-18-12.pl --data_file Micro_July.class.count_data.singletons_removed.normed.NAMES.12-13-12.txt --groups_list groups.list.new.txt --num_perm 1000 --dist_method euclidean --cleanup 2>error.log &