#Ryan Williams 1-5-13
#make sure you have these packages loaded
getwd() 
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

####check to see if you can get data in 


#long lists can be made quickly in excel from the MG-RAST Table.xls that is on the Hofmockel_Kbase Google Drive.  By concatenating different cells together you can generate lists that can be pasted directly into R

#also below in data.metagenome is all the WGS sequences
#"CC12-Micro-July2012"='mgm4511174.3' ,
#"CC21-Micro-July2012"='mgm4511175.3' ,
#"CC35-Micro-July2012"='mgm4511176.3' ,
#"CC43-Micro-July2012"='mgm4511173.3' ,


data.metagenome<-c(
"PF15-LM-July2012"='mgm4509402.3' ,"PF23-LM-July2012"='mgm4511167.3' ,"PF32-LM-July2012"='mgm4509406.3' ,"PF41-LM-July2012"='mgm4509400.3' ,"PF15-Micro-July2012"='mgm4509398.3' ,"PF23-Micro-July2012"='mgm4509405.3' ,"PF32-Micro-July2012"='mgm4509401.3' ,"PF41-Micro-July2012"='mgm4509399.3' ,"PF15-MM-July2012"='mgm4511171.3' ,"PF23-MM-July2012"='mgm4511168.3' ,"PF32-MM-July2012"='mgm4511170.3' ,"PF41-MM-July2012"='mgm4511169.3' ,"PF15-SM-July2012"='mgm4509396.3' ,"PF23-SM-July2012"='mgm4511172.3' ,"PF32-SM-July2012"='mgm4509397.3' ,"PF41-SM-July2012"='mgm4509404.3' ,"PF15-WS-July2012"='mgm4509403.3' ,"PF23-WS-July2012"='mgm4509407.3' ,"PF32-WS-July2012"='mgm4511177.3' ,"PF41-WS-July2012"='mgm4511166.3' 
)

#make your collection #m5rna is nucleotide database
#community.2<-collection(data.2, count = c (entry = "count", annot = "organism", source = "m5rna", level = "class"))
#entry="count" is count data, "normed" is normalized,standardized & scaled by total count
#level3 - level= "level3"  level 4 is level="function"; level 4 tierd KW, followed by t test with top candidates, for example; need some level of false discovery rate control for statisticians to be happy.
my_data <- collection(data.metagenome, count = c(entry = "normed", annot = "func", source = "Subsystem", level = "level3"))
# remove singletons
#community.singletons_removed <- remove_singletons(as.matrix(community.2$count))

# see ordering of columns in imported data
colnames(my_data)
colnames(my_data$count)
#see first few rows
 head(my_data$count)
# say you want cols in order 5,1,2,3,7,4,6
#my_data$count [ , c (5,1,2,3,7,4,6)]
# also works with metagenome ids:
#my_data$count [ , c (“4441166.3”, “444....etc...
# but may sometimes need “mgm” prefix..
organized_data = my_data$count [ , data.metagenome]
analysis_1.kw_stat <-sigtest(organized_data,c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4)), "Kruskal-Wallis")

write.csv(analysis_1.kw_stat$stat,"KWwgs.csv")

# create a data opject that has the data used in the statistical analysis as the stat output
data_out <- cbind(my_data$count, analysis_1.kw_stat$stat)

#test=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4))
#test=c(rep(1,4))
#test
#To cite this package in a publication, start R and enter:
#    citation("qvalue")

#this allows us to source into R from url

#source_https <- function(url, ...) {
	     # load package
#	     require(RCurl)

	     # parse and evaluate each .R script
#	     sapply(c(url, ...), function(u) {
#	     eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
#	     })
#	}

#source_https("https://raw.github.com/braithwaite/matR-apps/master/qvalue-function.R")




#command arrow & shift command arrow in excel to highlight column
my_p<-scan()
#highlight a column from excel, paste it, double enter
#fdr level is significance level - false discovery rate
my_q<-qvalue(my_p,fdr.level=0.001)

summary(my_q)

#make qwrite
qwrite <<- function (qobj, filename = "my-qvalue-results.txt") {    cat(c("pi0:", qobj$pi0, "\n\n"), file = filename, append = FALSE)    if (any(names(qobj) == "fdr.level")) {        cat(c("FDR level:", qobj$fdr.level, "\n\n"), file = filename,             append = TRUE)        cat(c("p-value q-value significant", "\n"), file = filename,             append = TRUE)        write(t(cbind(qobj$pval, qobj$qval, qobj$significant)),             file = filename, ncolumns = 3, append = TRUE)    }    else {        cat(c("p-value q-value", "\n"), file = filename, append = TRUE)        write(t(cbind(qobj$pval, qobj$qval)), file = filename,             ncolumns = 2, append = TRUE)    }}
qwrite(my_q,filename="myqresults")

mydataobject<-cbind(my_data$count,analysis_1.kw_stat$stat,my_q$qvalues, my_q$significant)

write.table(mydataobject, file = "analysis_1.kw_stat.txt", col.names=NA, row.names = TRUE, sep="\t", quote=FALSE)

#subset my_data for just micros & large macros
LM_micro<-c("PF15-LM-July2012"='mgm4509402.3' ,
"PF23-LM-July2012"='mgm4511167.3' ,
"PF32-LM-July2012"='mgm4509406.3' ,
"PF41-LM-July2012"='mgm4509400.3',
"PF15-Micro-July2012"='mgm4509398.3' ,
"PF23-Micro-July2012"='mgm4509405.3' ,
"PF32-Micro-July2012"='mgm4509401.3' ,
"PF41-Micro-July2012"='mgm4509399.3')

my_subset<-my_data[LM_micro]
colnames(my_data$count)
organized_subset = my_subset$count [ , LM_micro]
colnames(organized_subset)

analysis_2.kw_stat <-sigtest(organized_subset,c(rep(1,4),rep(2,4)), "t-test-un-paired")

write.csv(cbind(my_subset$count,analysis_2.kw_stat$stat),"ttestLMmicro.csv")

# create a data opject that has the data used in the statistical analysis as the stat output
##

