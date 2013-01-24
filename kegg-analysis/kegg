#need these libraries
library(ggplot2)
library(reshape)


#read in the data from wherever you have saved it on your computer
data<-read.csv(file.choose())

#I had assemblies wtih all 4 levels with all 4 aggregate size classes (Micro, SM, MM, LM)
#The melt function will keep all columns specified in the id=() the same.  The other columns will be condensed into two columns.  The column labeled variable will contain all names of the condensed columns.  The column labeled value will contain all values that were held in the columns that were condensed.  This function can be confusing at first so it may be useful to look at your data before and after the melt
data.melt<-melt(data, id=c("Level.1","Level.2","Level.3","Level.4"))

names(data.melt)
head(data.melt)

#made columns labeled Fraction and Abundance to hold the data made from the melt function and then removed the variable and value column
data.melt$Fraction<-data.melt$variable
data.melt$Abundance<-data.melt$value
names(data.melt)
data.melt<-data.melt[,c(1:4,7,8)]
names(data.melt)


#the ggplot function works in this manner: ggplot(dataset you are plotting from, aes(x,y))+geom_bar(aes(fill=color bars based on a column in the data set)).  coord_flip() rotates the figure so the x-axis is in place of the y-axis.  This rotates the text as well and is good for these kinds of figures.  theme_bw() makes the figure a little more clean looking and opts(aspect.ratio=1) keeps a square aspect ratio.  You may need to enlarge the window for your imag (green + sign) to make it fit better. 
ggplot(data.melt, aes(Level.1, Abundance))+geom_bar(aes(fill=Fraction),stat="identity")+coord_flip()+theme_bw()+opts(aspect.ratio=1)


