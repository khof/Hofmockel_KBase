# THis is a more generalized version of 10-15-12 -- created specifically to plot the EHFI data -- this does some more automated stuff
#
# This is a script to plot the PCoA outputs of the following functions:
#
#     plot_pco_with_stats,
#     plot_qiime_pco_with_stats, or
#     plot_OTU_pco_with_stats
#
#     Relevent output from these scripts *.PCoA
#
# These can be driven by the following master function
#     plot_pco_with_stats.MASTER
# May be nice to add the plotting as a standard output to this batch processing script
# and or the same functionality to each of the 3 scripts above
#


plot_pcoa <- function(
                      file_in,
                      #groups_in,
                      out_prefix = "Analysis_1.",
                      my_colors = "red",
                      legend_colors <- c ('red'),
                      names (legend_colors) <- c ('group_1'),
                      figure_width = 500,
                      figure_height = 500,
                      figure_res = NA,
                      debug = FALSE,
                      PC1 = 1,
                      PC2 = 2
                      )

  {

    # read file_in into a list
    con <- file(file_in)
    num_lines <- length(readLines(con))
    if(debug==TRUE){print(paste("num_lines:", num_lines))}
    open(con)
    results.list <<- list();
    current.line <- 1
    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
      results.list[[current.line]] <<- line
      current.line <- current.line + 1
    } 
    close(con)    

    # parse file_in into Eigen values and Eigen vectors

    num_eigen_values <- length(grep("PC", results.list)) # need to figure out how to do emty line regex in R ...

    eigen_value.list <<- vector("list", num_eigen_values)
    eigen_vector.list <<- vector("list", num_eigen_values)
    value.count <- 1
    vector.count <- 1

    for (i in 1:num_lines){
      my_line <<- noquote(results.list[i])
      if (any(grep("^#", my_line))==TRUE){
        # skip comment line
        
      }else{
      #if(any(my_line)==TRUE){

        #if(debug==TRUE){print(paste(i, my_line))}
        if (any(grep("PC", my_line))==TRUE){ # pull out the eigen value
          #if(debug==TRUE){print("made it here (3)")}
          eigen_value.list[value.count] <<- my_line 
          value.count <- value.count + 1
          if(debug==TRUE){print(paste("PC", my_line))}
        }else if (any(grep("[0-9]", my_line))==TRUE){ # pull out the eigen vectors
          eigen_vector.list[vector.count] <<- my_line
          vector.count <- vector.count + 1
          #if(debug==TRUE){print(paste("vector", line))}
        }else{print("encountered empty variable")}
        
      }
    }

    test_something <<- strsplit(as.character(eigen_vector.list), split="\t")

    my_names <<- vector("list", num_eigen_values)
    for (i in 1:num_eigen_values){
      my_names[i] <<- test_something[[i]][1]
    }

    my_data <<- matrix(ncol = num_eigen_values, nrow = num_eigen_values) 
    for (i in 1:num_eigen_values){
      for (j in 1:num_eigen_values){
        my_data[i,j] <<- test_something[[i]][j+1]
      }
    }

    rownames(my_data) <<- my_names

    eigen_value.names <- vector("list", num_eigen_values)
    for (i in 1:num_eigen_values){
      pcName_pcEigen <- unlist(strsplit(as.character(eigen_value.list[i]), split="\t"))
      pcName <- gsub("\"", "", (pcName_pcEigen[1]))
      pcEigen <- gsub(" ", "", paste(  (round( as.numeric(pcName_pcEigen[2]), digits = 4))*100, "%"))
      eigen_value.names[i] <- paste( pcName, "::", pcEigen, "of observed variation")
    } 

    sorted_my_data <<- my_data[order(rownames(my_data)), ]
    colnames(sorted_my_data) <<- eigen_value.names


    #write.table(sorted_my_data, file = "sorted_names.txt", row.names = FALSE, sep="\t", quote=FALSE)
    write.table(sorted_my_data, file = "sorted_names.txt", row.names = FALSE, sep="\t", quote=FALSE)
    #write.table(analysis.2.SM.raw.data.sans_sg.normed, file = "analysis.2.SM.data.sans_sg.normed.txt", col.names=NA, row.names = TRUE, sep="\t", quote=FALSE)
    
    min_PC1 <- min(sorted_my_data[,PC1])
    max_PC1 <- max(sorted_my_data[,PC1])
    min_PC2 <- min(sorted_my_data[,PC2])
    max_PC2 <- max(sorted_my_data[,PC2])


    # Now get the group information out of the groups_list and use that to generate automatic colors


    ###
    
    image_out = gsub(" ", "", (paste(out_prefix, ".plot.pdf")))
    legend_out = gsub(" ", "", (paste(out_prefix, ".legend.pdf")))
    
    pdf(file=image_out, width = 10, height = 10)
    plot(
         x<-sorted_my_data[,PC1],
         y<-sorted_my_data[,PC2],        
         type="n",
         xlab = eigen_value.names[PC1],
         ylab = eigen_value.names[PC2],
         cex = 0.8
         )
    my_cex <- 1.2

    # analysis 1
    ## my_colors <<- c (
    ##                  "CC12-LM-July2012"="red",
    ##                  "CC12-LM-October2012"="blue",
    ##                  "CC12-Micro-July2012"="red",
    ##                  "CC12-Micro-October2012"="blue",
    ##                  "CC12-MM-July2012"="red",
    ##                  "CC12-MM-October2012"="blue",
    ##                  "CC12-SM-July2012"="red",
    ##                  "CC12-SM-October2012"="blue",
    ##                  "CC12-WS-July2012"="red",
    ##                  "CC12-WS-October2012"="blue",
    ##                  "CC21-LM-July2012"="red",
    ##                  "CC21-LM-October2012"="blue",
    ##                  "CC21-Micro-July2012"="red",
    ##                  "CC21-Micro-October2012"="blue",
    ##                  "CC21-MM-July2012"="red",
    ##                  "CC21-MM-October2012"="blue",
    ##                  "CC21-SM-July2012"="red",
    ##                  "CC21-SM-October2012"="blue",
    ##                  "CC21-WS-July2012"="red",
    ##                  "CC21-WS-October2012"="blue",
    ##                  "CC35-LM-July2012"="red",
    ##                  "CC35-LM-October2012"="blue",
    ##                  "CC35-Micro-July2012"="red",
    ##                  "CC35-Micro-October2012"="blue",
    ##                  "CC35-MM-July2012"="red",
    ##                  "CC35-MM-October2012"="blue",
    ##                  "CC35-SM-July2012"="red",
    ##                  "CC35-SM-October2012"="blue",
    ##                  "CC35-WS-July2012"="red",
    ##                  "CC35-WS-October2012"="blue",
    ##                  "CC43-LM-July2012"="red",
    ##                  "CC43-LM-October2012"="blue",
    ##                  "CC43-Micro-July2012"="red",
    ##                  "CC43-Micro-October2012"="blue",
    ##                  "CC43-MM-July2012"="red",
    ##                  "CC43-MM-October2012"="blue",
    ##                  "CC43-SM-July2012"="red",
    ##                  "CC43-SM-October2012"="blue",
    ##                  "CC43-WS-July2012"="red",
    ##                  "CC43-WS-October2012"="blue",
    ##                  "P13-LM-July2012"="green",
    ##                  "P13-LM-October2012"="magenta",
    ##                  "P13-Micro-July2012"="green",
    ##                  "P13-Micro-October2012"="magenta",
    ##                  "P13-MM-July2012"="green",
    ##                  "P13-MM-October2012"="magenta",
    ##                  "P13-SM-July2012"="green",
    ##                  "P13-SM-October2012"="magenta",
    ##                  "P13-WS-July2012"="green",
    ##                  "P13-WS-October2012"="magenta",
    ##                  "P24-LM-July2012"="green",
    ##                  "P24-LM-October2012"="magenta",
    ##                  "P24-Micro-July2012"="green",
    ##                  "P24-Micro-October2012"="magenta",
    ##                  "P24-MM-July2012"="green",
    ##                  "P24-MM-October2012"="magenta",
    ##                  "P24-SM-July2012"="green",
    ##                  "P24-SM-October2012"="magenta",
    ##                  "P24-WS-July2012"="green",
    ##                  "P24-WS-October2012"="magenta",
    ##                  "P31-LM-July2012"="green",
    ##                  "P31-LM-October2012"="magenta",
    ##                  "P31-Micro-July2012"="green",
    ##                  "P31-Micro-October2012"="magenta",
    ##                  "P31-MM-July2012"="green",
    ##                  "P31-MM-October2012"="magenta",
    ##                  "P31-SM-July2012"="green",
    ##                  "P31-SM-October2012"="magenta",
    ##                  "P31-WS-July2012"="green",
    ##                  "P31-WS-October2012"="magenta",
    ##                  "P46-LM-July2012"="green",
    ##                  "P46-LM-October2012"="magenta",
    ##                  "P46-Micro-July2012"="green",
    ##                  "P46-Micro-October2012"="magenta",
    ##                  "P46-MM-July2012"="green",
    ##                  "P46-MM-October2012"="magenta",
    ##                  "P46-SM-July2012"="green",
    ##                  "P46-SM-October2012"="magenta",
    ##                  "P46-WS-July2012"="green",
    ##                  "P46-WS-October2012"="magenta",
    ##                  "PF15-LM-July2012"="orange",
    ##                  "PF15-LM-October2012"="purple",
    ##                  "PF15-Micro-July2012"="orange",
    ##                  "PF15-Micro-October2012"="purple",
    ##                  "PF15-MM-July2012"="orange",
    ##                  "PF15-MM-October2012"="purple",
    ##                  "PF15-SM-July2012"="orange",
    ##                  "PF15-SM-October2012"="purple",
    ##                  "PF15-WS-July2012"="orange",
    ##                  "PF15-WS-October2012"="purple",
    ##                  "PF23-LM-July2012"="orange",
    ##                  "PF23-LM-October2012"="purple",
    ##                  "PF23-Micro-July2012"="orange",
    ##                  "PF23-Micro-October2012"="purple",
    ##                  "PF23-MM-July2012"="orange",
    ##                  "PF23-MM-October2012"="purple",
    ##                  "PF23-SM-July2012"="orange",
    ##                  "PF23-SM-October2012"="purple",
    ##                  "PF23-WS-July2012"="orange",
    ##                  "PF23-WS-October2012"="purple",
    ##                  "PF32-LM-July2012"="orange",
    ##                  "PF32-LM-October2012"="purple",
    ##                  "PF32-Micro-July2012"="orange",
    ##                  "PF32-Micro-October2012"="purple",
    ##                  "PF32-MM-July2012"="orange",
    ##                  "PF32-MM-October2012"="purple",
    ##                  "PF32-SM-July2012"="orange",
    ##                  "PF32-SM-October2012"="purple",
    ##                  "PF32-WS-July2012"="orange",
    ##                  "PF32-WS-October2012"="purple",
    ##                  "PF41-LM-July2012"="orange",
    ##                  "PF41-LM-October2012"="purple",
    ##                  "PF41-Micro-July2012"="orange",
    ##                  "PF41-Micro-October2012"="purple",
    ##                  "PF41-MM-July2012"="orange",
    ##                  "PF41-MM-October2012"="purple",
    ##                  "PF41-SM-July2012"="orange",
    ##                  "PF41-SM-October2012"="purple",
    ##                  "PF41-WS-July2012"="orange",
    ##                  "PF41-WS-October2012"="purple"
    ##                  )

    ## # analysis 2.WS
    ## my_colors <<- c (
    ##                  "CC12.WS.July2012"="red",
    ##                  "CC12.WS.October2012"="red",
    ##                  "CC21.WS.July2012"="red",
    ##                  "CC21.WS.October2012"="red",
    ##                  "CC35.WS.July2012"="red",
    ##                  "CC35.WS.October2012"="red",
    ##                  "CC43.WS.July2012"="red",
    ##                  "CC43.WS.October2012"="red",
    ##                  "P13.WS.July2012"="blue",
    ##                  "P13.WS.October2012"="blue",
    ##                  "P24.WS.July2012"="blue",
    ##                  "P24.WS.October2012"="blue",
    ##                  "P31.WS.July2012"="blue",
    ##                  "P31.WS.October2012"="blue",
    ##                  "P46.WS.July2012"="blue",
    ##                  "P46.WS.October2012"="blue",
    ##                  "PF15.WS.July2012"="green",
    ##                  "PF15.WS.October2012"="green",
    ##                  "PF23.WS.July2012"="green",
    ##                  "PF23.WS.October2012"="green",
    ##                  "PF32.WS.July2012"="green",
    ##                  "PF32.WS.October2012"="green",
    ##                  "PF41.WS.July2012"="green",
    ##                  "PF41.WS.October2012"="green"
    ##                  )

    ## # analysis 2.LM
    ## my_colors <<- c (
    ##                  "CC12.LM.July2012"="red",
    ##                  "CC12.LM.October2012"="red",
    ##                  "CC21.LM.July2012"="red",
    ##                  "CC21.LM.October2012"="red",
    ##                  "CC35.LM.July2012"="red",
    ##                  "CC35.LM.October2012"="red",
    ##                  "CC43.LM.July2012"="red",
    ##                  "CC43.LM.October2012"="red",
    ##                  "P13.LM.July2012"="blue",
    ##                  "P13.LM.October2012"="blue",
    ##                  "P24.LM.July2012"="blue",
    ##                  "P24.LM.October2012"="blue",
    ##                  "P31.LM.July2012"="blue",
    ##                  "P31.LM.October2012"="blue",
    ##                  "P46.LM.July2012"="blue",
    ##                  "P46.LM.October2012"="blue",
    ##                  "PF15.LM.July2012"="green",
    ##                  "PF15.LM.October2012"="green",
    ##                  "PF23.LM.July2012"="green",
    ##                  "PF23.LM.October2012"="green",
    ##                  "PF32.LM.July2012"="green",
    ##                  "PF32.LM.October2012"="green",
    ##                  "PF41.LM.July2012"="green",
    ##                  "PF41.LM.October2012"="green"
    ##                  )

    ## # analysis 2.MM
    ## my_colors <<- c (
    ##                  "CC12.MM.July2012"="red",
    ##                  "CC12.MM.October2012"="red",
    ##                  "CC21.MM.July2012"="red",
    ##                  "CC21.MM.October2012"="red",
    ##                  "CC35.MM.July2012"="red",
    ##                  "CC35.MM.October2012"="red",
    ##                  "CC43.MM.July2012"="red",
    ##                  "CC43.MM.October2012"="red",
    ##                  "P13.MM.July2012"="blue",
    ##                  "P13.MM.October2012"="blue",
    ##                  "P24.MM.July2012"="blue",
    ##                  "P24.MM.October2012"="blue",
    ##                  "P31.MM.July2012"="blue",
    ##                  "P31.MM.October2012"="blue",
    ##                  "P46.MM.July2012"="blue",
    ##                  "P46.MM.October2012"="blue",
    ##                  "PF15.MM.July2012"="green",
    ##                  "PF15.MM.October2012"="green",
    ##                  "PF23.MM.July2012"="green",
    ##                  "PF23.MM.October2012"="green",
    ##                  "PF32.MM.July2012"="green",
    ##                  "PF32.MM.October2012"="green",
    ##                  "PF41.MM.July2012"="green",
    ##                  "PF41.MM.October2012"="green"
    ##                  )

    ## # analysis 2.SM
    ## my_colors <<- c (
    ##                  "CC12.SM.July2012"="red",
    ##                  "CC12.SM.October2012"="red",
    ##                  "CC21.SM.July2012"="red",
    ##                  "CC21.SM.October2012"="red",
    ##                  "CC35.SM.July2012"="red",
    ##                  "CC35.SM.October2012"="red",
    ##                  "CC43.SM.July2012"="red",
    ##                  "CC43.SM.October2012"="red",
    ##                  "P13.SM.July2012"="blue",
    ##                  "P13.SM.October2012"="blue",
    ##                  "P24.SM.July2012"="blue",
    ##                  "P24.SM.October2012"="blue",
    ##                  "P31.SM.July2012"="blue",
    ##                  "P31.SM.October2012"="blue",
    ##                  "P46.SM.July2012"="blue",
    ##                  "P46.SM.October2012"="blue",
    ##                  "PF15.SM.July2012"="green",
    ##                  "PF15.SM.October2012"="green",
    ##                  "PF23.SM.July2012"="green",
    ##                  "PF23.SM.October2012"="green",
    ##                  "PF32.SM.July2012"="green",
    ##                  "PF32.SM.October2012"="green",
    ##                  "PF41.SM.July2012"="green",
    ##                  "PF41.SM.October2012"="green"
    ##                  )

    ## # analysis 2.SM
    ## my_colors <<- c (
    ##                  "CC12.Micro.July2012"="red",
    ##                  "CC12.Micro.October2012"="red",
    ##                  "CC21.Micro.July2012"="red",
    ##                  "CC21.Micro.October2012"="red",
    ##                  "CC35.Micro.July2012"="red",
    ##                  "CC35.Micro.October2012"="red",
    ##                  "CC43.Micro.July2012"="red",
    ##                  "CC43.Micro.October2012"="red",
    ##                  "P13.Micro.July2012"="blue",
    ##                  "P13.Micro.October2012"="blue",
    ##                  "P24.Micro.July2012"="blue",
    ##                  "P24.Micro.October2012"="blue",
    ##                  "P31.Micro.July2012"="blue",
    ##                  "P31.Micro.October2012"="blue",
    ##                  "P46.Micro.July2012"="blue",
    ##                  "P46.Micro.October2012"="blue",
    ##                  "PF15.Micro.July2012"="green",
    ##                  "PF15.Micro.October2012"="green",
    ##                  "PF23.Micro.July2012"="green",
    ##                  "PF23.Micro.October2012"="green",
    ##                  "PF32.Micro.July2012"="green",
    ##                  "PF32.Micro.October2012"="green",
    ##                  "PF41.Micro.July2012"="green",
    ##                  "PF41.Micro.October2012"="green"
    ##                  )

    ## # analysis 3.CC
    ## my_colors <<- c (
    ##                  "CC12.WS.July2012"="orange",
    ##                  "CC12.WS.October2012"="orange",
    ##                  "CC21.WS.July2012"="orange",
    ##                  "CC21.WS.October2012"="orange",
    ##                  "CC35.WS.July2012"="orange",
    ##                  "CC35.WS.October2012"="orange",
    ##                  "CC43.WS.July2012"="orange",
    ##                  "CC43.WS.October2012"="orange",
    ##                  "CC12.LM.July2012"="red",
    ##                  "CC12.LM.October2012"="red",
    ##                  "CC21.LM.July2012"="red",
    ##                  "CC21.LM.October2012"="red",
    ##                  "CC35.LM.July2012"="red",
    ##                  "CC35.LM.October2012"="red",
    ##                  "CC43.LM.July2012"="red",
    ##                  "CC43.LM.October2012"="red",
    ##                  "CC12.MM.July2012"="blue",
    ##                  "CC12.MM.October2012"="blue",
    ##                  "CC21.MM.July2012"="blue",
    ##                  "CC21.MM.October2012"="blue",
    ##                  "CC35.MM.July2012"="blue",
    ##                  "CC35.MM.October2012"="blue",
    ##                  "CC43.MM.July2012"="blue",
    ##                  "CC43.MM.October2012"="blue",
    ##                  "CC12.SM.July2012"="magenta",
    ##                  "CC12.SM.October2012"="magenta",
    ##                  "CC21.SM.July2012"="magenta",
    ##                  "CC21.SM.October2012"="magenta",
    ##                  "CC35.SM.July2012"="magenta",
    ##                  "CC35.SM.October2012"="magenta",
    ##                  "CC43.SM.July2012"="magenta",
    ##                  "CC43.SM.October2012"="magenta",
    ##                  "CC12.Micro.July2012"="green",
    ##                  "CC12.Micro.October2012"="green",
    ##                  "CC21.Micro.July2012"="green",
    ##                  "CC21.Micro.October2012"="green",
    ##                  "CC35.Micro.July2012"="green",
    ##                  "CC35.Micro.October2012"="green",
    ##                  "CC43.Micro.July2012"="green",
    ##                  "CC43.Micro.October2012"="green"
    ##                  )

    ## # analysis 3.P
    ## my_colors <<- c (
    ##                  "P13.WS.July2012"="orange",
    ##                  "P13.WS.October2012"="orange",
    ##                  "P24.WS.July2012"="orange",
    ##                  "P24.WS.October2012"="orange",
    ##                  "P31.WS.July2012"="orange",
    ##                  "P31.WS.October2012"="orange",
    ##                  "P46.WS.July2012"="orange",
    ##                  "P46.WS.October2012"="orange",
    ##                  "P13.LM.July2012"="red",
    ##                  "P13.LM.October2012"="red",
    ##                  "P24.LM.July2012"="red",
    ##                  "P24.LM.October2012"="red",
    ##                  "P31.LM.July2012"="red",
    ##                  "P31.LM.October2012"="red",
    ##                  "P46.LM.July2012"="red",
    ##                  "P46.LM.October2012"="red",
    ##                  "P13.MM.July2012"="blue",
    ##                  "P13.MM.October2012"="blue",
    ##                  "P24.MM.July2012"="blue",
    ##                  "P24.MM.October2012"="blue",
    ##                  "P31.MM.July2012"="blue",
    ##                  "P31.MM.October2012"="blue",
    ##                  "P46.MM.July2012"="blue",
    ##                  "P46.MM.October2012"="blue",
    ##                  "P13.SM.July2012"="magenta",
    ##                  "P13.SM.October2012"="magenta",
    ##                  "P24.SM.July2012"="magenta",
    ##                  "P24.SM.October2012"="magenta",
    ##                  "P31.SM.July2012"="magenta",
    ##                  "P31.SM.October2012"="magenta",
    ##                  "P46.SM.July2012"="magenta",
    ##                  "P46.SM.October2012"="magenta",
    ##                  "P13.Micro.July2012"="green",
    ##                  "P13.Micro.October2012"="green",
    ##                  "P24.Micro.July2012"="green",
    ##                  "P24.Micro.October2012"="green",
    ##                  "P31.Micro.July2012"="green",
    ##                  "P31.Micro.October2012"="green",
    ##                  "P46.Micro.July2012"="green",
    ##                  "P46.Micro.October2012"="green"
    ##                  )

    ## # analysis 3.PF
    ## my_colors <<- c (
    ##                  "PF15.WS.July2012"="orange",
    ##                  "PF15.WS.October2012"="orange",
    ##                  "PF23.WS.July2012"="orange",
    ##                  "PF23.WS.October2012"="orange",
    ##                  "PF32.WS.July2012"="orange",
    ##                  "PF32.WS.October2012"="orange",
    ##                  "PF41.WS.July2012"="orange",
    ##                  "PF41.WS.October2012"="orange",
    ##                  "PF15.LM.July2012"="red",
    ##                  "PF15.LM.October2012"="red",
    ##                  "PF23.LM.July2012"="red",
    ##                  "PF23.LM.October2012"="red",
    ##                  "PF32.LM.July2012"="red",
    ##                  "PF32.LM.October2012"="red",
    ##                  "PF41.LM.July2012"="red",
    ##                  "PF41.LM.October2012"="red",
    ##                  "PF15.MM.July2012"="blue",
    ##                  "PF15.MM.October2012"="blue",
    ##                  "PF23.MM.July2012"="blue",
    ##                  "PF23.MM.October2012"="blue",
    ##                  "PF32.MM.July2012"="blue",
    ##                  "PF32.MM.October2012"="blue",
    ##                  "PF41.MM.July2012"="blue",
    ##                  "PF41.MM.October2012"="blue",
    ##                  "PF15.SM.July2012"="magenta",
    ##                  "PF15.SM.October2012"="magenta",
    ##                  "PF23.SM.July2012"="magenta",
    ##                  "PF23.SM.October2012"="magenta",
    ##                  "PF32.SM.July2012"="magenta",
    ##                  "PF32.SM.October2012"="magenta",
    ##                  "PF41.SM.July2012"="magenta",
    ##                  "PF41.SM.October2012"="magenta",
    ##                  "PF15.Micro.July2012"="green",
    ##                  "PF15.Micro.October2012"="green",
    ##                  "PF23.Micro.July2012"="green",
    ##                  "PF23.Micro.October2012"="green",
    ##                  "PF32.Micro.July2012"="green",
    ##                  "PF32.Micro.October2012"="green",
    ##                  "PF41.Micro.July2012"="green",
    ##                  "PF41.Micro.October2012"="green"
    ##                  )

    ## # analysis 4
    ## my_colors <<- c (
    ##                  "CC12-LM-July2012"="#0000FFFF",
    ##                  "CC12-LM-October2012"="#0000FFFF",
    ##                  "CC12-Micro-July2012"="#00FFFFFF",
    ##                  "CC12-Micro-October2012"="#00FFFFFF",
    ##                  "CC12-MM-July2012"="#0055FFFF",
    ##                  "CC12-MM-October2012"="#0055FFFF",
    ##                  "CC12-SM-July2012"="#00AAFFFF",
    ##                  "CC12-SM-October2012"="#00AAFFFF",
    ##                  "CC12-WS-July2012"="#5500FFFF",
    ##                  "CC12-WS-October2012"="#5500FFFF",
    ##                  "CC21-LM-July2012"="#0000FFFF",
    ##                  "CC21-LM-October2012"="#0000FFFF",
    ##                  "CC21-Micro-July2012"="#00FFFFFF",
    ##                  "CC21-Micro-October2012"="#00FFFFFF",
    ##                  "CC21-MM-July2012"="#0055FFFF",
    ##                  "CC21-MM-October2012"="#0055FFFF",
    ##                  "CC21-SM-July2012"="#00AAFFFF",
    ##                  "CC21-SM-October2012"="#00AAFFFF",
    ##                  "CC21-WS-July2012"="#5500FFFF",
    ##                  "CC21-WS-October2012"="#5500FFFF",
    ##                  "CC35-LM-July2012"="#0000FFFF",
    ##                  "CC35-LM-October2012"="#0000FFFF",
    ##                  "CC35-Micro-July2012"="#00FFFFFF",
    ##                  "CC35-Micro-October2012"="#00FFFFFF",
    ##                  "CC35-MM-July2012"="#0055FFFF",
    ##                  "CC35-MM-October2012"="#0055FFFF",
    ##                  "CC35-SM-July2012"="#00AAFFFF",
    ##                  "CC35-SM-October2012"="#00AAFFFF",
    ##                  "CC35-WS-July2012"="#5500FFFF",
    ##                  "CC35-WS-October2012"="#5500FFFF",
    ##                  "CC43-LM-July2012"="#0000FFFF",
    ##                  "CC43-LM-October2012"="#0000FFFF",
    ##                  "CC43-Micro-July2012"="#00FFFFFF",
    ##                  "CC43-Micro-October2012"="#00FFFFFF",
    ##                  "CC43-MM-July2012"="#0055FFFF",
    ##                  "CC43-MM-October2012"="#0055FFFF",
    ##                  "CC43-SM-July2012"="#00AAFFFF",
    ##                  "CC43-SM-October2012"="#00AAFFFF",
    ##                  "CC43-WS-July2012"="#5500FFFF",
    ##                  "CC43-WS-October2012"="#5500FFFF",
    ##                  "P13-LM-July2012"="#FF0000FF",
    ##                  "P13-LM-October2012"="#FF0000FF",
    ##                  "P13-Micro-July2012"="#FF00FFFF",
    ##                  "P13-Micro-October2012"="#FF00FFFF",
    ##                  "P13-MM-July2012"="#FF0055FF",
    ##                  "P13-MM-October2012"="#FF0055FF",
    ##                  "P13-SM-July2012"="#FF00AAFF",
    ##                  "P13-SM-October2012"="#FF00AAFF",
    ##                  "P13-WS-July2012"="#FF5500FF",
    ##                  "P13-WS-October2012"="#FF5500FF",
    ##                  "P24-LM-July2012"="#FF0000FF",
    ##                  "P24-LM-October2012"="#FF0000FF",
    ##                  "P24-Micro-July2012"="#FF00FFFF",
    ##                  "P24-Micro-October2012"="#FF00FFFF",
    ##                  "P24-MM-July2012"="#FF0055FF",
    ##                  "P24-MM-October2012"="#FF0055FF",
    ##                  "P24-SM-July2012"="#FF00AAFF",
    ##                  "P24-SM-October2012"="#FF00AAFF",
    ##                  "P24-WS-July2012"="#FF5500FF",
    ##                  "P24-WS-October2012"="#FF5500FF",
    ##                  "P31-LM-July2012"="#FF0000FF",
    ##                  "P31-LM-October2012"="#FF0000FF",
    ##                  "P31-Micro-July2012"="#FF00FFFF",
    ##                  "P31-Micro-October2012"="#FF00FFFF",
    ##                  "P31-MM-July2012"="#FF0055FF",
    ##                  "P31-MM-October2012"="#FF0055FF",
    ##                  "P31-SM-July2012"="#FF00AAFF",
    ##                  "P31-SM-October2012"="#FF00AAFF",
    ##                  "P31-WS-July2012"="#FF5500FF",
    ##                  "P31-WS-October2012"="#FF5500FF",
    ##                  "P46-LM-July2012"="#FF0000FF",
    ##                  "P46-LM-October2012"="#FF0000FF",
    ##                  "P46-Micro-July2012"="#FF00FFFF",
    ##                  "P46-Micro-October2012"="#FF00FFFF",
    ##                  "P46-MM-July2012"="#FF0055FF",
    ##                  "P46-MM-October2012"="#FF0055FF",
    ##                  "P46-SM-July2012"="#FF00AAFF",
    ##                  "P46-SM-October2012"="#FF00AAFF",
    ##                  "P46-WS-July2012"="#FF5500FF",
    ##                  "P46-WS-October2012"="#FF5500FF",
    ##                  "PF15-LM-July2012"="#00FF00FF",
    ##                  "PF15-LM-October2012"="#00FF00FF",
    ##                  "PF15-Micro-July2012"="#FFFF00FF",
    ##                  "PF15-Micro-October2012"="#FFFF00FF",
    ##                  "PF15-MM-July2012"="#55FF00FF",
    ##                  "PF15-MM-October2012"="#55FF00FF",
    ##                  "PF15-SM-July2012"="#AAFF00FF",
    ##                  "PF15-SM-October2012"="#AAFF00FF",
    ##                  "PF15-WS-July2012"="#00FF55FF",
    ##                  "PF15-WS-October2012"="#00FF55FF",
    ##                  "PF23-LM-July2012"="#00FF00FF",
    ##                  "PF23-LM-October2012"="#00FF00FF",
    ##                  "PF23-Micro-July2012"="#FFFF00FF",
    ##                  "PF23-Micro-October2012"="#FFFF00FF",
    ##                  "PF23-MM-July2012"="#55FF00FF",
    ##                  "PF23-MM-October2012"="#55FF00FF",
    ##                  "PF23-SM-July2012"="#AAFF00FF",
    ##                  "PF23-SM-October2012"="#AAFF00FF",
    ##                  "PF23-WS-July2012"="#00FF55FF",
    ##                  "PF23-WS-October2012"="#00FF55FF",
    ##                  "PF32-LM-July2012"="#00FF00FF",
    ##                  "PF32-LM-October2012"="#00FF00FF",
    ##                  "PF32-Micro-July2012"="#FFFF00FF",
    ##                  "PF32-Micro-October2012"="#FFFF00FF",
    ##                  "PF32-MM-July2012"="#55FF00FF",
    ##                  "PF32-MM-October2012"="#55FF00FF",
    ##                  "PF32-SM-July2012"="#AAFF00FF",
    ##                  "PF32-SM-October2012"="#AAFF00FF",
    ##                  "PF32-WS-July2012"="#00FF55FF",
    ##                  "PF32-WS-October2012"="#00FF55FF",
    ##                  "PF41-LM-July2012"="#00FF00FF",
    ##                  "PF41-LM-October2012"="#00FF00FF",
    ##                  "PF41-Micro-July2012"="#FFFF00FF",
    ##                  "PF41-Micro-October2012"="#FFFF00FF",
    ##                  "PF41-MM-July2012"="#55FF00FF",
    ##                  "PF41-MM-October2012"="#55FF00FF",
    ##                  "PF41-SM-July2012"="#AAFF00FF",
    ##                  "PF41-SM-October2012"="#AAFF00FF",
    ##                  "PF41-WS-July2012"="#00FF55FF",
    ##                  "PF41-WS-October2012"="#00FF55FF"
    ##                  )


    ## ## # analysis 5
    ## my_colors <<- c (
    ##                  "CC12-LM-July2012"="#00FF00FF",
    ##                  "CC12-LM-October2012"="#00FF00FF",
    ##                  "CC12-Micro-July2012"="#00FF00FF",
    ##                  "CC12-Micro-October2012"="#00FF00FF",
    ##                  "CC12-MM-July2012"="#00FF00FF",
    ##                  "CC12-MM-October2012"="#00FF00FF",
    ##                  "CC12-SM-July2012"="#00FF00FF",
    ##                  "CC12-SM-October2012"="#00FF00FF",
    ##                  "CC12-WS-July2012"="#00FF00FF",
    ##                  "CC12-WS-October2012"="#00FF00FF",
    ##                  "CC21-LM-July2012"="#00FFFFFF",
    ##                  "CC21-LM-October2012"="#00FFFFFF",
    ##                  "CC21-Micro-July2012"="#00FFFFFF",
    ##                  "CC21-Micro-October2012"="#00FFFFFF",
    ##                  "CC21-MM-July2012"="#00FFFFFF",
    ##                  "CC21-MM-October2012"="#00FFFFFF",
    ##                  "CC21-SM-July2012"="#00FFFFFF",
    ##                  "CC21-SM-October2012"="#00FFFFFF",
    ##                  "CC21-WS-July2012"="#00FFFFFF",
    ##                  "CC21-WS-October2012"="#0080FFFF",
    ##                  "CC35-LM-July2012"="#8000FFFF",
    ##                  "CC35-LM-October2012"="#8000FFFF",
    ##                  "CC35-Micro-July2012"="#8000FFFF",
    ##                  "CC35-Micro-October2012"="#8000FFFF",
    ##                  "CC35-MM-July2012"="#8000FFFF",
    ##                  "CC35-MM-October2012"="#8000FFFF",
    ##                  "CC35-SM-July2012"="#8000FFFF",
    ##                  "CC35-SM-October2012"="#8000FFFF",
    ##                  "CC35-WS-July2012"="#8000FFFF",
    ##                  "CC35-WS-October2012"="#8000FFFF",
    ##                  "CC43-LM-July2012"="#FF0080FF",
    ##                  "CC43-LM-October2012"="#FF0080FF",
    ##                  "CC43-Micro-July2012"="#FF0080FF",
    ##                  "CC43-Micro-October2012"="#FF0080FF",
    ##                  "CC43-MM-July2012"="#FF0080FF",
    ##                  "CC43-MM-October2012"="#FF0080FF",
    ##                  "CC43-SM-July2012"="#FF0080FF",
    ##                  "CC43-SM-October2012"="#FF0080FF",
    ##                  "CC43-WS-July2012"="#FF0080FF",
    ##                  "CC43-WS-October2012"="#FF0080FF",
    ##                  "P13-LM-July2012"="#FFFF00FF",
    ##                  "P13-LM-October2012"="#FFFF00FF",
    ##                  "P13-Micro-July2012"="#FFFF00FF",
    ##                  "P13-Micro-October2012"="#FFFF00FF",
    ##                  "P13-MM-July2012"="#FFFF00FF",
    ##                  "P13-MM-October2012"="#FFFF00FF",
    ##                  "P13-SM-July2012"="#FFFF00FF",
    ##                  "P13-SM-October2012"="#FFFF00FF",
    ##                  "P13-WS-July2012"="#FFFF00FF",
    ##                  "P13-WS-October2012"="#FFFF00FF",
    ##                  "P24-LM-July2012"="#FF8000FF",
    ##                  "P24-LM-October2012"="#FF8000FF",
    ##                  "P24-Micro-July2012"="#FF8000FF",
    ##                  "P24-Micro-October2012"="#FF8000FF",
    ##                  "P24-MM-July2012"="#FF8000FF",
    ##                  "P24-MM-October2012"="#FF8000FF",
    ##                  "P24-SM-July2012"="#FF8000FF",
    ##                  "P24-SM-October2012"="#FF8000FF",
    ##                  "P24-WS-July2012"="#FF8000FF",
    ##                  "P24-WS-October2012"="#FF8000FF",
    ##                  "P31-LM-July2012"="#80FF00FF",
    ##                  "P31-LM-October2012"="#80FF00FF",
    ##                  "P31-Micro-July2012"="#80FF00FF",
    ##                  "P31-Micro-October2012"="#80FF00FF",
    ##                  "P31-MM-July2012"="#80FF00FF",
    ##                  "P31-MM-October2012"="#80FF00FF",
    ##                  "P31-SM-July2012"="#80FF00FF",
    ##                  "P31-SM-October2012"="#80FF00FF",
    ##                  "P31-WS-July2012"="#80FF00FF",
    ##                  "P31-WS-October2012"="#80FF00FF",
    ##                  "P46-LM-July2012"="#FF0000FF",
    ##                  "P46-LM-October2012"="#FF0000FF",
    ##                  "P46-Micro-July2012"="#FF0000FF",
    ##                  "P46-Micro-October2012"="#FF0000FF",
    ##                  "P46-MM-July2012"="#FF0000FF",
    ##                  "P46-MM-October2012"="#FF0000FF",
    ##                  "P46-SM-July2012"="#FF0000FF",
    ##                  "P46-SM-October2012"="#FF0000FF",
    ##                  "P46-WS-July2012"="#FF0000FF",
    ##                  "P46-WS-October2012"="#FF0000FF",
    ##                  "PF15-LM-July2012"="#00FF80FF",
    ##                  "PF15-LM-October2012"="#00FF80FF",
    ##                  "PF15-Micro-July2012"="#00FF80FF",
    ##                  "PF15-Micro-October2012"="#00FF80FF",
    ##                  "PF15-MM-July2012"="#00FF80FF",
    ##                  "PF15-MM-October2012"="#00FF80FF",
    ##                  "PF15-SM-July2012"="#00FF80FF",
    ##                  "PF15-SM-October2012"="#00FF80FF",
    ##                  "PF15-WS-July2012"="#00FF80FF",
    ##                  "PF15-WS-October2012"="#00FF80FF",
    ##                  "PF23-LM-July2012"="#0080FFFF",
    ##                  "PF23-LM-October2012"="#0080FFFF",
    ##                  "PF23-Micro-July2012"="#0080FFFF",
    ##                  "PF23-Micro-October2012"="#0080FFFF",
    ##                  "PF23-MM-July2012"="#0080FFFF",
    ##                  "PF23-MM-October2012"="#0080FFFF",
    ##                  "PF23-SM-July2012"="#0080FFFF",
    ##                  "PF23-SM-October2012"="#0080FFFF",
    ##                  "PF23-WS-July2012"="#0080FFFF",
    ##                  "PF23-WS-October2012"="#0080FFFF",
    ##                  "PF32-LM-July2012"="#0000FFFF",
    ##                  "PF32-LM-October2012"="#0000FFFF",
    ##                  "PF32-Micro-July2012"="#0000FFFF",
    ##                  "PF32-Micro-October2012"="#0000FFFF",
    ##                  "PF32-MM-July2012"="#0000FFFF",
    ##                  "PF32-MM-October2012"="#0000FFFF",
    ##                  "PF32-SM-July2012"="#0000FFFF",
    ##                  "PF32-SM-October2012"="#0000FFFF",
    ##                  "PF32-WS-July2012"="#0000FFFF",
    ##                  "PF32-WS-October2012"="#0000FFFF",
    ##                  "PF41-LM-July2012"="#FF00FFFF",
    ##                  "PF41-LM-October2012"="#FF00FFFF",
    ##                  "PF41-Micro-July2012"="#FF00FFFF",
    ##                  "PF41-Micro-October2012"="#FF00FFFF",
    ##                  "PF41-MM-July2012"="#FF00FFFF",
    ##                  "PF41-MM-October2012"="#FF00FFFF",
    ##                  "PF41-SM-July2012"="#FF00FFFF",
    ##                  "PF41-SM-October2012"="#FF00FFFF",
    ##                  "PF41-WS-July2012"="#FF00FFFF",
    ##                  "PF41-WS-October2012"="#FF00FFFF"
    ##                  )
    
    
    points(x=((sorted_my_data[,PC1])), y=((sorted_my_data[,PC2])), pch=18, col=my_colors, bg =my_colors, cex=my_cex)
    
    #points(x=((sorted_my_data[,PC1])[1:60]), y=((sorted_my_data[,PC2])[1:60]), pch=23, col="green", bg = "green", cex=my_cex) #C
    #points(x=((sorted_my_data[,PC1])[61:120]), y=((sorted_my_data[,PC2])[61:120]), pch=21, col="blue", bg = "blue", cex=my_cex) #ch1
    #points(x=((sorted_my_data[,PC1])[11:15]), y=((sorted_my_data[,PC2])[11:15]), pch=22, col="purple", bg = "purple", cex=my_cex) #ch2
    #points(x=((sorted_my_data[,PC1])[16:20]), y=((sorted_my_data[,PC2])[16:20]), pch=23, col="purple", bg = "purple", cex=my_cex) #ch3
    #points(x=((sorted_my_data[,PC1])[21:25]), y=((sorted_my_data[,PC2])[21:25]), pch=24, col="purple", bg = "purple", cex=my_cex) #ch4
    #points(x=((sorted_my_data[,PC1])[26:30]), y=((sorted_my_data[,PC2])[26:30]), pch=25, col="purple", bg = "purple", cex=my_cex) #ch5
    #points(x=((sorted_my_data[,PC1])[31:35]), y=((sorted_my_data[,PC2])[31:35]), pch=21, col="blue", bg = "blue", cex=my_cex) #E1
    #points(x=((sorted_my_data[,PC1])[36:40]), y=((sorted_my_data[,PC2])[36:40]), pch="X", col="orange", bg = "orange", cex=0.7*my_cex, lwd=3) #R
    #points(x=((sorted_my_data[,PC1])[41:45]), y=((sorted_my_data[,PC2])[41:45]), pch="+", col="brown", bg = "brown", cex=my_cex, lwd=3) #W
 


    #title( (paste(file_in,"\n", "PC", PC1, "vs PC", PC2 )), cex.main = 0.8)
    title( (paste(out_prefix,"\n", "PC", PC1, "vs PC", PC2 )), cex.main = 0.8)
    dev.off()
    
    pdf (file=legend_out) #, width = 2, height = 4)

    ## # analysis 1
    ## legend_colors <- c ('red','green','orange','blue','magenta','purple')
    ## names (legend_colors) <- c ('Jul.CC','Jul.P','Jul.PF','Oct.CC','Oct.P','Oct.PF')

    ## # analysis 2
    ## group_colors <- c ('red','blue','green')
    ## names (group_colors) <- c ('CC','P','PF')

    ## # analysis 3
    ## group_colors <- c ('orange','red','blue','magenta','green')
    ## names (group_colors) <- c ('WS','LM','MM','SM','Micro')

    ## # analysis 5
    ## group_colors <- c ('#FF0000FF','#FF8000FF','#FFFF00FF','#80FF00FF','#00FF00FF','#00FF80FF','#00FFFFFF','#0080FFFF','#0000FFFF','#8000FFFF','#FF00FFFF','#FF0080FF')
    ## names (group_colors) <- c ('negsix','negfour','negthree','negone','twelve','fifteen','twentyone','twentythree','thirtytwo','thirtyfive','fourtyone','fourtythree')
    
    plot.new ()
    legend (0, 1, legend = names (legend_colors), pch = 18, col = group_colors )
    dev.off ()




  }

