
args <- commandArgs(trailingOnly = TRUE) #import arguments

if (length(args)<3) {
  stop("At least 3 arguments must be supplied (patient's name, input file, output file).n", call.=FALSE)
} else {
	patient=args[1] #$id 
	file=args[2] #$output.tsv 
	output=args[3] #$output.pdf 
	size=args[4] #$size of ROH
}

# chromosome names and positoins
chr_name <- paste("chr", seq(1, 22, 1), sep = "") #list de 1 à 22
#position of the middle of the chrm
chr_position <- c(124478211,370053186,590297730,784552788,970429194,1146601313,1311677289,1463919594,1605686270,1741782340,1876224362,2010405327,2134225146,2244929169,2349446622,2445611389,2532409282,2614224645,2683720096,2745250987,2800828062,2849592288)# nolint
#position of the chrm limits
chr_limits <- c(0,248956422, 491149951, 689445510, 879660065, 1061198324, 1232004303, 1391350276, 1536488912, 1674883629, 1808681051, 1943767673, 2077042982, 2191407310, 2298451028, 2400442217, 2490780562, 2574038003, 2654411288, 2713028904, 2777473071, 2824183054, 2875001522) # nolint

######chrm lenght(1 to 22)
chr_len <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468) # nolint

data <- read.table(file, sep = "\t", # import the homozygous region file
comment.char = "@", #delete lines with a #
fill = TRUE) #fill empty box

if (dim(data)[1] > 4) { #check if there's enought ROH
    data <- read.table(file, header = FALSE, #importe without header
    sep = "\t") 

    data <- subset(data, data$V1 != "chrX") # remove the chrX positions

    ## PRINT THE HOMOZYGOSITY DATA

    chr_pos <- as.numeric(substr(data$V1, 4, 5)) # extract the chromosome number
    #Vector of the chrcorresponding to each ROH (1 1 1 3 3 5 6 9 12..)

    # extract the left position (position on the chromosome + chromosome limit)
    xleft <- data$V2 + chr_limits[chr_pos]

    # extract the right position (position on the chromosome + chromosome limit)
    xright <- data$V3 + chr_limits[chr_pos]

##########################################################################
    floor_roh <- 0.20 #detect UPD when xx% of a chrm is composed of ROH###
#########################################################################
    moins <- xright - xleft #ROH lenght

    chr_len_pos <- chr_len[chr_pos]

    compil <- data.frame(chr_pos,moins)

    compil_tri <- tapply(compil$moins, compil$chr_pos, sum)
    compil_tri <- t(compil_tri)
    compil_tri <- data.matrix(compil_tri)
    compil_tri <- t(compil_tri) #ROH lenght on each chrm

    chr_len_tri <- data.frame(chr_pos, chr_len_pos)
    chr_len_tri <- data.matrix(t(unique(chr_len_tri)))
    chr_len_tri <- t(chr_len_tri) #lenght of the chrm that contain ROH

###ROH proportion on each chrm
    div_matrix <- compil_tri[,1] / chr_len_tri[,2]
    roh_max <- max(div_matrix) 

    if (roh_max > floor_roh) { 

    roh_max <- floor(roh_max * 100)
    output <- paste(output, "_Possible_disomie_", roh_max, "percent.pdf", sep = "")
    }
}


pdf(output, width = 10, height = 3)

par(mar = c(5.1, 2.1, 4.1, 2.1))

plot(
c(0, max(chr_limits)), c(0, 100), 
xlim = c(1, max(chr_limits)), 
yaxs = "i", xaxs = "i", type = "n", 
xlab = "", ylab = "", 
axes = FALSE, 
main = paste("Homozygous Regions for ", patient,
"\n Total = ", size, " Mb (autosomes)", " ROH_max : ",roh_max, " %", sep = ""))

box() 

## For the x axis
axis(1, tck = 0, at = chr_position,
labels = chr_name, las = 2, cex = 0.8) #chr name
axis(side = 1, tck = -.10, labels = NA, at = chr_limits) 

## For the y axis
# axis(side=2,tck=-.010, at=c(0, 25, 50, 75, 100), labels=c(0, 25, 50, 75, 100),las=1)

## PLOT THE CHROMOSOME LIMIT"S
abline(v = chr_limits[2:22], col = "grey", lty = 2) 

if(dim(data)[1]>4){
 
   if(roh_max > floor_roh ){
   rect(0, 80, 2875001522, 100, col = "red", border = NA) #red rectangle
   }

  rect(xleft, 0, xright, 100, col = "blue", border = NA) #blue rectangles
}

dev.off()
