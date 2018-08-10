#!/usr/bin/env Rscript

## CPCantalapiedra 2017-2018
##
## marker_LD_region
##
## An R script to define an interval around a given marker
## based on its LD decay with surrounding markers.
##

# Load R libraries
suppressMessages(library(methods))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

write("Running marker_LD_region.R...", file=stderr())

## Read parameters
args = commandArgs(trailingOnly=TRUE)

#marker_id <- args[1] # a list of markers to analyze
#marker_chrom <- args[2]
#marker_start <- as.numeric(args[3])
#marker_end <- as.numeric(args[4])
#max_interval <- as.numeric(args[5])
#MIN_INTERVAL <- as.numeric(args[6])
#hapmap_filename <- args[7]
#WINDOW_SIZE <- as.numeric(args[8]) # 10
#WINDOW_STEP <- as.numeric(args[9]) # 5
#LDTHRES <- as.numeric(args[10]) # 0.2

## mandatory parameters
marker_id <- args[1] # a list of markers to analyze
hapmap_filename <- args[2] # the genotyping file to be zoom in
marker_chrom <- args[2] # chromosome
marker_start <- as.numeric(args[3]) # start position
marker_end <- as.numeric(args[4]) # end position

## recommended parameters
LDTHRES <- as.numeric(args[5]) # LD threshold which will define the ending interval

## other parameters
max_interval <- as.numeric(args[6]) # max ending interval
MIN_INTERVAL <- as.numeric(args[7]) # min ending interval
WINDOW_SIZE <- as.numeric(args[8]) # number of markers to compute LD. def: 10
WINDOW_STEP <- as.numeric(args[9]) # number of markers to skip between each window. def: 5

cat("Parameters:\n", file=stderr())
cat(paste("\t", marker_id, "\n"), file=stderr())
cat(paste("\t", hapmap_filename, "\n"), file=stderr())
cat(paste("\t", marker_chrom, "\n"), file=stderr())
cat(paste("\t", marker_start, "\n"), file=stderr())
cat(paste("\t", marker_end, "\n"), file=stderr())
cat(paste("\t", LDTHRES, "\n"), file=stderr())
cat(paste("\t", max_interval, "\n"), file=stderr())
cat(paste("\t", MIN_INTERVAL, "\n"), file=stderr())
cat(paste("\t", WINDOW_SIZE, "\n"), file=stderr())
cat(paste("\t", WINDOW_STEP, "\n"), file=stderr())

########################################### Functions
###########################################

### Functions used to format markers as needed to compute LD

f_replace_alleles <- function(x, a, b, h1, h2){
    x <- gsub(a, 0, x);
    x <- gsub(b, 2, x);
    x <- gsub(h1, 1, x);
    x <- gsub(h2, 1, x);
    
    return(x);
}

f_ldcorsv_format_row <- function(x, ldcorsv_df){
    
    # retrieve marker data
    ldcorsv_df_row <- ldcorsv_df[x,]
    
    #cat(paste("Formatting marker ", unlist(as.character(ldcorsv_df_row$rs)), "\n", sep=""), file=stderr())
    
    # substitute "NN" to NA
    ldcorsv_df_row <- data.frame(lapply(ldcorsv_df_row, function(x) {gsub("NN", NA, x);}))
    
    datarow <- unname(unlist(ldcorsv_df_row))
    
    datarow <- datarow[! is.na(datarow)]
    
    # obtain the different genotypes, only once, to obtain the alleles
    datarow <- unique(datarow)
    
    alleles = as.character(datarow[[2]]);
    
    allelessplit <- unlist(strsplit(alleles, "/"))
    allele_a <- paste(allelessplit[[1]], allelessplit[[1]], sep="")
    allele_b <- paste(allelessplit[[2]], allelessplit[[2]], sep="")
    allele_ab <- paste(allelessplit[[1]], allelessplit[[2]], sep="")
    allele_ba <- paste(allelessplit[[2]], allelessplit[[1]], sep="")
    
    # reformat the alleles to new format
    ldcorsv_df_row <- data.frame(lapply(ldcorsv_df_row,
                        f_replace_alleles, allele_a, allele_b, allele_ab, allele_ba))
    
    return(ldcorsv_df_row);
}

f_ldcorsv_format <- function(hapmap_df){
    
    ldcorsv_df <- hapmap_df %>% select(-c, -pos, -strand, -assembly, -center,
                                    -protLSID, -assayLSID, -panelLSID, -QCCode)
    
    ldcorsv_df <- unique(ldcorsv_df)
    
    ldcorsv_df <- do.call("rbind", lapply(seq(1:nrow(ldcorsv_df)), f_ldcorsv_format_row, ldcorsv_df))
    
    #cat("Creating formatted data frame...\n", file=stderr())
    
    rownames(ldcorsv_df) <- ldcorsv_df$rs
    ldcorsv_df <- ldcorsv_df %>% select(-alleles, -rs)
    
    # check.names FALSE is very important to avoid
    # changes in marker names by R automatically (e.g. BOPA1_1283-332 --> BOPA1_1283.332)
    ldcorsv_df <- data.frame(t(ldcorsv_df), check.names = FALSE)
    
    return(ldcorsv_df)
}

### Functions used to compute LD

f_compute_marker_LD <- function(x, ldcorsv_df, marker_id){
    #cat(paste("Computing LD of marker ", x, "\n", sep=""), file=stderr())
    markerdf <- ldcorsv_df %>% select(x, marker_id)
    
    # call to LDcorSV function LD.Measures
    markerLD = LD.Measures(markerdf,V=NA,S=NA,data ="G",supinfo=TRUE,na.presence=TRUE)
    
    return(markerLD)
}

f_compute_LD <- function(ldcorsv_df, marker_id){
    markers_LD <- colnames(ldcorsv_df)
    markers_LD <- markers_LD[markers_LD != marker_id]
    LDdf <- do.call("rbind", lapply(markers_LD, f_compute_marker_LD, ldcorsv_df, marker_id))
    
    return(LDdf)
}

### Functions used to compute LD decay in windows

f_compute_window_LD <- function(x, hapmap_df, LD, WINDOW_SIZE){
    start = ifelse(x - WINDOW_SIZE > 0, x - WINDOW_SIZE, 1)
    end = ifelse(x + WINDOW_SIZE <= nrow(hapmap_df), x + WINDOW_SIZE, nrow(hapmap_df))
    markers_df <- hapmap_df[start:end,]
    
    markers_r2 <- LD[LD$loc1 %in% markers_df$rs,]$r2
    markers_r2_90percen <- quantile(markers_r2, probs = c(0.9))
    
    markers_pos <- markers_df$pos
    
    avg_pos = round(mean(markers_pos))
    
    window_df <- data.frame(pos = avg_pos, r2 = markers_r2_90percen, row.names = NULL)
    
    return(window_df);
}

f_compute_windows_LD <- function(hapmap_df, LD, WINDOW_SIZE, WINDOW_STEP){
    LDwindows <- do.call("rbind", lapply(seq(WINDOW_SIZE+1, nrow(hapmap_df)-WINDOW_SIZE, WINDOW_STEP),
                        f_compute_window_LD, hapmap_df, LD, WINDOW_SIZE));
    
    return(LDwindows);
}

########################################### BEGIN
###########################################

gmat <- NA
# This is done here because add.expr of heatmap.2 only searches variables
# in the global environment, and not in the local function from where heatmap.2 is called
# This was known as a bug already in 2015

# Reading markers list
cat("Reading hapmap file... this could take a while...\n", file=stderr())
hapmap_df <- read.table(hapmap_filename, header = TRUE)

# Remove markers in other chromosomes
hapmap_df <- hapmap_df[hapmap_df$c==marker_chrom,] %>% arrange(c, pos)

#print(head(hapmap_df))

## Transform the data to LDcorSV format
## columns: markers
## rows: genotypes
## values: 0,1,2
##

cat("Formatting markers... please be patient...\n", file=stderr())
ldcorsv_df <- f_ldcorsv_format(hapmap_df)

#print(head(ldcorsv_df))

# I need to change the type of data of the data.frame
# from "factor" (as can be seen running sapply(lodcorsv_df, class))
# to "numeric". I do this with the library "varhandle":
library(varhandle)
ldcorsv_df <- unfactor(ldcorsv_df)

## Compute LD for all the markers
cat("Computing LD... just a bit more of patience... ;)\n", file=stderr())
library(LDcorSV)

LD <- f_compute_LD(ldcorsv_df, marker_id)

# LD fields are: loc1, loc2, r2, ... others not so important for this script

# Compute LD for windows

cat("Computing LD decay... we are almost done... ;)\n", file=stderr())
LDwindows <- f_compute_windows_LD(hapmap_df, LD, WINDOW_SIZE, WINDOW_STEP)

#print(LDwindows)

# markers to plot along with the LD windows

#hapmap_pos = hapmap_df %>% select(pos) %>% mutate(r2 = 0)

LDoutfile=paste("LD/", marker_id, ".", marker_chrom, ".", marker_start, ".r2.png", sep="")
png(LDoutfile, width=400, height=400)
plot(LDwindows, pch = 19, cex = 0.5, ylim = c(0,1))
#points(hapmap_pos, pch = 19, cex = 0.5, col = "red")
devnull <- dev.off()

# write the tables also, so that the plot could be done later without running
# everything again

LDtsvfile=paste("LD/", marker_id, ".", marker_chrom, ".", marker_start, ".r2.tsv", sep="")
write.table(LDwindows, file=LDtsvfile, sep="\t")

#############################
cat("Defining LD threshold region...\n", file=stderr())

## Look for LD windows threshold

# upstream

upstream_windows <- LDwindows[LDwindows$pos < marker_start,] %>% arrange(desc(pos))

prev_pos = -1
for (i in seq(1, nrow(upstream_windows))){
    curr_window = upstream_windows[i,]
    if (curr_window$r2 >= LDTHRES) {
        prev_pos = curr_window$pos;
    } else {
        break;
    }
}

MIN_INTERVAL = 100000
if (prev_pos == -1){
    upstream_pos = ifelse(marker_start - MIN_INTERVAL > 0, marker_start - MIN_INTERVAL, 1);
    
} else if (marker_start - prev_pos < MIN_INTERVAL) {
    upstream_pos = ifelse(marker_start - MIN_INTERVAL > 0, marker_start - MIN_INTERVAL, 1);
    
} else {
    upstream_pos = prev_pos;
}

# downstream

downstream_windows <- LDwindows[LDwindows$pos > marker_start,] %>% arrange(pos)

prev_pos = -1
for (i in seq(1, nrow(downstream_windows))){
    curr_window = downstream_windows[i,]
    if (curr_window$r2 >= LDTHRES) {
        prev_pos = curr_window$pos;
    } else {
        break;
    }
}

if (prev_pos == -1){
    downstream_pos = marker_start + MIN_INTERVAL
} else if (prev_pos - marker_start < MIN_INTERVAL) {
    downstream_pos = marker_start + MIN_INTERVAL;
} else {
    downstream_pos = prev_pos
}

cat(paste("The final interval defined by LD will span: ", upstream_pos, " - ", downstream_pos, "\n", sep=""), file=stderr())

LDregion <- LDwindows[LDwindows$pos >= upstream_pos &
                        LDwindows$pos <= downstream_pos,]

#hapmap_pos = hapmap_pos[hapmap_pos$pos >= upstream_pos &
#                        hapmap_pos$pos <= downstream_pos,]

LDoutfile=paste("LD/", marker_id, ".", marker_chrom, ".", marker_start, ".r2.LDregion.png", sep="")
png(LDoutfile, width=400, height=400)
plot(LDregion, pch = 19, cex = 0.5, ylim = c(0, 1))
#points(hapmap_pos, pch = 19, cex = 0.5, col = "red")
devnull <- dev.off()

LDtsvfile=paste("LD/", marker_id, ".", marker_chrom, ".", marker_start, ".r2.LDregion.tsv", sep="")
write.table(LDregion, file=LDtsvfile, sep="\t")

cat("Finished!\n", file=stderr())

q()

## END
