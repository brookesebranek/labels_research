peaks_filename<-"/Users/brookesebranek/H3K4me3_TDH_ENCODE/samples/skeletalMuscle/ENCFF280HQO/problems/chr10:60000-17974675/coverage.bedgraph_penalty=4303.13667437155_segments.bed"
labels_filename<-"https://raw.githubusercontent.com/tdhock/peaklearner-api/main/labels/H3K4me3_TDH_ENCODE/2021-09-21_14_55_10/labels.csv"
coverage_filename<-sub("_penalty.*", "", peaks_filename)
library(data.table) #must install every session
peaks_table<-fread(peaks_filename, col.names = c("chrom", "chromStart", "chromEnd", "status", "penalty"))
labels_table<-fread(labels_filename, col.names = c("labels","chrom", "chromEnd", "chromStart", "createdBy", "modifiedBy", "lastModified", "track", "dupe"))
coverage_table<-fread(coverage_filename, col.names = c("chrom", "chromStart", "chromEnd", "matches"))
window_table<- data.table::data.table(chrom="chr10", windowStart=15200000, windowEnd=15300000)
data.table::setkey(peaks_table, chrom, chromStart, chromEnd)
data.table::setkey(labels_table, chrom, chromStart, chromEnd)
data.table::setkey(coverage_table, chrom, chromStart, chromEnd)
data.table::setkey(window_table, chrom, windowStart, windowEnd)
peaks_and_window_table<-data.table::foverlaps(peaks_table, window_table, nomatch=NULL)
labels_and_window_table<-data.table::foverlaps(labels_table, window_table, nomatch=NULL)
coverage_and_window_table<-data.table::foverlaps(coverage_table, window_table, nomatch=NULL)
library(ggplot2) #must run every time
ggplot2::ggplot()+geom_path(aes(x=chromStart, y=matches), data=coverage_and_window_table)
ggplot2::ggplot()+geom_segment(aes(x=chromStart, y=status, xend=chromEnd, yend=status), data= peaks_and_window_table)


