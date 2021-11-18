peaks_filename<-"/Users/brookesebranek/H3K4me3_TDH_ENCODE/samples/skeletalMuscle/ENCFF280HQO/problems/chr10:60000-17974675/coverage.bedgraph_penalty=4303.13667437155_segments.bed"
contig_directory<-dirname(peaks_filename)
problem_directory<-dirname(contig_directory)
encode_directory<-dirname(problem_directory)
sample_type_directory<-dirname(encode_directory)
sample_type_basename<- basename(sample_type_directory)
encode_basename<-basename(encode_directory)
sample_encode<-paste0(sample_type_basename,"_", encode_basename)
labels_filename<-"https://raw.githubusercontent.com/tdhock/peaklearner-api/main/labels/H3K4me3_TDH_ENCODE/2021-09-21_14_55_10/labels.csv"
coverage_filename<-sub("_penalty.*", "", peaks_filename)
library(data.table) #must install every session
peaks_table<-fread(peaks_filename, col.names = c("chrom", "chromStart", "chromEnd", "status", "mean"))
labels_table<-fread(labels_filename, col.names = c("annotation","chrom", "chromEnd", "chromStart", "createdBy", "modifiedBy", "lastModified", "track", "dupe"))
track_labels_table<-labels_table[sample_encode==track]
track_labels_peaks_table<-track_labels_table[, .(chrom, chromStart, chromEnd, track, annotation=sub("noPeak", "noPeaks", annotation))]
coverage_table<-fread(coverage_filename, col.names = c("chrom", "chromStart", "chromEnd", "matches"))
window_table<- data.table::data.table(chrom="chr10", windowStart=15200000, windowEnd=15300000)
errors_frame<-PeakError::PeakError(peaks_table[status=="peak"], track_labels_peaks_table)
errors_table<-setDT(errors_frame)
data.table::setkey(peaks_table, chrom, chromStart, chromEnd)
data.table::setkey(track_labels_table, chrom, chromStart, chromEnd)
data.table::setkey(coverage_table, chrom, chromStart, chromEnd)
data.table::setkey(window_table, chrom, windowStart, windowEnd)
data.table::setkey(errors_table, chrom, chromStart, chromEnd)
