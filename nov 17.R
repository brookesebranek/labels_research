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
errors_frame<-PeakError::PeakError(peaks_table[status=="peak"], track_labels_peaks_table)
errors_table<-setDT(errors_frame)
errors_table[status!="correct"]
data.table::setkey(peaks_table, chrom, chromStart, chromEnd)
data.table::setkey(track_labels_table, chrom, chromStart, chromEnd)
data.table::setkey(coverage_table, chrom, chromStart, chromEnd)
data.table::setkey(errors_table, chrom, chromStart, chromEnd)
for(to do title){
window_table<- data.table::data.table(chrom="chr10", windowStart=chromStart-bases_error, windowEnd=chromEnd+bases_around_error))
data.table::setkey(window_table, chrom, windowStart, windowEnd)
peaks_and_window_table<-data.table::foverlaps(peaks_table, window_table, nomatch=NULL)
track_labels_and_window_table<-data.table::foverlaps(track_labels_table, window_table, nomatch = NULL)
coverage_and_window_table<-data.table::foverlaps(coverage_table, window_table, nomatch=NULL)
errors_and_window_table<-data.table::foverlaps(errors_table, window_table, nomatch=NULL)
incorrect_errors_table<-errors_and_window_table[status!="correct"]
library(ggplot2) #must run every time
bases_around_error<-1000
for(row_number in 1:nrow(errors_and_window_table)){
  gg<-ggplot2::ggplot()+geom_path(
    aes(x=chromStart, y=matches), color="gray50", data=coverage_and_window_table)+
    geom_segment(
      aes(x=chromStart, y=mean, xend=chromEnd, yend=mean), color="orange", size=3, data= peaks_and_window_table)+
    geom_rect(aes(
      xmin=chromStart, ymin=-Inf, xmax=chromEnd, ymax=Inf, fill=annotation), alpha= 0.25, data=track_labels_and_window_table)+
    geom_rect(aes(
      xmin=chromStart, ymin=-Inf, xmax=chromEnd, ymax=Inf,  linetype=status), fill="transparent", color="black", size=1,
      data=incorrect_errors_table)+
    scale_linetype_manual(
      "error type",
      values=c(
        correct=0,
        "false negative"=3,
        "false positive"=1))+
    errors_and_window_table[row_number ,coord_cartesian(xlim=c(chromStart-bases_around_error, chromEnd+bases_around_error))]
  png(paste0("nov11_errors_", row_number, ".png"), width=10, height = 4, units = "in", res = 200)   
  print(gg)
  dev.off()
}
