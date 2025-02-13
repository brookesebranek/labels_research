peaks_filename<-"/Users/brookesebranek/H3K4me3_TDH_ENCODE/samples/skeletalMuscle/ENCFF280HQO/problems/chr10:60000-17974675/coverage.bedgraph_penalty=4303.13667437155_segments.bed"
contigStart<-60000
contigEnd<- 17974675
contig_table<- data.table::data.table(chrom="chr10", contigStart, contigEnd)
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
labels_table<-fread(labels_filename, col.names = c(
  "annotation","chrom", "chromEnd", "chromStart", "createdBy", "modifiedBy", "lastModified", "track", "dupe"))
track_labels_table<-labels_table[sample_encode==track]
track_labels_peaks_table<-track_labels_table[, .(chrom, chromStart, chromEnd, track, annotation=sub("noPeak", "noPeaks", annotation))]
coverage_table<-fread(coverage_filename, col.names = c("chrom", "chromStart", "chromEnd", "matches"))
data.table::setkey(peaks_table, chrom, chromStart, chromEnd)
data.table::setkey(track_labels_peaks_table, chrom, chromStart, chromEnd)
data.table::setkey(coverage_table, chrom, chromStart, chromEnd)
data.table::setkey(contig_table, chrom, contigStart, contigEnd)
labels_and_contig_table<-foverlaps(contig_table, track_labels_peaks_table)
errors_frame<-PeakError::PeakError(peaks_table[status=="peak"], labels_and_contig_table)
errors_and_correct_frame<-setDT(errors_frame)
errors_table<-errors_and_correct_frame[status!="correct"]
data.table::setkey(errors_table, chrom, chromStart, chromEnd)
#ALL TABLES BEFORE HERE OK
#
bases_around_error<-100000
for(row_of_window in 1:nrow(errors_table)){
  current_error<-errors_table[row_of_window]
  window_table<-current_error[,data.table::data.table(
    chrom, windowStart=chromStart-bases_around_error, windowEnd=chromEnd+bases_around_error)]
  data.table::setkey(window_table, chrom, windowStart, windowEnd)
  peaks_and_window_table<-data.table::foverlaps(peaks_table, window_table, nomatch=NULL)
  track_labels_and_window_table<-data.table::foverlaps(labels_and_contig_table, window_table)
  coverage_and_window_table<-data.table::foverlaps(coverage_table, window_table, nomatch=NULL)
  errors_and_window_table<-data.table::foverlaps(errors_table, window_table, nomatch=NULL)
  library(ggplot2) #must run every time
    gg<-ggplot2::ggplot()+geom_path(aes(x=chromStart, y=matches), color="gray50", data=coverage_and_window_table)+
      geom_segment(
        aes(
          x=chromStart, y=mean, xend=chromEnd, yend=mean), color="orange", 
        size=3, data= peaks_and_window_table)+
      geom_rect(aes(
        xmin=chromStart, ymin=-Inf, xmax=chromEnd, ymax=Inf, fill=annotation), 
        alpha= 0.25, data=track_labels_and_window_table)+
      geom_rect(aes(
        xmin=chromStart, ymin=-Inf, xmax=chromEnd, ymax=Inf,  linetype=status), 
        fill="transparent", color="black", size=1,
        data=errors_and_window_table)+
      scale_linetype_manual(
        "error type",
        values=c(
          correct=0,
          "false negative"=3,
          "false positive"=1))+
      window_table[,coord_cartesian(xlim=c(windowStart, windowEnd))]
    png(paste0("jan12_errors2.0_", row_of_window, ".png"), width=10, height = 4, units = "in", res = 200)   
    print(gg)
    dev.off()
}

# wikipedia regular expressions for more info.
int.pattern <- list("[0-9]*", as.integer)
start.end.pattern <- list(
  chrom="chr.*", ":", start=int.pattern, "-", end=int.pattern)
nc::capture_first_vec(basename(contig_directory), start.end.pattern)
setDT(x=start.end.pattern, )

