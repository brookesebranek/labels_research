models_rds_vec<-Sys.glob("~/H3K4me3_TDH_ENCODE/samples/*/*/problems/chr10:60000-17974675/models.rds")
labels_filename<-"https://raw.githubusercontent.com/tdhock/peaklearner-api/main/labels/H3K4me3_TDH_ENCODE/2021-09-21_14_55_10/labels.csv"
for(file_number_of_rds in 1:length(models_rds_vec)){
  peaks_filename<-(models_rds_vec[[file_number_of_rds]])
  rds_table<-readRDS(peaks_filename)
  if(nrow(rds_table)> 0){
    contig_table<- nc::capture_first_vec(
      peaks_filename, chrom="chr.*", ":", 
      contigStart="[0-9]+", as.integer,
      "-",
      contigEnd="[0-9]+", as.integer)
    contig_directory<-dirname(peaks_filename)
    problem_directory<-dirname(contig_directory)
    encode_directory<-dirname(problem_directory)
    sample_type_directory<-dirname(encode_directory)
    sample_type_basename<- basename(sample_type_directory)
    encode_basename<-basename(encode_directory)
    sample_encode<-paste0(sample_type_basename,"_", encode_basename)
    coverage_filename<-sub("models.rds", "coverage.bedgraph", peaks_filename)
    library(data.table)
    best_model_info<-rds_table[which.min(fp+fn),]
    peaks_table<-best_model_info[, segments.dt[[1]] ]
    labels_table<-fread(labels_filename, col.names = c(
      "annotation","chrom", "chromEnd", "chromStart", "createdBy", "modifiedBy", "lastModified", "track", "dupe"))
    track_labels_table<-labels_table[sample_encode==track]
    track_labels_peaks_table<-track_labels_table[, .(chrom, chromStart, chromEnd, track, annotation=sub("noPeak", "noPeaks", annotation))]
    coverage_table<-fread(coverage_filename, col.names = c("chrom", "chromStart", "chromEnd", "matches"))
    data.table::setkey(peaks_table, chromStart, chromEnd)
    data.table::setkey(track_labels_peaks_table, chromStart, chromEnd)
    data.table::setkey(coverage_table,chromStart, chromEnd)
    data.table::setkey(contig_table, contigStart, contigEnd)
    labels_and_contig_table<-foverlaps(contig_table, track_labels_peaks_table)
    errors_before_correct<-best_model_info[, errors.dt[[1]] ]
    errors_table<-errors_before_correct[status!="correct"]
    data.table::setkey(errors_table, chromStart, chromEnd)
  todays_date.<-Sys.Date()
  setwd("/Users/brookesebranek/Desktop/labels_research./pngs_for_5_samples")
  bases_around_error<-10000
  for(row_of_window in 1:nrow(errors_table)){
    current_error<-errors_table[row_of_window]
    window_table<-current_error[,data.table::data.table(
      windowStart=chromStart-bases_around_error, windowEnd=chromEnd+bases_around_error)]
    data.table::setkey(window_table, windowStart, windowEnd)
    peaks_and_window_table<-data.table::foverlaps(peaks_table, window_table, nomatch = NULL)
    track_labels_and_window_table<-data.table::foverlaps(labels_and_contig_table, window_table)
    coverage_and_window_table<-data.table::foverlaps(coverage_table, window_table, nomatch = NULL)
    errors_and_window_table<-data.table::foverlaps(errors_table, window_table, nomatch=NULL)
    library(ggplot2) 
    gg<-ggplot2::ggplot()+geom_path(aes(x=chromStart, y=matches), color="gray50", data=coverage_and_window_table)+
      geom_segment(
        aes(
          x=chromStart, y=mean, xend=chromEnd, yend=mean), color="orange", 
        size=3, data= peaks_and_window_table)+
      geom_rect(aes(
        xmin=chromStart, ymin=-Inf, xmax=chromEnd, ymax=Inf, fill=annotation), 
        alpha= 0.10, data=track_labels_and_window_table)+
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
    png(paste0(todays_date., sample_encode, row_of_window, ".png"), width=10, height = 4, units = "in", res = 200)   
    print(gg)
    dev.off()
  }
  }
}

