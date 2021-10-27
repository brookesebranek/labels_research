bedgraph_contig <-fread(bedgraph_data, col.names = c( "chrom" ,"chromStart", "chromEnd", "coverage"))   
bedgraph_window<-data.table::setkey(bedgraph_contig, chrom, chromStart, chromEnd)
window_and_bedgraph<-data.table::foverlaps(bedgraph_window, window, nomatch=NULL)

bedgraph4303_table<-data.table::fread(interesting_bedfile, col.names =c( "chrom" ,"chromStart", "chromEnd", "status", "mean"))
interesting_errors<-PeakError::PeakError(bedgraph4303_table[status=="peak"], labels_track)
bedgraph4303_table <-fread(segment, col.names = c( "chrom" ,"chromStart", "chromEnd", "status", "penalty"))
bedgraph4303_window<-data.table::setkey(bedgraph4303_table, chrom, chromStart, chromEnd)
labels_window<-data.table::setkey(labels, chrom, chromStart, chromEnd)
segment_window<-data.table::setkey(segment_table, chrom, chromStart, chromEnd)
window_and_labels<-data.table::foverlaps(labels_window, window)
window_and_segment<-data.table::foverlaps(segment_window, window)
window_and_bedgraph4303<-data.table::foverlaps(bedgraph4303_window, window)

ggplot()+geom_segment(aes(x=chromStart, y=chrom, xend=chromEnd, yend=chrom), data = segment_window)
ggplot()+geom_rect(aes(x =chromStart, y =chromEnd, xmin=15200000, xmax=15300000, ymin=0, ymax=20), data =labels_window)
ggplot()+geom_line(aes(x=chromStart, y=chromEnd), data =bedgraph4303_window)

ggplot()+geom_segment(aes(x=chromStart, y=chrom, xend=chromEnd, yend=chrom), data = segment_window)+ggplot()+geom_rect(aes(x =chromStart, y =chromEnd, xmin=15200000, xmax=15300000, ymin=0, ymax=20), data =labels_window)


