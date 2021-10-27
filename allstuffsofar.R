(contigs<- data.table::fread("hg19_contigs.bed", col.names=c("chrom", "contig_start", "contig_end")))
labels<-data.table::fread("https://raw.githubusercontent.com/tdhock/peaklearner-api/main/labels/H3K4me3_TDH_ENCODE/2021-09-21_14_55_10/labels.csv")
data.table::setkey(labels,chrom,chromStart,chromEnd)
data.table::setkey(contigs,chrom,contig_start,contig_end)
labels_and_contigs<-data.table::foverlaps(labels,contigs)
labels_per_contig_track <- labels_and_contigs[, .N, by = .(contig_start,contig_end, chrom, track)]
labels_per_contig_track[order(N)]


library(data.table)
labelcounts <- fread ("wc -l Downloads/data/*/samples/*/*/problems/*/labels.bed",col.names = c("number_of_labels" ,"file"))
labelcounts[order(number_of_labels)]
int.pattern <- list("[0-9]*", as.integer)
start.end.pattern <- list(":", start=int.pattern, "-", end=int.pattern)
label_counts_with_positions = nc::capture_first_df(labelcounts, file=start.end.pattern, nomatch.error = FALSE)
label_counts_with_positions[, bases := end-start]
 str(label_counts_with_positions)
 label_counts_with_positions[, density := number_of_labels/bases]
 str(label_counts_with_positions)
 label_counts_with_positions[order(density)]
 label_counts_with_positions[order(number_of_labels)]
 
 
 
 
(contigs<- data.table::fread("hg19_contigs.bed", col.names=c("chrom", "contig_start", "contig_end")))
bad_labels<-data.table::fread("https://raw.githubusercontent.com/tdhock/peaklearner-api/main/labels/H3K4me3_TDH_ENCODE/2021-08-24/labels.csv")
bad_labels<-data.table::fread("https://raw.githubusercontent.com/tdhock/peaklearner-api/main/labels/H3K4me3_TDH_ENCODE/2021-09-21_14_55_10/labels.csv")
labels<-bad_labels[, .(chrom, chromStart, chromEnd, track, annotation=sub("noPeak", "noPeaks", annotation))]
data.table::setkey(labels,chrom,chromStart,chromEnd)
data.table::setkey(contigs,chrom,contig_start,contig_end)
labels_and_contigs<-data.table::foverlaps(labels,contigs)
labels_per_contig_track <- labels_and_contigs[, .N, by = .(contig_start,contig_end, chrom, track)]
labels_per_contig_track[order(N)]
biggest_contigs_tracks<-labels_per_contig_track[N==max(N)] 
better_contigs<-contigs[!grepl("_", chrom)]
ggplot()+geom_segment(aes(x=contig_start, y=chrom, xend= contig_end, yend=chrom),data=better_contigs)+geom_point(aes(x=chromEnd, y=chrom, color=annotation),data=labels)

for(track.contig.i in 1:nrow(biggest_contigs_tracks)){
        current_contig_track<-biggest_contigs_tracks[track.contig.i]
        contig_labels<-labels_and_contigs[current_contig_track, on=.NATURAL]
        print(current_contig_track)
        track_slash<-sub("_", "/", current_contig_track$track)
        in_big_wig<-paste0(
                "https://rcdata.nau.edu/genomic-ml/PeakSegFPOP/labels/H3K4me3_TDH_ENCODE/samples/",
                track_slash, 
                "/coverage.bigWig")
        contig_directory<-paste0("H3K4me3_TDH_ENCODE/samples/", track_slash,"/problems/chr10:60000-17974675")
        coverage.bedGraph <-paste0(contig_directory, "/coverage.bedGraph")
        labels.bed<-paste0(contig_directory, "/labels.bed")
        dir.create(contig_directory, showWarnings = F, recursive = T)
        data.table::fwrite (renamed_labels, 
                labels.bed,sep = "\t")
        PeakSegPipeline::bigWigToBedGraphNoGaps(
                in.bigWig = in_big_wig, 
                out.bedGraph = coverage.bedGraph, 
                chrom = current_contig_track$chrom, 
                start = current_contig_track$contig_start, 
                end = current_contig_track$contig_end)
        target_results<-PeakSegPipeline::problem.target(contig_directory)
        segment_files<-Sys.glob(paste0(contig_directory, "/*segments.bed"))
        for(segment in segment_files){
               segment_table <-fread(segment, col.names = c( "chrom" ,"chromStart", "chromEnd", "status", "mean"))
               area_to_view<-PeakError::PeakErrorChrom(segment_table[status=="peak"], renamed_labels)
               
        }
window<- data.table::data.table(chrom="chr10", window_start=15200000, window_end=15300000)  
data.table::setkey(window, chrom, window_start, window_end)
segment_window<-segment_table[chromStart>15200000 & chromEnd<15300000]
interesting_track<-"ENCFF280HQO"
interesting_cell_type<-"skeletalMuscle"
interesting_bedfile<-paste0(
        "H3K4me3_TDH_ENCODE/samples/",
        interesting_cell_type, 
        "/",
        interesting_track, 
        "/problems/chr10:60000-17974675/coverage.bedgraph_penalty=4303.13667437155_segments.bed")
labels_track<- labels[paste0(interesting_cell_type,
                             "_",
                             interesting_track),
                      on="track"]

bedgraph_data<- paste0(
        "H3K4me3_TDH_ENCODE/samples/",
        interesting_cell_type, 
        "/",
        interesting_track, 
        "/problems/chr10:60000-17974675/coverage.bedgraph")
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

usethis::create_github_token()
install.packages("usethis")
