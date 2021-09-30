data.table::fread("hg19_contigs.bed")
contigs<- data.table::fread("hg19_contigs.bed", col.names=c("chrom", "contig_start", "contig_end"))
labels<-data.table::fread("https://raw.githubusercontent.com/tdhock/peaklearner-api/main/labels/H3K4me3_TDH_ENCODE/2021-08-24/labels.csv")
labels<-data.table::fread("https://raw.githubusercontent.com/tdhock/peaklearner-api/main/labels/H3K4me3_TDH_ENCODE/2021-09-21_14_55_10/labels.csv")
data.table::setkey(labels,chrom,chromStart,chromEnd)
data.table::setkey(contigs,chrom,contig_start,contig_end)
labels_and_contigs<-data.table::foverlaps(labels,contigs)
labels_per_contig_track <- labels_and_contigs[, .N, by = .(contig_start,contig_end, chrom, track)]
labels_per_contig_track[order(N)]
biggest_contigs_tracks<-labels_per_contig_track[N==max(N)] 


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
        data.table::fwrite (
                contig_labels[, .(chrom, chromStart, chromEnd, annotation=sub("noPeak", "noPeaks", annotation))],
                labels.bed,sep = "\t")
        PeakSegPipeline::bigWigToBedGraphNoGaps(
                in.bigWig = in_big_wig, 
                out.bedGraph = coverage.bedGraph, 
                chrom = current_contig_track$chrom, 
                start = current_contig_track$contig_start, 
                end = current_contig_track$contig_end)
        PeakSegPipeline::problem.target(contig_directory)
}
