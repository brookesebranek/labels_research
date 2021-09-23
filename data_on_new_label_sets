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
