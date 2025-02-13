for(track.contig.i in 1:nrow(biggest_contigs_tracks)){
        current_contig_track<-biggest_contigs_tracks[track.contig.i]
        print(current_contig_track)
        track_slash<-sub("_", "/", current_contig_track$track)
        in_big_wig<-paste0(
                "https://rcdata.nau.edu/genomic-ml/PeakSegFPOP/labels/H3K4me3_TDH_ENCODE/samples/",
                track_slash, 
                "/coverage.bigWig")
        contigs.bedgraph <-paste0("H3K4me3_TDH_ENCODE/", "samples/", track_slash,"/problems", "/chr10:60000-17974675", "/coverage", ".bedgraph")
        sample_type_directory<-dirname(contigs.bedgraph)
        dir.create(sample_type_directory, showWarnings = F)
        PeakSegPipeline::bigWigToBedGraph(
                in.bigWig = in_big_wig, 
                out.bedGraph = contigs.bedgraph, 
                chrom = current_contig_track$chrom, 
                start = current_contig_track$contig_start, 
                end = current_contig_track$contig_end)
}
