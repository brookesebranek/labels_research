library(data.table)
meta_data_filename<-"https://rcdata.nau.edu/genomic-ml/2022-03-24-brooke-data/meta_data.csv"
meta_data_table<-fread(meta_data_filename)  
interesting_samples_table<-meta_data_table[antibody=="H3K4me3"]
interesting_samples_table[,status:=ifelse(disease=="None","healthy","cancer")]
interesting_samples_table[, sampleID:=paste0(status,.I),by=status]

show_coverage_list<-list()
for(sample_number in 1:nrow (interesting_samples_table)){
    current_sample_row<-interesting_samples_table[sample_number,]
    print(current_sample_row)
    in_big_wig<-current_sample_row[, paste0("https://rcdata.nau.edu/genomic-ml/2022-03-24-brooke-data/",
      "ihec.chipseq.ihec-chipseq-containerv1.1.4.", epirr_id, ".", uuid, ".raw.bigwig")]
    setwd("/Users/brookesebranek/Desktop/labels_research./H3K4me3_analysis")
    disease_comparison_directory<-current_sample_row[, paste0("/Users/brookesebranek/Desktop/labels_research./", 
                                      antibody, "/", sampleID, "/chr11:118,186,295-118,270,296")]
    coverage.bedGraph <-paste0(disease_comparison_directory, "/coverage.bedGraph")
    dir.create(disease_comparison_directory, recursive = T)
      PeakSegPipeline::bigWigToBedGraph(
        in.bigWig = in_big_wig, 
        out.bedGraph = coverage.bedGraph, 
        chrom = "chr11", 
        start = 118186295, 
        end = 118270296 )
      coverage_table<-fread(coverage.bedGraph, col.names = c( "chrom" ,"chromStart", "chromEnd", "matches"))
      show_coverage_list[[paste0(sample_number)]]<-current_sample_row[, data.table(sampleID, coverage_table)]
      }
install.packages('R.utils')
genes_table<-fread("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz", col.names = c("name", "chrom",
"strand", "txStart", "txEnd", "cdsStart","cdsEnd", "exonCount", "exonStarts", "exonEnds", "proteinID", "alignID"))
chrom_11_genes_table<-genes_table[chrom=="chr11"]
JAML1_gene_table<-chrom_11_genes_table[name=="ENST00000356289.10"]
ggplot_table<-do.call(rbind, show_coverage_list)
library(ggplot2)
#read in four files as data tables
gg<-ggplot2::ggplot()+geom_path(aes(x=chromStart, y=matches), color="gray50", data=ggplot_table)+facet_grid(sampleID~., scales = "free"
            )+geom_segment(aes(x=txStart, xend=txEnd, y=250, yend=250), color="green3",  lwd=3, data=JAML1_gene_table)
print(gg)
