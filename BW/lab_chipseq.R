library(rtracklayer)
library(GenomicAlignments)
library(Rsamtools)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPQC)


bamfil = "BAM/00_0CWB_01S6JHU_LN95-Pooled_Input_hg38_i85.bam"
bedfil = "BW/00_0CWB_01S6JHU_LN95-Pooled_Input_hg38_i85_uniqnorm_signal.bw"

bam <- readGAlignments(bamfil)
peaks <- import.bw(con = bedfil);

# Quality metrics with ChIPQC
BlackListFile <- ("BW/blacklist/ENCFF001TDO.bed")
exp <- ChIPQCsample(bamfil,peaks = bedfil, annotation = "hg38",blacklist=BlackListFile, verbose = FALSE)


# Generate Ideogram
starting =   13121902 
ending = 13801106
ideogram <- IdeogramTrack("chrY", "hg38")
axis <- GenomeAxisTrack()
tx <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, chromosome = "chrY", start = starting, end = ending, name = "Genes")

# Creating data tracks - these are pairs (supposedly numeric, genome interval), with other useful things we care about defined.
# I am guessing it does filtering based on the other arguments by itself.
peaksDT <- DataTrack(peaks, chromosome = "chrY", start = starting, end = ending, name = "Peaks", type = "h")
peaksDT <- AnnotationTrack(bamfil, genome = "hg38", chromosome = "chrY", start = starting, end = ending, name = "Reads", type = "")

#covDT <- DataTrack(coverage(reads)[["chrY"]], chromosome = "chrY", start = starting, end = ending, name = "Coverage", type = "h") - Figure how to plot coverage
blacklist_go <- AnnotationTrack(blacklist, name = "Blacklist", id = blacklist@elementMetadata$name)
plotTracks(list(ideogram, axis, peaksDT, tx), chromosome = "chrY", from = starting, to = ending, featureAnnotation = "id")
