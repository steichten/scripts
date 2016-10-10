#making gene to gene TAIR10 annotation
#SRE
#oct11th, 2016
#

mkdir $(date +"%Y-%m-%d")_TAIR10_gene_to_gene_annotation
cd *_TAIR10_gene_to_gene_annotation

#get TAIR10 gene annotations
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff



#make bed file of all TAIR10 genes with proper strand information


###########################
R

getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

gene=gffRead('TAIR10_GFF3_genes.gff')

#I am subsetting to annotated 'gene's which there are 28,775 in total for TAIR10. This may be modified if we are looking for other things.
gene=subset(gene,gene$feature=='gene')
gene$Name=getAttributeField(gene$attributes, 'Name')
gene$ID=getAttributeField(gene$attributes, 'ID')

gene.out=gene[,c('seqname','start','end','Name','score','strand')]
write.table(gene.out,'TAIR10_genes.bed',sep='\t',row.names=F,col.names=F,quote=F)
quit()
n

########################


#ID closest gene to each gene
sort -k1,1 -k2,2n TAIR10_genes.bed > TAIR10_genes.sorted.bed

#flags for bedtools (v2.25.0)
# -N cannot match same name of gene (i.e. you don't match yourself)
# -iu ignore upstream features
# -s require feature on same strand
# -D a report distance in relationship to the orientation of file 'a'

bedtools closest -N -iu -s -D a -a TAIR10_genes.sorted.bed -b TAIR10_genes.sorted.bed > TAIR10_gene_to_gene.samestrand.bed


#parse output and create final annotation file

########################

R

input=read.delim('TAIR10_gene_to_gene.samestrand.bed',head=F)
colnames(input)=c('chr1','start1','stop1','gene1','score1','strand1','chr2','start2','stop2','gene2','score2','strand2','distance')

input.flt <- input[,c('chr1','start1','stop1','gene1','strand1','chr2','start2','stop2','gene2','strand2','distance')]
input.flt$chr2 <- ifelse(input.flt$chr2=='.',NA,as.character(input.flt$chr2))
input.flt$chr2 <- as.factor(input.flt$chr2)
input.flt$start2 <- ifelse(input.flt$start2=='-1',NA,input.flt$start2)
input.flt$stop2 <- ifelse(input.flt$stop2=='-1',NA,input.flt$stop2)
input.flt$gene2 <- ifelse(input.flt$gene2=='.',NA,as.character(input.flt$gene2))
input.flt$gene2 <- as.factor(input.flt$gene2)
input.flt$strand2 <- ifelse(input.flt$strand2=='.',NA,as.character(input.flt$strand2))
input.flt$strand2 <- as.factor(input.flt$strand2)
input.flt$distance <- ifelse(is.na(input.flt$chr2)==T,NA,input.flt$distance)

input.flt$flag <- ifelse(input.flt$distance==0,'overlap','nonoverlap')
#there are 568 records in which the distance equals 0. This indicates overlapping gene annotations on the same strand (???) not sure what that is all about.

write.table(input.flt,'TAIR10_gene_to_gene.samestrand.anno',sep='\t',row.names=F,quote=F)

quit()
n

#


