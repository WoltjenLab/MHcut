args = commandArgs(TRUE)
CVFILE = args[1]
DBFILE = args[2]
GENCFILE = args[3]
OUTFILE = args[4]

## TO DEBUG
## CVFILE = 'clinvar-grch38-deletion.tsv.gz'
## DBFILE = 'dbsnp-nsf-div-del.tsv'
## OUTFILE = 'testMerge.tsv'
## GENCFILE = '/Users/jeanmonlong/Documents/GenomeCenter/data/gencode.v28.basic.annotation.gtf.gz'

## Read ClinVar file
message('Reading the ClinVar file...')
cv = read.table(CVFILE, as.is=TRUE, header=TRUE, sep='\t', quote='')
cv.coords = paste(cv$chr, cv$start, cv$end, sep='-')

## Import gencode annotation
message('Reading the Gencode file...')
genc = read.table(GENCFILE, as.is=TRUE, sep='\t')
colnames(genc) = c('chr','source','type','start','end','score', 'strand', 'phase', 'attributes')
suppressMessages(library(GenomicRanges))
gene.gr = makeGRangesFromDataFrame(subset(genc, type=='gene'))
exon.gr = makeGRangesFromDataFrame(subset(genc, type=='exon'))
utr.gr = makeGRangesFromDataFrame(subset(genc, type=='UTR'))

## Read dbSNP file by chunks
message('Reading the dbSNP file...')
CHUNK.SIZE = 1e5
con = file(DBFILE, 'r')
lines = readLines(con, n=1)
db.head = unlist(strsplit(lines, '\t'))
final.cols = c(db.head, setdiff(colnames(cv), c('chr', 'start', 'end', 'RS')))
firstChunk = TRUE
while(length((lines = readLines(con, n=CHUNK.SIZE)))>0){
  db = read.table(textConnection(lines), as.is=TRUE, sep='\t', quote='')
  colnames(db) = db.head
  ## Merge by rsID first
  rsol = which(db$RS %in% cv$RS)
  if(length(rsol)){
    dbrs = merge(db[rsol,], cv[, setdiff(colnames(cv), c('chr', 'start', 'end'))])
    db = db[-rsol,]
    rs.notol = which(!(cv$RS %in% dbrs$RS))
    cv = cv[rs.notol,]
    cv.coords = cv.coords[rs.notol]
  }
  ## Then merge by coord
  db.coords = paste(db$chr, db$start, db$end, sep='-')
  coordol = which(db.coords %in% cv.coords)
  if(length(coordol)>0){
    dbcoord = merge(db[coordol, setdiff(colnames(db), c('RS'))], cv)
    db = db[-coordol,]
    coord.notol = which(!(cv.coords %in% db.coords))
    cv = cv[coord.notol,]
    cv.coords = cv.coords[coord.notol]
  }
  ## Add columns from ClinVar
  for(coln in setdiff(colnames(cv), c('chr', 'start', 'end', 'RS'))){
    db[,coln] = NA
  }
  db = db[, final.cols]
  ## Bind the three subsets
  if(length(rsol)>0){
    db = rbind(db, dbrs[, final.cols])
  }
  if(length(coordol)>0){
    db = rbind(db, dbcoord[, final.cols])
  }
  ## Gene annotation
  db.gr = makeGRangesFromDataFrame(db)
  db$geneloc = 'intergenic'
  db$geneloc[overlapsAny(db.gr, gene.gr)] = 'intronic'
  db$geneloc[overlapsAny(db.gr, utr.gr)] = 'UTR'
  db$geneloc[overlapsAny(db.gr, exon.gr)] = 'exonic'
  ##
  write.table(db, OUTFILE, quote=FALSE, sep='\t', append=!firstChunk, col.names=firstChunk, row.names=FALSE)
  firstChunk = FALSE
}
close(con)

## Write what remains in ClinVar
if(nrow(cv)>0){
  message('Adding remaining ClinVar variants')
  for(coln in setdiff(colnames(db), c('chr', 'start', 'end', 'RS'))){
    cv[,coln] = NA
  }
  cv = cv[, final.cols]
  ## Gene annotation
  cv.gr = makeGRangesFromDataFrame(cv)
  cv$geneloc = 'intergenic'
  cv$geneloc[overlapsAny(cv.gr, gene.gr)] = 'intronic'
  cv$geneloc[overlapsAny(cv.gr, utr.gr)] = 'UTR'
  cv$geneloc[overlapsAny(cv.gr, exon.gr)] = 'exonic'
  ##
  write.table(cv, OUTFILE, quote=FALSE, sep='\t', append=!firstChunk, col.names=firstChunk, row.names=FALSE)
}
