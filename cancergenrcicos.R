library("RCircos")
data("UCSC.HG38.Human.CytoBandIdeogram")
cyto.info = UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
chr_order = unique(cyto.info$Chromosome)

cyto.info = UCSC.HG38.Human.CytoBandIdeogram
cyto.info$Name = NA
cyto.info$Stain = NA
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

chr_order = unique(cyto.info$Chromosome)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

ideo = RCircos.Get.Plot.Ideogram()
ideo$BandColor = 'salmon'
num = which(ideo$Chromosome == 'chrX')
ideo[num, 'BandColor'] = 'chartreuse'

num = which(ideo$Chromosome == 'chrY')
ideo[num, 'BandColor'] = 'purple'

RCircos.Reset.Plot.Ideogram(ideo)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

num = which(ideo$Chromosome == 'chr1')
ideo[num, 'ChrColor'] = 'goldenrod2'

RCircos.Reset.Plot.Ideogram(ideo)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

library(biomaRt)

mat = read.csv('C:/Users/awahi/Desktop/gnjnerror/deg.csv', row.names = 1)
#df= data.frame(mat) 
library(org.Hs.eg.db)
library(AnnotationDbi)

head(mat)
dim(mat)
a=mat$gene_id
mat$gene_symbol<-mapIds(org.Hs.eg.db,keys=a,keytype = "ENSEMBL",column = "SYMBOL")
head(mat)
library(dplyr)
data=mat %>% select(gene_symbol,everything())
head(data)
library(tidyr)
data1=data %>% drop_na()
head(data1)
dim(data1)
df=data1[-2]
head(df)
#drop duplictes
#genes = c('01-Mar', '02-Mar', 'ALDOA', 'ALG1L9P', 'ARHGAP11B', 'ATXN7', 'BMS1P4', 'CCDC39', 'COG8', 'EMG1', 'GOLGA8M', 'HSPA14', 'LINC00484', 'LINC01238', 'LINC01422', 'LINC01481', 'MATR3', 'POLR2J3', 'POLR2J4', 'RF00001', 'RF00003', 'RF00012', 'RF00015', 'RF00017', 'RF00019', 'RF00066', 'RF00156', 'RF00181', 'RF00322', 'RF00432', 'RF00554', 'RF00561', 'RF00564', 'RF01210', 'RF01225', 'RF01241', 'RF02271', 'TBCE', 'TMSB15B', 'ZNF883')
#num = which(df[,1] %in% genes)
#df = df[-num, ]

rownames(df) = df[,1]
df = df[,-1]


saveRDS(df,"deg.rds")

m = useMart('ensembl', dataset='hsapiens_gene_ensembl')

coords = getBM(attributes=c('chromosome_name', 'start_position', 
                            'end_position', 'hgnc_symbol'),
               filters = c('hgnc_symbol'),
               values = list(rownames(df)),
               mart = m)

write.csv(coords, file = 'coords.csv')
coords$chromosome_name = paste0('chr', coords$chromosome_name)
coords$chromosome_name = factor(coords$chromosome_name, levels = chr_order)
num = which(is.na(coords$chromosome_name))
coords = coords[-num, ]

up = which((df$pvalue<1) &
             (df$log2FC>-2))
upmat = df[up, ]

num = which(coords$hgnc_symbol %in% rownames(upmat))
coords1 = coords[num, ]
library(RCircos)
RCircos.Gene.Name.Plot(coords1, name.col=4, track.num = 2, side = "in",
                       is.sorted = F)


genes = intersect(rownames(df), coords$hgnc_symbol)

mat2 = df[genes, ]
dff = cbind.data.frame(row.names(df), df[, c(1,2,3)])
colnames(dff)[1] = 'hgnc_symbol'

data = merge(coords, dff, by = 'hgnc_symbol', all.x = T)
library(tidyr)
data=data %>% drop_na()
data=data[,c(2,3,4,5,6,1,7)]

data = data[, c('chromosome_name', 'start_position',
                'end_position', 
                'meanTumor', 'meanControl', 'hgnc_symbol','log2FC')]
library(tibble)

RCircos.Heatmap.Plot(data, data.col = 7, track.num = 6, side = "in",
                     min.value = -3, max.value = 3.5, genomic.columns = 3,
                     is.sorted = F)


RC.param = RCircos.Get.Plot.Parameters()
RC.param['heatmap.color'] = "GreenWhiteRed"
RCircos.Reset.Plot.Parameters(RC.param)

RCircos.Heatmap.Plot(data, data.col = 7, track.num = 10, side = "in",
                     min.value = -2, max.value = -1,
                     is.sorted=F)