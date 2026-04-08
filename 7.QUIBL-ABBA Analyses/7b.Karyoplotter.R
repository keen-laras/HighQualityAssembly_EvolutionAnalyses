library(karyoploteR)
args = commandArgs(T)

lenFile = args[1]
datFile = args[2]
pdfFile = args[3]

pdf(pdfFile)

genomeLen = read.table(lenFile)
colnames(genomeLen) = c('Chr', 'Len')
n = nrow(genomeLen)
custom.genome <- toGRanges(data.frame(chr=genomeLen$Chr, start=rep(1, n), end=genomeLen$Len))

kp <- plotKaryotype(genome = custom.genome)

dat = read.table(datFile, h=T)
dat_discord = dat[dat$value == 0,]
dat_concord = dat[dat$value == 1,]

custom.genes = toGRanges(data.frame(chr=dat_discord$Chr, start=dat_discord$Start, end=dat_discord$End))

kp <- kpPlotDensity(kp, custom.genes, col='blue')
