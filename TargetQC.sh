for i in F 6 Bo
do
echo "
library(\"TarSeqQC\")
library(\"BiocParallel\")
pdf(\"$i.QC.pdf\")
bedFile<-\"ampl.hg19.bed\"
bamFile<-\"$i.sorted.bam\"
fastaFile<-\"/BIC/Pipeline/SMTONE/database/hg19/hg19.fa\"
BPPARAM<-bpparam()
myPanel<-TargetExperiment(bedFile, bamFile, fastaFile, feature=\"amplicon\", attribute=\"coverage\", BPPARAM=BPPARAM)
setFeature(myPanel)<-\"amplicon\"
setAttribute(myPanel)<-\"coverage\"
scanBamP<-ScanBamParam()
bamWhich(scanBamP)<-getBedFile(myPanel)
setScanBamP(myPanel)<-scanBamP
setPileupP(myPanel)<-PileupParam(max_depth=2000)
#setFeaturePanel(myPanel)<-buildFeaturePanel(myPanel, BPPARAM)
setGenePanel(myPanel)<-summarizePanel(myPanel, BPPARAM)
g<-plotAttrExpl(myPanel,level=\"feature\",join=TRUE, log=FALSE, color=\"blue\")
g
plotMetaDataExpl(myPanel, \"length\", log=FALSE, join=FALSE, color=\"blueviolet\")
plotMetaDataExpl(myPanel, \"gene\", abs=FALSE)
readFrequencies(myPanel)
plotInOutFeatures(readFrequencies(myPanel))

attributeThres<-c(0,1,50,200,500,1000,Inf)
plot(myPanel, attributeThres=attributeThres, chrLabels =TRUE)

g<-plotFeatPerform(myPanel, attributeThres, complete=TRUE, log=FALSE,featureLabs=TRUE, sepChr=TRUE, legend=TRUE)
g
biasExploration(myPanel, source=\"gc\", dens=TRUE)
summaryIntervals(myPanel, attributeThres)
plotAttrPerform(myPanel, attributeThres)

getLowCtsFeatures(myPanel, level=\"gene\", threshold=50)
for (gene in c(\"DLG1\",\"GFAP\",\"GRB7\",\"NCKAP1L\",\"SLFN14\",\"SRCIN1\",\"THSD1\")){
        g<-plotGeneAttrPerFeat(myPanel, geneID=gene)
        g<-g+theme(title=element_text(size=16), axis.title=element_text(size=16),legend.text=element_text(size=14))
        g
        g<-plotGeneAttrPerFeat(myPanel, geneID=gene, overlap = TRUE, level=\"both\")
        g
}

" >tmp/$i.R
R CMD BATCH tmp/$i.R
done
