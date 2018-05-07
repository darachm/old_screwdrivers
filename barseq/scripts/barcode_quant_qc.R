---
title: "Preliminary analysis and QC"
author: "darach miller"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
---

# HPC work

Previously, the fastq files of Bill Chirico's were downloaded onto
the NYU HPC mercer cluster, on dhm267's account. dhm267 went ahead
and ran BarNone with the following script, each time setting the 
`$MM` variable to different mismatches:

    #!/bin/bash
    #PBS -l nodes=1:ppn=1,walltime=72:00:00
    #PBS -N barNoneChirico
    #PBS -M dhm267@nyu.edu
    #PBS -m abe
    #PBS -e localhost:/scratch/dhm267/${PBS_JOBNAME}.${PBS_JOBID}.e
    #PBS -o localhost:/scratch/dhm267/${PBS_JOBNAME}.${PBS_JOBID}.o
    #by Darach Miller
    
    # When are we?
    
    DATE=$(date +%y%m%d:%H:%M:%S)
    
    # What is the experiment designator?
    
    EXP="chiricoR2"
    
    # Where are you working?
    
    DIR="/scratch/dhm267/"$EXP
    
    # Where is your output GZip'd FASTQ file?
    
    INPUTFASTQ="./lane8_NoIndex_L008_R1_001.fastq"
    
    # And let's try with these number of allowed mismatches
    
    MM=$MM
    
    # Let's begin
    
    module load barnone/intel
    
    cd $DIR
    
    echo "JobID is $PBS_JOBID, running on $DATE
    
    You are procing $EXP in $DIR with $MM mismatches
    
    Good luck
    "
    
    BarNone \
      -f fastq \
      --multiplexfile /home/dhm267/ref/demuxing_barseq_index.txt \
      --multiplexstart 1 --multiplexlength 5 \
      --tagstart 18 --taglength 3 \
      --start 21 \
      --mismatches $MM \
      --mismatchfile ${DIR}/${EXP}mismatch${MM}MM${DATE}.txt \
      --revisedcatalog ${DIR}/${EXP}revisedBarcodes${MM}MM${DATE}.txt \
      -p 100000 \
      $INPUTFASTQ \
      ${DIR}/${EXP}counts${MM}MM${DATE}.txt \
      /home/dhm267/ref/nislow_revised.txt
    
    # More argument info
    #	--uptag specifies what the uptag is, default "TCT"
    #	--downtag specifies what the uptag is, default "TAG"
    #	--tagstart tells BarNone where to look
    #	--taglength tells how long the tag indicator is
    #	--multiplexstart tells the start of that index
    #	--multiplexlength is of course the length of that index
    #	--multiplexfile is what maps the that index to your sample
    #	--mismatches tells how many mismatches to tolerate
    #	--mismatchfile is a nice record of the mismatches
    #	--start is where the barcode starts
    #	--length is it's length
    #	--revisedcatalog outputs a reannotation of mutant to its most common close
    #		barcode, and use that catalog
    #	-p reports every p reads
    #	-n just limits the run to first n reads, good for testing?
    #	-f specifies which input file type it is, namely txt or fastq
    #	arguments end with "infile outfile strain_barcode_file.txt"

On the first batch of runs, mismatches of 4 and 5 failed to complete
after 72 hours, due to timing out. 
Re-ran these on April 19th with 96 hour limits.  
On April 20th they were done.

# Local work

```{r,libs,cache=F}
library(tidyverse)
library(pcaMethods)
library(ggrepel)
set.seed(2345)
```

## Reading in

The counts tables, mismatch files, etc. were sftp'd to local laptop.
Bill's index was edited in localc to make it easier to use in R,
and saved as a csv. We read that in below.

Importantly, we keep track of samples 33-43, as those were test 
replicates for which it took way too long to sort, so Bill doesn't 
trust those.
We keep track, and later I'll show why we're excluding those
(not a lot of experimental variation, by PCA).

The index is convienently in order by row.

```{r index_read,cache=T}

index <- read_csv("../rawdata/repairScreenIndex.csv") 

head(index)

```


```{r data_read,cache=T}

rawdat <- tibble(FileName=list.files(path="../rawdata/",
    pattern="chirico.*counts.*txt",full.names=T)) %>%
  group_by(FileName) %>%
  mutate(NumberMismatches=sub(".*(\\d)MM.*","\\1",FileName),
    RawFile=list(read_tsv(FileName))) %>%
  unnest() %>% 
  gather(SampleName,Counts,-FileName,-NumberMismatches,-Strain) %>%
  separate(SampleName,c("SampleNumber","Tag"),sep="_") %>%
  mutate(SampleNumber=as.numeric(sub("Sample","",SampleNumber)))

save(rawdat,file="../tmp/rawData.RData")

```

## Which mismatch parameter to use?

So which mismatch tolerance to use? The assumption is that some
barcodes ain't perfect, and that Barnone will find the right one
within a mismatch parameter. But search too far and we'll pick
up strains that shouldn't be in our collection.

Which one to use?

Let's look at what happens if we expect it in the library or not.
Here's a list of all the strain supposed to be in the 
invitrogen diploid homo deletion collection library.
The `Homo_diploids_041902.txt` file is from Invitrogen.

```{r,realstrains,cache=T,warning=F}

includedStrains <- (tibble(RawLines=
    read_lines("../rawdata/Homo_diploids_041902.txt",skip=2))%>%
  separate(RawLines,c("record","ORF","strain1","strain2",
    "batch","plate","row","column"),"\\s") %>%
  filter(ORF!=""))$ORF

rawdat <- rawdat %>% mutate(Expected=Strain%in%includedStrains)

```

And below is a plot, on the left are strains we don't expect in our
library (hetero KO of essential genes in a haploid collection 
for example).

```{r,realstrainsPlot,cache=T,fig.height=10}

rawdat %>% ggplot()+
  aes(x=Counts,col=factor(SampleNumber))+
  facet_wrap(~NumberMismatches+Expected+Tag,scales="free_y",ncol=4)+
  scale_x_log10()+
  stat_bin(bins=50,geom="line",aes(y=..count..),position="identity")+
  guides(col=F)

```

There's a problem here, because some things on the left look
like there's authentically things we do not expect in the data.
That's odd. Maybe this is a different library? Need to check later.
I don't think this is an error, since we see it in the `0` 
mismatches analysis.

But we do see that with more mismatches we start to pick up the tail
of 1-off detected strains, so we definitely don't want above `3`
mismatches. Let's look another way below.

Each plot below shows on one axis the total sum of counts for a strain
with a certain mismatch parameter in BarNone (so 3MM means three
mismatches are allowed), and the next highest mismatch on the y.
So if they're the same value, it falls on the diagonal.

If increasing the mismatch gets more counts for that Strain (good)
then the point goes to the right. If that mismatch loses counts,
the point goes down from the diagonal. This would be a case where
the mismatch is so high that some other strain is canabalizing 
counts from that strain.

Above diagonal good, below diagonal bad. Let's begin.

```{r,sumCounts,cache=T}

countsPerMM <- rawdat %>% 
  group_by(Strain,Tag,NumberMismatches) %>%
  summarize(TotalCounts=sum(Counts)) %>%
  spread(NumberMismatches,TotalCounts)

g<-countsPerMM %>% ggplot()+
  geom_point()+facet_wrap(~Tag)+
  scale_x_log10()+scale_y_log10()

g+aes(x=`0`,y=`1`)
g+aes(x=`1`,y=`2`)
g+aes(x=`2`,y=`3`)
g+aes(x=`3`,y=`4`)
g+aes(x=`4`,y=`5`)
g+aes(x=`5`,y=`6`)

```

Well, `3` looks pretty good, but let's be a bit conservative and
run with `2` mismatches. We don't get much more with `3`.
This is based on UPTAGs (the DOWNTAGs always seem odd...) .

Next, how about our total counts in each sample?

```{r,canaBetweenSamples,cache=T,fig.height=15}

tmp <- full_join(rawdat,index%>%mutate(SampleNumber=index),
  by="SampleNumber") %>%
  mutate(WasUsed=!is.na(Replicate)) %>%
  group_by(SampleNumber,NumberMismatches,Tag,WasUsed) %>%
  summarize(TotalCounts=sum(Counts))

tmp %>%
  ggplot()+
  facet_grid(NumberMismatches~Tag,scales="free")+
  aes(x=SampleNumber,y=TotalCounts,col=WasUsed)+
  geom_point()+
  geom_text_repel(data=tmp%>%filter(WasUsed,TotalCounts<1e5),
    aes(label=SampleNumber))+
  scale_y_log10(limits=c(100,NA),breaks=10^(2:7))+
  scale_color_discrete("Is a sample\nnumber we\nmeasured")+
  theme(axis.text.x=element_text(angle=90))+
  ylab("Total counts for index")

```

We see that samples 13 and 25 uptag are odd ones. Look like they 
failed, so we should exclude them.

```{r,totals,cache=T}

full_join(rawdat,index%>%mutate(SampleNumber=index),by="SampleNumber") %>%
  mutate(InDataSet=!is.na(Replicate)) %>%
  filter(!(SampleNumber%in%c(13,25))) %>%
  group_by(SampleNumber,Tag,InDataSet,NumberMismatches) %>%
  summarize(TotalCounts=sum(Counts)) %>%
  group_by(Tag,InDataSet,NumberMismatches) %>%
  summarize(MeanTotalCounts=mean(TotalCounts)) %>%
  ggplot()+
  aes(x=factor(NumberMismatches):factor(Tag),fill=InDataSet,
      y=MeanTotalCounts) +
    geom_point()+
    scale_y_log10(breaks=c(0.2*10^c(1:7),0.4*10^c(1:7),0.6*10^c(1:7),0.8*10^c(1:7),10^c(1:7)))+
    theme(axis.text.x=element_text(angle=90))+
    facet_wrap(~InDataSet,ncol=1,scales="free_y")

```

More mismatches, more detection. Most of the wanted gains happen by
`2MM`, then it's more noise after that.

```{r,lookingatMMdist,cache=T,fig.height=15}
full_join(rawdat,index%>%mutate(SampleNumber=index),
  by="SampleNumber") %>% filter(!is.na(Replicate)) %>%
  ggplot()+aes(x=Counts,col=factor(SampleNumber):factor(Tag))+
  stat_bin(geom="line",bins=30,position="identity")+
  scale_x_log10()+
  guides(col=F)+
  facet_wrap(~NumberMismatches+Expected+(SampleNumber>32),
    ncol=4,scales="free")
```

Looks like there's still plenty of counts on things that we don't
expect to be in our library.
With 1MM, it looks like very few things that we don't expect in our
library are counted. But, there still is!

Let's go with 2 mismaches for now, and drop sample 13 and 25 UP.

```{r pick_one,cache=T}
useMM <- "2"
datar <- full_join(rawdat,
  index%>%mutate(SampleNumber=index),by="SampleNumber") %>% 
    filter(NumberMismatches==useMM,!is.na(Replicate),Expected) %>%
  filter(!(SampleNumber%in%c(13,25)))

tmp <- datar %>% 
  group_by(Damaged,Repaired,Expanded,Gate,Replicate,Tag) %>%
  summarize(TotalCounts=sum(Counts))
datar <- full_join(datar,tmp,
  by=c("Damaged","Repaired","Expanded","Gate","Replicate","Tag")) %>%
  mutate(ProportionCounts=Counts/TotalCounts)

```


```{r,checkingNumberdetected,cache=T}

datar %>% group_by(Damaged,Repaired,Expanded,Gate,Replicate,Tag) %>%
  summarize(NumDetected=sum(Counts>10)) %>%
  ggplot()+
    aes(y=NumDetected,x=factor(Gate):factor(Replicate):factor(Tag))+
    geom_point()+facet_wrap(~Expanded+Damaged+Repaired,scales="free_x")+
    theme(axis.text.x=element_text(angle=90))+
    scale_y_continuous(limits=c(0,4000))
```

We detect more in the ones we expect.

```{r,checkingInfocontent,cache=T}

datar %>% group_by(Damaged,Repaired,Expanded,Gate,Replicate,Tag) %>%
  mutate(Info=ProportionCounts*log(ProportionCounts)) %>%
  summarize(TotalInfo=-sum(Info[!is.na(Info)&is.finite(Info)])) %>%
  ggplot()+
    aes(y=TotalInfo,x=factor(Gate):factor(Replicate):factor(Tag))+
    geom_point()+facet_wrap(~Expanded+Damaged+Repaired,scales="free_x")+
    theme(axis.text.x=element_text(angle=90))
```

We're using information content as a metric for how diverse and 
complex the sample is on a strain by strain basis. Those expanded
populations look very selected, which may be informative. To
continue on this, let's turn to PCA.

## Let's look at the variation with PCA, see if there's weird samples

Do replicates look like each other?

Do the major PCA axis capture some of the variation?

```{r pca,cache=T}

ddatar <- data.frame(datar %>% ungroup() %>%
    unite(SampleName,SampleNumber,Tag,Damaged,Repaired,Gate,Expanded,Replicate) %>% 
    select(Strain,SampleName,Counts) %>%
    spread(SampleName,Counts))
ddatar[,-1] <- apply(ddatar[,-1],2,as.numeric)
rownames(ddatar) <- ddatar[,"Strain"]
ddatar <- ddatar[,-1]
colnames(ddatar) <- sub("X","",colnames(ddatar))
ddatar <- as.matrix(ddatar)

pcAnalysis <- pca(t(ddatar),method="svd",scale="uv",center=T,nPcs=10)
plotdf <- scores(pcAnalysis)
plotdf <- data.frame(plotdf,reshape2::colsplit(rownames(plotdf),"_",
  c("SampleNumber","Tag","Damaged","Repaired","Gate","Expanded","Replicate")))

g<-ggplot(plotdf)+geom_point()+
  theme(legend.position="bottom")+
  aes(col=factor(Expanded):factor(Gate):factor(Damaged):
      factor(Repaired):factor(SampleNumber>32),
    label=SampleNumber)+
  facet_wrap(~Tag+(SampleNumber>32))+
  guides(col=guide_legend(title="Sample",nrow = 4))
g+aes(x=PC1,y=PC2)
```

Also note we still haven't excluded the untrusted 
longer-than-reasonably sorted samples (33-43). How do they look?
On the right are those quesitonable samples, left is the rest.

```{r,pcaplotz,cache=T}

g+aes(x=PC1,y=PC2)+geom_text_repel(size=3)
g+aes(x=PC3,y=PC2)+geom_text_repel(size=3)
g+aes(x=PC3,y=PC4)+geom_text_repel(size=3)

```

This looks like PC1 captures a lot of the variation based on 
selection, in that the expanded ones look really far right.
PC2 looks like DN/UP variation, and the other PCs look like they
capture variation reproducibly between up and down tags.
The replicate to replicate variation of the expanded populations
makes me think they're just jackpotting like crazy. These will be 
hard to incoproate.

The untrusted samples on the right show less variation, so I think
it's right to exclude them.
Then we return to the PCA to take a good harder look within the Tags.

```{r,exclude,cache=T}

datar <- datar %>% filter(SampleNumber<33)

ddatar <- data.frame(datar %>% ungroup() %>%
    unite(SampleName,SampleNumber,Tag,Damaged,Repaired,Gate,Expanded,Replicate) %>% 
    select(Strain,SampleName,Counts) %>%
    spread(SampleName,Counts))
ddatar[,-1] <- apply(ddatar[,-1],2,as.numeric)
rownames(ddatar) <- ddatar[,"Strain"]
ddatar <- ddatar[,-1]
colnames(ddatar) <- sub("X","",colnames(ddatar))
ddatar <- as.matrix(ddatar)

```


```{r,pcaUP,cache=T}

pcAnalysisUP <- pca(t(ddatar[,grepl("UP",colnames(ddatar))]),
  method="svd",scale="uv",center=T,nPcs=10)
plotdfUP <- scores(pcAnalysisUP)
plotdfUP <- data.frame(plotdfUP,reshape2::colsplit(rownames(plotdfUP),"_",
  c("SampleNumber","Tag","Damaged","Repaired","Gate","Expanded","Replicate")))

g<-ggplot(plotdfUP)+geom_point()+
  theme(legend.position="bottom")+
  aes(col=factor(Expanded):factor(Gate):factor(Damaged):
      factor(Repaired),
    label=SampleNumber)+
  guides(col=guide_legend(title="Sample",nrow = 4))

g+aes(x=PC1,y=PC2)+geom_text_repel(size=3)
g+aes(x=PC3,y=PC2)+geom_text_repel(size=3)
g+aes(x=PC3,y=PC4)+geom_text_repel(size=3)

```

```{r,pcaDOWN,cache=T}

pcAnalysisDOWN <- pca(t(ddatar[,grepl("DOWN",colnames(ddatar))]),
  method="svd",scale="uv",center=T,nPcs=10)
plotdfDOWN <- scores(pcAnalysisDOWN)
plotdfDOWN <- data.frame(plotdfDOWN,reshape2::colsplit(rownames(plotdfDOWN),"_",
  c("SampleNumber","Tag","Damaged","Repaired","Gate","Expanded","Replicate")))

g<-ggplot(plotdfDOWN)+geom_point()+
  theme(legend.position="bottom")+
  aes(col=factor(Expanded):factor(Gate):factor(Damaged):
      factor(Repaired),
    label=SampleNumber)+
  guides(col=guide_legend(title="Sample",nrow = 4))

g+aes(x=PC1,y=PC2)+geom_text_repel(size=3)
g+aes(x=PC3,y=PC2)+geom_text_repel(size=3)
g+aes(x=PC3,y=PC4)+geom_text_repel(size=3)

```

Huh. Well, let's take out sample 1 from the uptags

```{r,finalsubset,cache=F}
nrow(datar)
datar <- datar %>% filter(!(SampleNumber=="Sample1"&Tag=="UP")) %>%
  ungroup() %>%
  select(-FileName,-NumberMismatches,-index,-SampleNum,-TotalCounts,-ProportionCounts,-Expected)
nrow(datar)
```

# Writing QCd output for analysis

```{r output_qcd_counts,cache=F}
dim(datar)
sort(unique(factor(datar$SampleNumber):factor(datar$Tag)))
save(datar,file="../tmp/chirico2MMcountsMelty.RData")
write.csv(file="../tmp/chirico2MMcountsWide.csv",
  x=reshape2::dcast(datar,
    Strain~SampleNumber+Damaged+Repaired+Gate+Expanded+Replicate+Tag,
    value.var="Counts")
  ) 
```

---
```{r}
sessionInfo()
```
