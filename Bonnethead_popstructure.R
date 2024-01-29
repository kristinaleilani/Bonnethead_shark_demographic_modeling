library(dplyr)
library(vegan)
library(ggplot2)
library(plyr)
setwd("~/Documents/Bonnets") 

# Import sample names
samples=read.table("bams.qc") 
names(samples)[1]<-"Sample"
samples$Sample=paste(sub(".bam","",samples$Sample),sep="")
rownames(samples)<-samples$Sample

# Import IBS matrix
IBS=as.matrix(read.table("Bonnet2.ibsMat")) 
dimnames(IBS)=list(rownames(samples),rownames(samples))

# Initial clustering of samples by IBS
# plot dendrogram
hc=hclust(as.dist(IBS),"ave")
plot(hc) 
# plot PCoA
ord.all=capscale(as.dist(IBS)~1)
plot(ord.all,scaling=1) 


# Import metadata
meta=read.delim("Bonnets_meta2.txt")[1:127,]
env=join(samples, meta, by="Sample")
rownames(env)<-env$Sample
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}
env= env %>% mutate_each(funs(empty_as_na))
env$Maturity<-as.factor(env$Maturity)
env$Sex<-as.factor(env$Sex)
env$Maturity= factor(env$Maturity, levels = c("immature", "mature"))
env[,14] # Check levels of maturity: smallest to largest
env$Sex= factor(env$Sex, levels = c("M", "F"))
env[,5] # Check levels of sex: small(male) to large(female)



####---- Evaluate population structure -----####
# Plot ordinations colored by metadata
ord=capscale(IBS~1)
summary(ord)
ords=scores(ord,display="sites")
axes2plot=c(1:4) 
scores=data.frame(ord$CA$u[,axes2plot])

# "broken stick" model: expected random distribution of eigenvalues
plot(ord$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(ord$CA$eig/sum(ord$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(ord$CA$eig)),col="red") 

# Check significance of variables
adonis2(IBS ~ Bay, data=env) # p=0.007
adonis2(IBS ~ Location, data=env) # p=0.192
adonis2(IBS ~ Sex, data=env) # p=0.56
adonis2(IBS ~ Maturity, data=env) # p=0.309

# Plot colored by bay
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=env$Bay)) + 
  geom_point(aes(size=1, colour = env$Bay)) +
  theme_bw()+
  coord_equal()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  #geom_label(label=env$Sample)+
  guides(size = "none")

# Plot colored by sex
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=env$Sex)) + 
  geom_point(aes(size=1, colour = env$Sex)) +
  theme_bw()+
  coord_equal()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  geom_label(label=env$Sample)+
  guides(size = "none")

# Plot colored by maturity
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=env$Maturity)) + 
  geom_point(aes(size=1, colour = env$Maturity)) +
  theme_bw()+
  coord_equal()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  geom_label(label=env$Sample)+
  guides(size = "none")






###---- Admixture ----###
site<-env[,1:2]
names(site)<-c("INDIVIDUALS", "SITE")
bams<-read.table("bams.qc")
names(bams)<-"INDIVIDUALS"
bams$INDIVIDUALS<- gsub(".bam","",bams$INDIVIDUALS)
pop<-join(bams, site, by="INDIVIDUALS")
pop$SITE<- gsub("Area","",pop$SITE)

# K=2
q2<-read.table("Bonnet2_2.qopt") # name of the input file to plot, output of ngsAdmix run
names(pop)=c("ind","pop")
rownames(q2)=pop$ind
CSgt = colorRampPalette(colors = c("grey20", "grey40", "grey60", "grey80")) # create a simple colorscale
# Initial barplot
barplot(t(q2),
        col=CSgt(2),
        names=rownames(q2),
        las=2,
        space=0,
        border=NA,
        ylab="Admixture proportions for K=2")
# Barplot re-ordered by subpopulation
ord<-order(q2$V1, q2$V2) # re-ordered by subpopulation
barplot(t(q2)[,ord],
        col=CSgt(2),
        names=rownames(q2)[ord],
        las=2,
        space=0,
        border="black",
        ylab="Admixture proportions for K=2")
# Make color labels to go under barplot
cord=rownames(q2)[ord]
cord=substr(cord, 1, 1)
cordo = ifelse(cord == "B", "#F8766D", "#00BFC4")
plot(NULL, xlim=c(0,length(cordo)), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n")
rect(0:(length(cordo)-1), 0, 1:length(cordo), 1, col=cordo)
# Barplot grouped by site
TB=q2[70:89,]
TBord<-order(TB$V1, TB$V2) 
BB=q2[1:69,]
BBord<-order(BB$V1, BB$V2) 
par(mfrow = c(1,2))
barplot(t(TB)[,TBord],
        width=7,
        border= "black",
        col=CSgt(2),
        names=rownames(TB)[TBord],
        las=2,
        space=0,
        ylab="Admixture proportions for K=2")
barplot(t(BB)[,BBord],
        width=2,
        border= "black",
        col=CSgt(2),
        names=rownames(BB)[BBord],
        las=2,
        space=0)
# Plot PCoA colored by assignment to admixture subpopulations (k=2)
ggplot(scores,aes(scores[,1],scores[,2], asp=1, fill=q2$V1)) + 
  geom_point(aes(size=1, colour = q2$V1)) +
  theme_bw()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  geom_label(label=env$Sample)+
  guides(size = "none")

npops=ncol(q2)
cluster.admix=apply(q2[,1:npops],1,function(x) {return(paste(which(x>0.5),collapse=".")) })
# Get list of admixture groups
admix1 <- subset(rownames(scores), scores$cluster.admix == 1) #37 samples
admix1 <- paste0(admix1, ".bam")
#write.table(admix1, "ngsadmix_cluster1", sep="\t", col.names = F, row.names = F)
admix2 <- subset(rownames(scores), scores$cluster.admix == 2) #52 samples
admix2 <- paste0(admix2, ".bam")
#write.table(admix2, "ngsadmix_cluster2", sep="\t", col.names = F, row.names = F)


# K=3:
par(mfrow = c(2,1))
q3<-read.table("Bonnet2_3.qopt")
# Make a barplot (ordered by population)
ord<-order(q3$V1, q3$V2, q3$V3) #order it by population
barplot(t(q3)[,ord],
        col=CSgt(3),
        #names=pop$INDIVIDUALS[ord],
        las=2,
        space=0,
        border="black",
        ylab="Admixture proportions for K=3")
# K=4:
q4<-read.table("Bonnet2_4.qopt")
# Make a barplot (ordered by population)
ord<-order(q4$V1, q4$V2, q4$V3, q4$V4) #order it by population
barplot(t(q4)[,ord],
        col=CSgt(4),
        #names=pop$INDIVIDUALS[ord],
        las=2,
        space=0,
        border="black",
        ylab="Admixture proportions for K=4")




###----- Relatedness -----####
par(mfrow = c(1,1))
# reading long relatedness table, output of ngsRelate
rel=read.table("Bonnet.res",sep="\t",header=T)
# creating an empty square matrix
relm=matrix(0,nrow=length(unique(rel$a))+1,ncol=length(unique(rel$a))+1)
# filling up the square matrix with entries from "rab" column
for (a in unique(rel$a)) {
  for (b in unique(rel$b)) {
    if (b<=a) {next}
    relm[a+1,b+1]=relm[b+1,a+1]=rel[rel$a==a & rel$b==b,"rab"]
  }
}
diag(relm)=1
# add names to columns and rows 
bams=scan("bams.qc",what="character") # list of bam files
bams=sub(".bam","",bams)
dimnames(relm)=list(bams,bams)
# Plot dendrogram
hc_rel=(hclust(as.dist(1-relm),method="ave"))
plot(hc_rel)
# Remove outlier sample
goods=rownames(relm)
goods=goods[ !goods == 'BB18']
# subsetting all data for only the retained samples
relm=relm[goods,goods] 
hc_rel=(hclust(as.dist(1-relm),method="ave"))
plot(hc_rel)





###----- Genotype associations with morphology -----####
morph0<-env[,5:14]
morph0$Location=NULL
rownames(morph0)<-env$Sample
morph0<-na.omit(morph0)
goods<-rownames(morph0)
morph <-as.data.frame(morph0[goods,]) # subset only samples with sufficient morphological data
morph=morph0[,2:4]
Xscale <- as.data.frame(apply(morph,2,function(x){return(scale(rank(x)))})) # scale variables for RDA
IBS=IBS[goods,goods] # also subset genetic matrix for samples with morphological data

# Partial RDA conditional on sex and maturity
prda <-capscale(as.dist(IBS) ~. + Condition(morph0$Maturity)+Condition(morph0$Sex),  data=Xscale)
RsquareAdj(prda) 
anova(prda, perm=999) # Model not significant
# Forward variable selection
rda0<- capscale(as.dist(IBS) ~ 1 + Condition(morph0$Maturity)+Condition(morph0$Sex), Xscale) # null w/ just intercept
rdaG<- capscale(as.dist(IBS) ~ . + Condition(morph0$Maturity)+Condition(morph0$Sex), Xscale) # w/ all explanatory variables
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="forward") # Forward selection
Sel$anova # no variables retained

# Subset only mature sharks and rerun pRDA
mature0=env[grepl("mature", env$Maturity),]
mature0=env[!grepl("immature", env$Maturity),] # 67 samples retained
goods<-rownames(mature0)
mature <-as.data.frame(mature0[goods,])
mature$Sex=NULL
mature$Maturity=NULL
mature=mature[,c(1:3)]
Xscale <- as.data.frame(apply(mature,2,function(x){return(scale(rank(x)))})) # scale variables for RDA
IBS=as.matrix(read.table("Bonnet2.ibsMat")) # Need to re-subset IBS matrix
dimnames(IBS)=list(rownames(samples),rownames(samples))
IBS=IBS[goods,goods] 
prda <-capscale(as.dist(IBS) ~. + Condition(mature0$Sex),  data=Xscale)
RsquareAdj(prda) 
anova(prda, perm=999) # Model not significant

# Subset mature only, Plot colored by sex
ord2=capscale(IBS~1)
summary(ord2)
ords2=scores(ord2,display="sites")
axes2plot=c(1:4) 
scores2=data.frame(ord2$CA$u[,axes2plot])
ggplot(scores2,aes(scores2[,1],scores2[,2], asp=1, fill=mature0$Sex)) + 
  geom_point(aes(size=1, colour = mature0$Sex)) +
  theme_bw()+
  coord_equal()+
  xlab(names(scores2)[1])+
  ylab(names(scores2)[2])+
  #geom_label(label=mature0$Sample)+
  guides(size = "none")
adonis2(IBS ~ Sex, data=mature0) # p=0.798


# Subset only mature female sharks and rerun pRDA
fem0=mature0[grepl("F", mature0$Sex),] # 36 mature females retained
goods<-rownames(fem0)
fem <-as.data.frame(fem0[goods,])
fem=fem[,c(2:4)]
Xscale <- as.data.frame(apply(fem,2,function(x){return(scale(rank(x)))})) # scale variables for RDA
IBS=as.matrix(read.table("Bonnet2.ibsMat")) # Need to re-subset IBS matric
dimnames(IBS)=list(rownames(samples),rownames(samples))
IBS=IBS[goods,goods] 
prda <-capscale(as.dist(IBS) ~. ,  data=Xscale)
RsquareAdj(prda) 
anova(prda, perm=999) # Model not significant
plot(prda)

# Subset only mature male sharks and rerun pRDA
male0=mature0[!grepl("F", mature0$Sex),] # 31 mature males retained
goods<-rownames(male0)
male <-as.data.frame(male0[goods,])
male=male[,c(2:4)]
Xscale <- as.data.frame(apply(male,2,function(x){return(scale(rank(x)))})) # scale variables for RDA
IBS=as.matrix(read.table("Bonnet2.ibsMat")) # Need to re-subset IBS matric
dimnames(IBS)=list(rownames(samples),rownames(samples))
IBS=IBS[goods,goods] 
prda <-capscale(as.dist(IBS) ~. ,  data=Xscale)
RsquareAdj(prda) 
anova(prda, perm=999) # Model not significant



