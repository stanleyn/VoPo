#purpose: generate a distribution of AUC
source('VoPo_main/runClassif.R')
library('ggplot2')
library('reshape2')

###please update!!!!
Build_Preg=readRDS('Processed/Build_Pregnancy.rds')

#Build_Preg=readRDS('~/Clean_BClust/Pregnancy/B2_March312020')
Meta_Preg=readRDS('Processed/Meta_Pregnancy.rds')

FuncDF=Build_Preg$FProp
IterNumClus=Build_Preg$IterNumClus
#rename to not perturb RF
ForC=1:ncol(FuncDF)
NewName=paste(colnames(FuncDF),ForC,sep='_')
colnames(FuncDF)=NewName
Y=Meta_Preg$Class

#do baseline distribution first
Base=c()
for(i in 1:100){
#choose an iteration to get features from and get relavent columns
rndIter=sample(1:length(IterNumClus),1)
NumClus=IterNumClus[1]
start=((rndIter-1)*NumClus)+1
end=NumClus*rndIter
SubDF=FuncDF[,start:end]
AUC=runClassif(SubDF,Meta_Preg$Class,IterNumClus[1],IterNumClus[1],0.5,1,as.character(Meta_Preg$Subject),35)
Base=c(Base,AUC)
}

#get repeated distribution
Boot=runClassif(FuncDF,Meta_Preg$Class,40,IterNumClus,0.5,100,as.character(Meta_Preg$Subject),35)


#plotting
AllBL=rbind(Base,Boot)
rownames(AllBL)=c('Baseline','BootStrap')
DF=melt(AllBL)
names(DF)=c('Method','perm','AUC')

ValVec=c('black','#FF0266')
p26 <- ggplot(DF, aes(x=Method, y=AUC,color=Method)) +
	  geom_boxplot(notch=FALSE,lwd=1)+theme(text = element_text(size=15))+xlab('')+ylab('')+scale_color_manual(values = ValVec)+ggtitle('')
p26=p26+geom_point(size=.7,position=position_jitterdodge(),alpha=.5)+theme_classic()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size=14))      
p26=p26+theme(legend.position='bottom')+ggtitle('')+ylab('AUC')
ggsave('OutDir/Preg_Dist.pdf',p26,width=4,height=4)

