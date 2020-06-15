#purpose: generate a distribution of AUC
library('ggplot2')
library('reshape2')

source('VoPo_main/runClassif.R')
Build_Stroke=readRDS('Processed/Build_Stroke.rds')
Meta_Stroke=readRDS('Processed/Meta_Stroke.rds')

FuncDF=Build_Stroke$FProp
IterNumClus=Build_Stroke$IterNumClus
#rename to not perturb RF
ForC=1:ncol(FuncDF)
NewName=paste(colnames(FuncDF),ForC,sep='_')
colnames(FuncDF)=NewName
Y=factor(Meta_Stroke$Class)

#do baseline distribution first
Base=c()
for(i in 1:100){
#choose an iteration to get features from and get relavent columns
rndIter=sample(1:length(IterNumClus),1)
NumClus=IterNumClus[1]
start=((rndIter-1)*NumClus)+1
end=NumClus*rndIter
SubDF=FuncDF[,start:end]
AUC=runClassif(SubDF,Y,IterNumClus[1],IterNumClus[1],0.5,1,as.character(Meta_Stroke$Subject),35)
Base=c(Base,AUC)
}

#get repeated distribution
Boot=runClassif(FuncDF,Y,40,IterNumClus,0.5,100,as.character(Meta_Stroke$Subject),35)

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
ggsave('OutDir/Stroke_Dist.pdf',p26,width=2,height=4)

