library(dplyr)
library(ggplot2)
tmp <- read.csv("comma_deliminated_PositionScan_for_Aurora_First_frame_whole_mutations_for_heatmap.txt", sep=" ", header=T, row.names = 1)

#replace rows if necessary
#new data to be replaced
# df_rep <-read.csv("positionscan2.txt", sep="\t", header=FALSE)

# df[match(df_rep$V1, df$V1), ] <- df_rep
# df

#transpose and reformat

#tmpa <- data.frame(
  #X=tmp$V2,
  #ind=rep(1:21, nrow(tmp)/21))

df2 <-read.csv("aa_prop.txt", sep="\t")
rownames(df2)<- df2$X
df2[1]<-NULL
df2=df2[,order(colnames(df2))]

#combine data with amino acid properties

s<- rbind (df2,tmp)

#order according to mass
f=s[order(s[1,])]
z1=f[-c(1:nrow(df2)),]

#write.table(z1, file="mass_ordered_apo_np_298_0ns_PS.txt")


#order according to occurence
#f=s[order(s[2,])]
#z2=f[-c(1:nrow(tmp2)),]
#write.table(z2, file="occurence_ordered_apo_np_298_0ns_PS.txt")


#order according to hyropathy
#f=s[order(s[9,])]
#z3=f[-c(1:nrow(tmp2)),]
#write.table(z3, file="hydropathy_ordered_apo_np_298_0ns_PS.txt")


#order according to alpha
#f=s[order(s[10,])]
#z4=f[-c(1:nrow(tmp2)),]
#write.table(z4, file="alphaprop_ordered_apo_np_298_0ns_PS.txt")


#order according to beta
#f=s[order(s[12,])]
#z5=f[-c(1:nrow(tmp2)),]
#write.table(z5, file="betaprop_ordered_apo_np_298_0ns_PS.txt")


#order according to coil
#f=s[order(s[13,])]
#z6=f[-c(1:nrow(tmp2)),]
#write.table(z6, file="coilprop_ordered_apo_np_298_0ns_PS.txt")

#ss <- rbind (z1, z2, z3, z4, z5, z6)
#write.table(ss, file="ordered_MOHABC_apo_np_298_0ns_PS.txt")


tmp2 <- reshape2::melt(as.matrix(z1))
colnames(tmp2) <- c("pos", "aminoacid", "count")

mat1factor <- tmp2 %>%
  # convert state to factor and reverse order of levels
  mutate(state=factor(aminoacid,levels=rev(sort(unique(aminoacid)))))  %>%
  
  # create a new variable from count
  mutate(countfactor=cut(count,breaks=c("-100","-5","-0.5","0.5","5","10","100",max(count,na.rm=T)),
                         labels=c("<-10","(-10,-5)","(-5,-0.5)","(-0.5,0.5)","(0.5,5)","(5,10)",">10")))  %>%
  # change level order
  mutate(countfactor=factor(as.character(countfactor),levels=rev(levels(countfactor))))

textcol <- "grey40"

# further modified ggplot
p <- ggplot(mat1factor,aes(x = pos ,y = aminoacid,fill = countfactor))
p <- p + geom_tile(colour = "white", size = 0.25)
p <- p + guides(fill=guide_legend(title="????G kcal/mol"))
p <- p + labs(x="",y="",title="PositionScan Aurora First Frame")
p <- p + scale_y_discrete(expand=c(0,0))
p <- p + scale_x_discrete(expand=c(0,0), breaks=c("R126" ,"E134" ,"G136" ,"R137" ,"L139" ,"G145" ,"E152" ,"S155" ,"G173" ,"V174" ,"E175" ,"R180" ,"R189" ,"P191" ,"L194" ,"H201" ,"D202" ,"R205" ,"V206" ,"G216" ,"V218" ,"R220" ,"L222" ,"K227" ,"L240" ,"S245" ,"Y246" ,"K250" ,"R251" ,"H254" ,"K258" ,"E260" ,"A267" ,"D274" ,"W277" ,"V279" ,"A281" ,"R285" ,"T287" ,"A290" ,"G291" ,"D294" ,"P297" ,"I301" ,"R304" ,"M305" ,"D311" ,"L312" ,"G316" ,"L318" ,"L323" ,"G325" ,"P328" ,"T333" ,"I341" ,"R343" ,"G355" ,"R357" ,"D358" ,"L359" ,"S361" ,"P368" ,"L378" ,"T384" ,"A385" ,"S387" ,"S388"))
#breaks for nmyc -> breaks=c("L61" ,"S62" ,"P63" ,"S64" ,"R65" ,"G66" ,"F67" ,"A68" ,"E69" ,"H70" ,"S71" ,"S72" ,"E73" ,"P74" ,"P75" ,"S76" ,"W77" ,"V78" ,"T79" ,"E80" ,"M81" ,"L82" ,"L83" ,"E84" ,"N85" ,"E86" ,"L87" ,"W88" ,"G89")
p <- p + scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4","#ddf1da"), na.value = "grey90")
p <- p + theme_grey(base_size=10)
p <- p + theme(legend.position="right",legend.direction="vertical",
               legend.title=element_text(colour=textcol, size = "14"),
               legend.margin=margin(grid::unit(0,"cm")),
               legend.text=element_text(colour=textcol,size=14,face="bold"),
               legend.key.height=grid::unit(3,"cm"),
               legend.key.width=grid::unit(1,"cm"),
               axis.text.x=element_text(size=11,colour=textcol, angle = 60, hjust = 1),
               axis.text.y=element_text(vjust=0.2,colour=textcol, size = 11),
               axis.ticks=element_line(size=0.4),
               plot.background=element_blank(),
               panel.border=element_blank(),
               plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
               plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"))

#the figure downloaded by this code but it does not look well, so I take a snapshot after zooming output "p".

ggsave(p,filename="Aurora_First_frame_whole_mutations_size.png", height=5.5,width=5.5,units="in",dpi=200)

p




