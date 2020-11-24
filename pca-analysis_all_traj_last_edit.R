library(bio3d)
library(labdsv)
dcd <- read.dcd("unwrapped_5g1x_complex_all_namd3.dcd")

pdb <- read.pdb("unwrapped_5g1x_complex_all_namd3.pdb")

#scree plot values
df2 <- data.frame(prop.table(pc2$L))
write.csv(df2[1:20,],file = "unwrapped_dissociated_complex_namd3_pc-L.csv")

#rescont
data2 <- pc2$au
write.csv(data2,"unwrapped_dissociated_complex_namd3_pc-au.csv")

#pc
data_z2 <- data.frame(pc2$z)
write.csv(data_z2[,1:5],"unwraped_5g1x_complex_namd3_pc-z.csv")


ca.inds <- atom.select(pdb, elety="CA")

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)

dim(xyz) == dim(dcd)
#RMSD
rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
write.csv(rd,"unwrapped_dissociated_only_nmyc_RMSD.csv")

plot(rd, typ="l", ylab="RMSD", xlab="Frame No.", main="5g1x no-nmyc")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)
#quick histogram
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram of 5g1x no-nmyc", xlab="RMSD")
lines(density.default(rd), col="blue", lwd=2)

summary(rd)

#RMSF
rf <- rmsf(xyz[,ca.inds$xyz])
write.csv(rf,"unwrapped_dissociated_only_nmyc_whole_RMSF.csv")
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l", lwd=2, main="5g1x_no-nmyc")

plot
  #PCA (PRINCIPLE COMPONENT ANALYSIS) (i.e. trajectory frames) colored from blue to red in order of time
pc2 <- pca.xyz(xyz[,ca.inds$xyz])


plot(pc, col=bwr.colors(nrow(xyz))) #, xlim=c(-30,40), ylim=c(-20,30))
#for representation of each pc'es separately pc.axes=1:2
pca12 <- plot(pc, col=bwr.colors(nrow(xyz)), pc.axes=1:2, xlim=c(-40,35), ylim=c(-20,30))  #xlim = c(-30,20),  ylim = c(-30,20) )
pca32 <- plot(pc, col=bwr.colors(nrow(xyz)), pc.axes=3:2,  xlim=c(-40,35), ylim=c(-20,30))
pca13 <- plot(pc, col=bwr.colors(nrow(xyz)), pc.axes=c(1,3),  xlim=c(-40,35), ylim=c(-20,30))


plot(pc$U)

11111#HIERHARCAL CLUSTERING
hc <- hclust(dist(pc$z[,1:3]))
grps <- cutree(hc, k=4)
plot(pc, col=grps) #, xlim=c(-40,35), ylim=c(-40,35))
pc_axes12 <- plot(pc, col=grps,pc.axes=1:2, xlim=c(-40,35), ylim=c(-20,30))
pc_axes32 <- plot(pc, col=grps,pc.axes=3:2,  xlim=c(-40,35), ylim=c(-20,30)) 
pc_axes13 <- plot(pc, col=grps,pc.axes=c(1,3), xlim=c(-40,35), ylim=c(-20,30))

#plot PC all  
par(mfrow = c(4, 1), cex = 0.75, mar = c(3, 4, 1, 1) )
plot.bio3d(pc1$au[,1], ylab="PC_all (A)", xlab="Residue Position", typ="l", col="red",lwd=2) + points(pc$au[,2], typ="l", col="blue",lwd=2) + points(pc$au[,3], typ="l", col="green",lwd=2.5)
#plot PC1 only
plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l", col="red",lwd=2)
#plot PC2 only
plot.bio3d(pc$au[,2], ylab="PC2 (A)", xlab="Residue Position", typ="l", col="blue", lwd=2)
#plot PC3 only
plot.bio3d(pc$au[,3], ylab="PC3 (A)", xlab="Residue Position", typ="l", col="green", lwd=2)

#plot both PC1 and PC2
p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file="pc1_5g1x_no-nmyc.pdb")
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file="pc2_5g1x_no-nmyc.pdb")
p3 <- mktrj.pca(pc, pc=3, b=pc$au[,3], file="pc3_5g1x_no-nmyc.pdb")


#write trajectory's as AMBER NetCDF format from low to high atomic displacement
#convert first cdf  then mdcrd to open in vmd in windows version
library("ncdf4")
write.ncdf(p1, "trj_pc1_5g1x_no-nmyc_amber.nc")
write.ncdf(p2, "trj_pc2_5g1x_no-nmyc_amber.nc")
write.ncdf(p3, "trj_pc3_5g1x_no-nmyc,_amber.nc")
#cross-correlation analysis
cij<-dccm(xyz[,ca.inds$xyz])
plot(cij, main="Residue Cross Correlation of 5g1x_no-nmyc")

# View the correlations in pymol (pymol u sildim)
pymol.dccm(cij, pdb, type="launch")

write.ncdf(p1, "trj_pc1.nc")
write.ncdf(p2, "trj_pc2.nc")
write.ncdf(p3, "trj_pc3.nc")

#########################################