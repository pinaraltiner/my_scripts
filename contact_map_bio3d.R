library(bio3d)
library(ggplot2)
library(plot.matrix)
#read pdb 
pdb <- read.pdb("5g1x_rep1_first_frame.pdb")
## Atom Selection indices
inds <- atom.select(pdb, "calpha")
## Reference contact map
ref.cont <- cmap( pdb$xyz[inds$xyz], dcut=6, scut=3 )
plot.cmap(ref.cont)
##- Read Traj file
trj <- read.dcd( "new_5g1x.dcd") 
ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=trj,
              fixed.inds=ca.inds$xyz,
              mobile.inds=ca.inds$xyz)
dim(xyz) == dim(trj)

#contact map for one frame
pdb <- read.pdb("5g1x_rep1_first_frame.pdb")
pdb2 <- read.pdb("5g1x_rep1_last_frame.pdb")
#set indices
inds <- atom.select(aurora, "calpha")
inds2 <- atom.select(pdb2, "calpha")
#pdb$xyz[inds$xyz] inds lerin belirttiği atomların pdb deki xyz si olmalı
first <- cmap(pdb$xyz[inds$xyz], dcut=4, scut=3)
last <- cmap( pdb2$xyz[inds2$xyz], dcut=4, scut=3)
#write.csv(cm_first, file="first_frame_cmp_mat_5g1x.csv")



first1 <- cmap(pdb$xyz[mycn$xyz] ,dcut=6, scut=3)
last1 <- cmap( pdb2$xyz[mycn2$xyz], dcut=6, scut=3 )
plot.cmap(first1, col="red", sse= dssp(pdb, exefile = "/usr/bin/dssp"))
plot.cmap(t(last1), col="green", sse= dssp(pdb, exefile = "/usr/bin/dssp"),add = T)




plot.cmap(mycn, col="blue", sse= dssp(pdb, exefile = "/usr/bin/dssp"), add = FALSE)


dssp(pdb, exefile = "/usr/bin/dssp")


## For each frame of trajectory
sum.cont <- NULL
for(i in 1:nrow(trj)) {
  
## Contact map for frame 'i'
  cont <- cmap(trj[i,inds$xyz], dcut=6, scut=3)

  ## Product with reference
  prod.cont <- ref.cont * cont
  sum.cont <- c(sum.cont, sum(prod.cont,na.rm=TRUE))
}

plot(sum.cont, typ="l")

## 
# DM
d <- dm(as.array(pdb, inds=atom.select()))
write.csv(d, file="dist_mat_5g1x.csv")
# Plot DM
##filled.contour(d, nlevels = 4)
##plot(d)
plot(d,resnum.1 = pdb$atom[pdb$calpha,"resno"],
     color.palette = mono.colors, xlab="Residue Number", 
     ylab="Residue Number")

rx <- range(x<-pdb$atom[pdb$calpha,"resno1"])
ry <- range(y<-pdb$atom[pdb$calpha,"resno"])
##-- Residue-wise distance matrix based on the minimal distance between all available atoms
inds <- atom.select(pdb, "calpha")
inds2 <- atom.select(pdb2, "calpha")
a <- dm.xyz(pdb$xyz[inds$xyz], inds=pdb2$atom[inds2$atom,"resno"])
colnames(a) <- c("126","127","128","129","130","131","132","133","134","135","136","137","138","139",
                 "140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155",
                 "156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171",
                 "172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187",
                 "188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203",
                 "204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219",
                 "220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235",
                 "236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251",
                 "252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267",
                 "268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283",
                 "284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299",
                 "300","301","302","303","304","305","306","307","308","309","310","311","312","313","314","315",
                 "316","317","318","319","320","321","322","323","324","325","326","327","328","329","330","331",
                 "332","333","334","335","336","337","338","339","340","341","342","343","344","345","346","347",
                 "348","349","350","351","352","353","354","355","356","357","358","359","360","361","362","363",
                 "364","365","366","367","368","369","370","371","372","373","374","375","376","377","378","379",
                 "380","381","382","383","384","385","386","387","388","389","61","62","63","64","65","66","67",
                 "68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86",
                 "87","88","89")
rownames(a) <- c("126","127","128","129","130","131","132","133","134","135","136","137","138","139",
                 "140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155",
                 "156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171",
                 "172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187",
                 "188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203",
                 "204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219",
                 "220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235",
                 "236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251",
                 "252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267",
                 "268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283",
                 "284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299",
                 "300","301","302","303","304","305","306","307","308","309","310","311","312","313","314","315",
                 "316","317","318","319","320","321","322","323","324","325","326","327","328","329","330","331",
                 "332","333","334","335","336","337","338","339","340","341","342","343","344","345","346","347",
                 "348","349","350","351","352","353","354","355","356","357","358","359","360","361","362","363",
                 "364","365","366","367","368","369","370","371","372","373","374","375","376","377","378","379",
                 "380","381","382","383","384","385","386","387","388","389","61","62","63","64","65","66","67",
                 "68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86",
                 "87","88","89")
plot.dmat(a, resnum.1 = row.names(a), resnum.2= row.names(a))
b <- dm.xyz(pdb2$xyz[inds2$xyz], grpby=pdb2$atom[inds2$atom,"resno"])
a1 <- cmap(a)
pl
plot.dmat(a,grid.nx = seq(1,ncol(a),2) , grid.ny =  seq(1,ncol(a),2))
range <- pdb2$atom[inds2$atom,"resno"]
plot((a - b),grid.col=TRUE, xlab="first", ylab="last",grid = TRUE, resnum.1 = row.names(a), resnum.2 = row.names(a))
plot.bio3d(a-b)
l <- dm.xyz(pdb2$xyz[inds2$xyz], scut=3, dcut=4)
plot(l)
dim <- ncol(a)
image(pdb2$atom[inds2$atom,"resno"], pdb2$atom[inds2$atom,"resno"], a, axes = FALSE, xlab="", ylab="")
axis(1, pdb2$atom[inds2$atom,"resno"], cex.axis = 0.5, las=3)
axis(2, pdb2$atom[inds2$atom,"resno"], cex.axis = 0.5, las=1)
# Load the library
library("lattice")

# Dummy data
data <- matrix(runif(100, 0, 5) , 10 , 10)
colnames(a) <- letters[c(1:10)]
rownames(data) <- paste( rep("row",10) , c(1:10) , sep=" ")

# plot it flipping the axis
levelplot(as.matrix(a), col.regions=heat.colors(100))


