library(bio3d)
library(dplyr)

entry <- as.vector(read.table(file = "pdb_entry_aurora.txt", header = T))
all_strct <- get.pdb(entry$pdb_entry)
split_chains <- pdbsplit(all_strct, ids = "_A")
pdbs_aurora <- pdbaln(split_chains, fit = T,  web.args = list(email = "altinerpinar@gmail.com"))
# Calculate sequence identity
pdbs_aurora$id <- substr(basename(pdbs$id),1,6)
fasta2 <- read.fasta("aurora-aln.fa", rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)
pim <-seqidentity(fasta2)
xyz2 <- pdbfit(pdbs_aurora)
rmsd(xyz2)

ind2 <- grep("5g1x", pdbs_aurora$id) 

core <- core.find(pdbs_aurora)
col=rep("black", length(core$volume))
col[core$volume<2]="pink"; col[core$volume<1]="red"
plot(core, col=col)
core.inds <- print(core, vol=1.0)
#write.pdb(xyz=pdbs$xyz[1,core.inds$xyz], file="quick_core_sod1_2v0a.pdb")

#xyz_new2 <- pdbfit(pdbs_aurora, core.inds )

rd <- rmsd(xyz)
rd

hist(rd, xlab="RMSD (�)", main="Histogram of RMSD", breaks = 72)

# RMSD clustering
hc.rd <- hclust(as.dist(rd))

pdbs_aurora$id <- substr(basename(pdbs_aurora$id), 1, 6)
hclustplot(hc.rd, labels=pdbs_aurora$id, k=3 ,cex=0.5,ylab="RMSD (�)", main="RMSD Cluster Dendrogram", fillbox=FALSE)

#blue Aurora B, yellow Aurora C, green like 5g1x, red inactive loop, black missing loop , gray missing loop Aurora B
# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs_aurora$ali)
gaps.pos <- gap.inspect(pdbs_aurora$xyz)


# Tailor the PDB structure to exclude gap positions for SSE annotation
id = grep("5g1x", pdbs_aurora$id)
ref.pdb = trim.pdb(pdb, inds=atom.select(pdb, resno = pdbs_aurora$resno[id]))



# Plot RMSF with SSE annotation and labeled with residue numbers (Figure 8.) residue
#rf <- rmsf(pdbs$xyz)
rf2 <- rmsf(xyz[, gaps.pos$f.inds])
plot.bio3d(rf2, ylab="RMSF (�)", xlab="Residue No.", typ="l")
ref_pdb <- read.pdb("5g1x")
tor <- torsion.pdb(ref_pdb)

# Basic Ramachandran plot (Figure 9)
plot(tor$phi, tor$psi, xlab="phi", ylab="psi")


############
# Locate the two structures in pdbs
# For comparison the best BLAST result pdb and 2wtv are chosen.
ind.a <- grep("5G1X_A", pdbs_aurora$id)
ind.b <- grep("2WTV_A", pdbs_aurora$id)

# Exclude gaps in the two structures to make them comparable
gaps.xyz2 <- gap.inspect(pdbs_aurora$xyz[c(ind.a, ind.b), ])
a.xyz <- pdbs_aurora$xyz[ind.a, gaps.xyz2$f.inds]
b.xyz <- pdbs_aurora$xyz[ind.b, gaps.xyz2$f.inds]

# Compare CA based pseudo-torsion angles between the two structures
a <- torsion.xyz(a.xyz, atm.inc=1)
b <- torsion.xyz(b.xyz, atm.inc=1)
d.ab <- wrap.tor(a-b)
d.ab[is.na(d.ab)] <- 0

# Plot results with SSE annotation
plot.bio3d(abs(d.ab), resno=pdb, sse=pdb, typ="h", xlab="Residue No.", ylab = "Difference Angle")

a <- dm.xyz(a.xyz)
b <- dm.xyz(b.xyz)

plot.dmat( (a - b), nlevels=10, grid.col="gray", xlab="5AAD", ylab="2WTV")


# Do PCA
grps <- cutree(hc.rd, k=3) #additional hier. clust based on rmsd
pc.xray <- pca.xyz(xyz2[, gaps.pos$f.inds])

pc.xray
        
plot(pc.xray, col = grps)
#coloring method u active ve inactive olarak d�zelt.
plot(pc.xray$z[,1:2], col=grps, type = "p", lwd=5, xlab = "PC1", ylab = "PC2")

text(pc.xray$z[,1], pc.xray$z[,2], pos= 1,label=pdbs_aurora$id <- substr(basename(pdbs_aurora$id), 1, 6))


# Left-click on a point to label and right-click to end
identify(pc.xray$z[,1:2], labels=basename.pdb(pdbs$id))

par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
plot.bio3d(pc.xray$au[,1], resno=ref.pdb, ylab="PC1")
plot.bio3d(pc.xray$au[,2], resno=ref.pdb, sse=ref.pdb, ylab="PC2")
plot.bio3d(pc.xray$au[,3], resno=ref.pdb, sse=ref.pdb, ylab="PC3")
# See pc1.pdb
mktrj.pca(pc.xray, pc=3, file="true-pc3.pdb")

#Compared only aurora A
hc <- hclust(dist(pc.xray$z[,1:2]))
new_label <- substr(basename(hc$labels), 1, 6)
grps <- cutree(hc, h=30)
cols <- c("black", "red", "green")[grps]
plot(pc.xray, pc.axes=1:2, col=cols)



#Dendrogram plot
names(cols) <- pdbs_aurora$id
hclustplot(hc, colors=cols, ylab="Distance in PC Space", label=new_label ,main="PC1-2", cex=0.5, fillbox=FALSE)

#for closer point identification
identify(pc.xray$z[,1], pc.xray$z[,2], labels=pdbs$id)

#rescont
data <- pc.xray$au
write.csv(data,"true-aurora-allpc-au.csv")

#pc
data_z <- data.frame(pc.xray$z)
write.csv(data_z[,1:5],"true-aurora_pc-z_col1-5.csv")

#???scree plot
scree <- pc.xray$U
scree1 <- prop.table(scree)
write.csv(scree1, "true-scree_plot_aurora.csv")
