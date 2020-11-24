mol delete all

cd ../../../Users/altin/Google\ Drive/MD_AURKA/2wtv/
mol load pdb complexx.pdb xtc md_0_1_noPBC.xtc
set mol0 2wtv-alisertib  

cd ../4j8m
mol load pdb complexx.pdb xtc md_0_1_noPBC.xtc
set mol1 4j8m-CD

cd ../../../rmsf-rmsd-salt-brigdes-2wtv-4j8m

set nf [molinfo 0 get numframes]	

#RMSF Alignment
for {set j 0} {$j<2} {incr j} {
set reference$j [atomselect $j "protein" frame 0]
# the frame being compared
set compare$j [atomselect $j "protein"]
set sel$j [atomselect $j "name CA"]
}

for {set frame 0} {$frame < $nf} {incr frame} {
		$compare0 frame $frame
		set trans_mat0 [measure fit $compare0 $reference0]
		$compare0 move $trans_mat0
		
		$compare1 frame $frame
		set trans_mat1 [measure fit $compare1 $reference1]
		$compare1 move $trans_mat1
		}

# Measure RMSF
set rmsf1w [measure rmsf $sel0 first 1 last 5000 step 1]
set rmsf1l [measure rmsf $sel0 first 4000 last 5000 step 1]

set rmsf2w [measure rmsf $sel1 first 1 last 5000 step  1]
set rmsf2l [measure rmsf $sel1 first 4000 last 5000 step 1]

set outfile [open rmsf_alisertib-CD_Whole-last.dat w]
for {set i 0} {$i < [$sel0 num]} {incr i} {
  puts $outfile "$mol0 [expr {$i+1}] whole [lindex $rmsf1w $i] last [lindex $rmsf1l $i]  $mol1 [expr {$i+1}] whole [lindex $rmsf2w $i] last [lindex $rmsf2l $i]"  
  }

close $outfile
	
# Alignment and RMSD	
set outfile [open Rmsd_AurA-$mol0.dat w]	

for {set frame 0} {$frame < $nf} {incr frame} {
		$compare0 frame $frame
		set trans_mat0 [measure fit $compare0 $reference0]
		$compare0 move $trans_mat0
		
		$compare1 frame $frame
		set trans_mat1 [measure fit $compare1 $reference1]
		$compare1 move $trans_mat1
		
		set rmsd0 [measure rmsd $compare0 $reference0]
		set rmsd1 [measure rmsd $compare1 $reference1]
        puts $outfile "RMSD $mol0 $frame $rmsd0 $mol1 $frame $rmsd1"
		}
close $outfile
						

# RGYR							
set outfile [open rgyr_2wtv-5g1x.dat w]

for {set j 0} {$j<2} {incr j} {
set sela$j [atomselect $j "protein and noh"]
# set selb$j [atomselect $j "serial 737 to 1178 1743 to 1875 and noh"]
# set selc$j [atomselect $j "serial 2926 to 3366 3932 to 4064 and noh"]
}

for {set i 0} {$i<$nf} {incr i} {
$sela0 frame $i
$sela0 update

$sela1 frame $i
$sela1 update

set rgyr0 [measure rgyr $sela0]
set rgyr1 [measure rgyr $sela1]

puts $outfile "$mol0 $i $rgyr0 $mol1 $i $rgyr1"
}  
close $outfile