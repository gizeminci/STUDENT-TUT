# Test case to control the tcl functions of center-of-mass and radius-of-gyration calculation

# What you should do is;
# 1. setup missing parameters
# 2. printout agglomeration list, center of mass of the agglomerate, radius of gyration of the agglomerate
# 3. print out particle positions 
# 4. calculate center of mass and Rg from the positions of the particles
# 5. compare the results


#----------------------------------------------------------------------------------------
# HALF OF THE BOX LENGTH MUST BE LARGER THAN THE PARTICLE NUMBER for small test cases
require_feature "COLLISION_DETECTION"
require_feature "VIRTUAL_SITES_RELATIVE"
cellsystem domain_decomposition -no_verlet_list

# Setting some parameters 
#----------------------------------------------------------------------------------------
# number of particles in the system (lets only use 6 particles)


# simulation box parameters (as you wish)


# time step and skin


# thermostat (use Langevin)


# lennard-jones


# bonds (use angular potential and harmonic bond)


# vmd setup (you do not need to do anything here)
set pidlist [list]
set fl [open line-test.vtf a]

# Place Particles (define the particle positions as they create a line shhape)



# vmd setup (you do not need to do anything here)
writevsf $fl

# Periodic boundary condition


# Collision setup


# Particle control (you do not need to do anything here)
set min 0
set cap 20
while { $min < 0.85 } {
    # set ljforcecap
    inter ljforcecap $cap
    # integrate a number of steps, e.g. 20
    integrate 40
    # check the status of the sytem
    set min [analyze mindist 0 0]
    # this is a shortcut for 'set cap [expr $cap+10]'
    incr cap 20
}
#puts "Warmup finished. Minimal distance now $min"
# turn off the ljforcecap, which is done by setting the 
# force cap value to zero:
inter ljforcecap 0
#---------------------------------------------------------#

# vmd setup
for {set i 0} {$i <$particle} {incr i} {
    lappend pidlist $i
    vtfpid $i
}

# TO FIND FRACTAL DIMENSION

set initial_radius $sigma
set increase [expr 0.75*$sigma]
set ListForDf(0) [list]

# compare list a with b to find same elements
 #----------------
 proc slistcomp {a b} {
  set same {}
  foreach i $a {
    if {[lsearch $b $i]>-1} {
      lappend same $i
    }
  }
  return $same
 }
 #----------------
 # compare list a with b to find different elements
 #----------------
 proc dlistcomp {a b} {
  set diff {}
  foreach i $a {
    if {[lsearch -exact $b $i]==-1} {
      lappend diff $i
    }
  }
  return $diff
 }
#----------------------------------------------------------

proc min_img_bond_length { pos1 pos2 box } {
    set dim [llength $pos1]
    set res 0.0
    for {set j 0} { $j < $dim } {incr j} {
	set dist [expr abs([lindex $pos1 $j]-[lindex $pos2 $j])]
	while { $dist > [expr $box/2.0] } { 
	    set dist [expr $dist - $box] 
	}
	set res [expr $res+ [sqr $dist] ]
    }
    return [expr sqrt($res)]
}
#----------------------------------------------------------
proc com_for_new_positions {coordinatenum box aggregatelist} {
	set c [llength $aggregatelist]
	set procpart0 [lindex $aggregatelist 0]
	set procpart0_coor [lindex [part $procpart0 print folded_position] $coordinatenum]
	set total_coordinate $procpart0_coor
	if {$procpart0_coor < [expr $box/2.0] || $procpart0_coor == [expr $box/2.0]} {
	for {set procpart1 1} {$procpart1<$c} {incr procpart1} {
		if {[lindex [part [lindex $aggregatelist $procpart1] print folded_position] $coordinatenum] < [expr $box/2.0] || [lindex [part [lindex $aggregatelist $procpart1] print folded_position] $coordinatenum] == [expr $box/2.0]} {
			set total_coordinate [expr $total_coordinate + [lindex [part [lindex $aggregatelist $procpart1] print folded_position] $coordinatenum]]
		} else {
			set new_coordinate [expr [lindex [part [lindex $aggregatelist $procpart1] print folded_position] $coordinatenum] - $box ]
			set total_coordinate [expr $total_coordinate + $new_coordinate]
		}
	}
	set com [expr $total_coordinate/$c]
	} else {
	for {set procpart1 1} {$procpart1<$c} {incr procpart1} {
		if {[lindex [part [lindex $aggregatelist $procpart1] print folded_position] $coordinatenum] > [expr $box/2.0]} {
			set total_coordinate [expr $total_coordinate + [lindex [part [lindex $aggregatelist $procpart1] print folded_position] $coordinatenum]]
		} else {
			set new_coordinate [expr [lindex [part [lindex $aggregatelist $procpart1] print folded_position] $coordinatenum] + $box ]
			set total_coordinate [expr $total_coordinate + $new_coordinate]
		}
	}
	set com [expr $total_coordinate/$c]		
	}
return $com
}

proc com_for_own_positions {coordinatenum box aggregatelist} {
	set c [llength $aggregatelist]
	set total_coordinate 0
	for {set procpart 0} {$procpart<$c} {incr procpart} {
		set coordinate [lindex [part [lindex $aggregatelist $procpart] print folded_position] $coordinatenum]
		set total_coordinate [expr $total_coordinate + $coordinate]
	}
	set com [expr $total_coordinate/$c]
	return $com
}

#---------------INTEGRATION---------------------------#

set n_cycle 3000
set n_steps 100

set i 0 
set realnumofagglists 0
set totalcollidedparticles 0

for {set q 0} {$q<$particle} {incr q} {
set BondList($q) [list]
set RCheckList($q) [list]
}
set AggList(0) [list]
set ListForDf(0) [list]
set numofagglist 1
set na 1
#------------------------------------------------------------------------------------------------------------
while { $i<$n_cycle } {   
     integrate $n_steps

set realnumofagglists 0
set totalcollidedparticles 0

for {set q 0} {$q<$particle} {incr q} {
set BondList($q) [list]
}

# Now Create BondLists

for {set n 0} {$n<$particle} {incr n} {

	set initialbond [part $n print bonds]
	set bondlist [lindex $initialbond 0]

	foreach bond $bondlist {
		set qq [lindex $bond 1]
		lappend BondList($n) $qq
		lappend BondList($qq) $n
	}
}

# Now Create CheckList including the target particle and all bonded particles

set CheckList [list]
set AggList(0) [list]
for {set a 0} {$a<$particle} {incr a} {

	set BondsOfParticle0 [lindex [list $BondList($a)] 0]	
	set numofbond1 [llength $BondsOfParticle0]
if {$numofbond1>0} {
	set CheckList [lappend CheckList $a]
	set criter [llength $CheckList]
	set b 0
	while {$b<$criter} {
		set check [lindex $CheckList $b]
		set BondsOfParticle [lindex [list $BondList($check)] 0]	
		set numofbond [llength $BondsOfParticle]
		if {$numofbond>0} {
			for {set c 0} {$c<$numofbond} {incr c} {
				set cc [lindex $BondsOfParticle $c]
				set kc [lsearch $CheckList $cc]
				if {$kc<0} {
					set CheckList [lappend CheckList $cc]
				}
			}			
			set criter [llength $CheckList]
			set b [expr $b+1]
		} else {
			set criter [llength $CheckList]
			set b [expr $b+1]
		}
	}
set RCheckList($a) [list $CheckList]
}
	set CheckList [list]
}
for {set p 0} {$p<$particle} {incr p} {
set aakk 0
set numofchecklist [llength $RCheckList($p)]
if {$numofchecklist>0} {
	for {set d 0} {$d<$numofagglist} {incr d} {
		set ak [slistcomp [lindex $RCheckList($p) 0] [lindex $AggList($d) 0]]
set nak [llength $ak]
if {$nak<[llength [lindex $RCheckList($p) 0]] && $nak>0} {
set notinlist [dlistcomp [lindex $RCheckList($p) 0] [lindex $AggList($d) 0]]
set updatedlist0 [lindex $AggList($d) 0]
set AggList($d) [list [concat $updatedlist0 $notinlist]]
#puts "THE AGGLOMERATION LIST($d) IS UPDATED--->$AggList($d)"
}
if {$nak<1} {
set akk 0
} else {
set akk 1
}
set aakk [expr $aakk || $akk]
}
if {$aakk==0} {
set AggList($na) [lappend AggList($na) [lindex $RCheckList($p) 0]]
#puts "NEW AGGLOMERATION LIST($na) IS CREATED--->$AggList($na)"
set na [expr $na+1]
set numofagglist [expr $numofagglist+1]
} else {
}

}

}
	for {set al 1} {$al<$numofagglist} {incr al} {
		for {set all [expr $al+1]} {$all<$numofagglist} {incr all} {
			set compareAggLists [slistcomp [lindex $AggList($al) 0] [lindex $AggList($all) 0]]
			set numofcomp [llength $compareAggLists]
			if {$numofcomp>0} {
				set notinagglist [dlistcomp [lindex $AggList($all) 0] [lindex $AggList($al) 0]]
				set updatedagglist [lindex $AggList($al) 0]
				set AggList($al) [list [concat $updatedagglist $notinagglist]]
				set AggList($all) [list]
			}
		}
		
	}

for {set kp 1} {$kp<$numofagglist} {incr kp} {
	set agglistlength [llength [lindex $AggList($kp) 0]]
	if {$agglistlength>0} {
		set realnumofagglists [expr $realnumofagglists+1]
		set totalcollidedparticles [expr $totalcollidedparticles+$agglistlength]
#puts "Agglist($kp) is $AggList($kp)"
	}
}

#*---------------------------------------------------------------*

for {set kp 1} {$kp<$numofagglist} {incr kp} {

set distancelist2 [list]
set alist [lindex $AggList($kp) 0]
set control [llength $alist]

if {$control>0} {

set ab 0

while {$ab<$control} {

set particle1 [lindex $alist $ab]
set b [expr $ab+1]

while {$b<$control} {

	set particle2 [lindex $alist $b]

set dist2 [min_img_bond_length [part $particle1 print folded_position] [part $particle2 print folded_position]	$box_length]

lappend distancelist2 $dist2
	incr b
}
incr ab
}

set max1 0
set kkk1 0
set dlist1 [llength $distancelist2]
set ma1 0

while {$kkk1<$dlist1} {

set ma1 [lindex $distancelist2 $kkk1]

if {$ma1>$max1} {

set max1 $ma1
incr kkk1

} else {

incr kkk1

}

}

for {set dl 0} {$dl<$control} {incr dl} {
}

	    set ppp 1

	    set cx 0
	    set cy 0
	    set cz 0
	    set aggcx 0
	    set aggcy 0
	    set aggcz 0
set coordinate_p0 [part [lindex $alist 0] print folded_position]
set j 0
for {set c 0} {$c<3} {incr c} {
set coordinate($c) 0
set new_coordinate($c) 0
set total_coordinate($c) [lindex $coordinate_p0 $c]
}
set center [list]

set distlists(0) [list]
set distlists(1) [list]
set distlists(2) [list]
set bydlists 0

	    while {$ppp<$control} {
		set p1 [lindex $alist $ppp]
		set coordinate_p1 [part $p1 print folded_position]
	
	set dist [expr [lindex $coordinate_p0 0]-[lindex $coordinate_p1 0]]
	lappend distlists(0) $dist
	set dist [expr [lindex $coordinate_p0 1]-[lindex $coordinate_p1 1]]
	lappend distlists(1) $dist
	set dist [expr [lindex $coordinate_p0 2]-[lindex $coordinate_p1 2]]
	lappend distlists(2) $dist

incr ppp
    }

for {set coorj 0} {$coorj < 3} {incr coorj} {
	set l_dlists [llength $distlists($coorj)]
set coordinates $coorj
set kj 0
set result_logic 1
while {$kj<$l_dlists} {
	set bydlists [lindex $distlists($coorj) $kj]
	if {abs($bydlists)>[expr $box_length/2]} {
		set logic 0
		set result_logic [expr $logic && $result_logic]
	} else {
		set logic 1
		set result_logic [expr $logic && $result_logic]
		}
if {$result_logic == 1} {
	incr kj
} else {
	set kj $l_dlists
}
}
if {$result_logic==0} {
	set aggc($coorj) [com_for_new_positions $coordinates $box_length $alist]
} else {
	set aggc($coorj) [com_for_own_positions $coordinates $box_length $alist]
}
}

for {set jjj 0} {$jjj < 3} {incr jjj} {
	if {$aggc($jjj) > $box_length} {
		set aggc($jjj) [expr $aggc($jjj)-$box_length]
	} elseif {$aggc($jjj) < 0} {
		set aggc($jjj) [expr $aggc($jjj)+$box_length]
		} else {
		set taggc($jjj) $aggc($jjj)
		}

}

set center "$aggc(0)	$aggc(1)	$aggc(2)"

	    set pr 0
	    set difcenterplist [list]

	    set ri 0

	    while {$pr<$control} {

		set p2 [lindex $alist $pr]

set difcenterp [min_img_bond_length $center [part $p2 print folded_position] $box_length]

set t_length 0
set n_incr 0 
while {$n_incr<[expr (($max1+$sigma)/$increase)-1]} {

set t_length [expr $initial_radius+($increase*$n_incr)]
if {$difcenterp<$t_length} {

lappend ListForDf($t_length) $p2

} else {
set ListForDf($t_length) [list]
}

incr n_incr 
}

		set difsquare [expr $difcenterp*$difcenterp]
		set ri [expr $ri+$difsquare]
		incr pr
	    }

	    set s [expr 1.0/$control]
	    set sc [expr $s*$ri]
	    set radius_of_gyration [expr pow($sc,1.0/2.0)]

}
}

    writevcf $fl pids $pidlist

    incr i


    # CLOSE MAIN LOOP
}

exit

