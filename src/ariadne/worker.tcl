# ariadne worker.
#
# Runs the stock Pathways plugin over an explicit list of trajectory frames.
# Makes no decisions and computes no derived quantities: it drives the plugin
# and preserves the raw artifacts for Python to parse.
#
# Invoked as: vmd -dispdev text -e worker.tcl -args <chunkdir>

# `-dispdev text` removes DISPLAY from VMD's Tcl env array outright, whatever
# the parent environment holds -- so that flag, not the parent env, is what
# actually forces the plugin's headless path. This guard therefore fires only
# when the worker is launched in graphics mode, which is the misuse worth
# catching: there the plugin renders per frame (pathways.tcl:833-835, 979-991).
if { [info exists env(DISPLAY)] } {
    puts "WORKER) ERROR: running in graphics mode; refusing to run."
    puts "WORKER)        Launch with -dispdev text. In graphics mode the plugin"
    puts "WORKER)        renders every frame, which is slow and unreliable."
    exit 1
}

set chunkdir [lindex $argv 0]
if { $chunkdir eq "" } {
    puts "WORKER) ERROR: no chunk directory given"
    exit 1
}
source [file join $chunkdir "params.tcl"]

# Record the donor and acceptor atom indices (1-based) within the pruned
# subsystem, so Python can replicate the plugin's terminal trimming.
proc write_da_indices { pdb outfile dstr astr } {
    set m [mol new $pdb type pdb waitfor all]
    set ds [atomselect $m $dstr]
    set as [atomselect $m $astr]
    set dlist {}
    foreach i [$ds get index] { lappend dlist [expr {$i + 1}] }
    set alist {}
    foreach i [$as get index] { lappend alist [expr {$i + 1}] }
    $ds delete
    $as delete
    mol delete $m
    set fh [open $outfile "w"]
    puts $fh $dlist
    puts $fh $alist
    close $fh
}

package require pathways

# Slice covering exactly the frames this worker owns.
set first [lindex $PT_FRAMES 0]
set last  [lindex $PT_FRAMES end]

mol new $PT_PRMTOP type parm7 waitfor all
set molid [molinfo top]
mol addfile $PT_NC type netcdf first $first last $last waitfor all molid $molid
puts "WORKER) Loaded frames $first-$last ([molinfo $molid get numframes] in memory)"

# With withh = 0 the plugin converts its hydrogen-bond list through selection
# strings of the form "segid $s and resid $r and name $n" (pathways.tcl:796,
# aip2dp/dp2aip). AMBER topologies carry no segid, so those become
# "segid {} and ..." which VMD cannot parse and every frame dies.
#
# setseg is the plugin's own documented remedy: it copies chain into segid
# wherever segid is empty (pathways.tcl:470-502). Applied only when segids are
# actually missing, so systems that already define them are left untouched.
# This relabels atoms; it does not touch coordinates, bonds or parameters.
set segsel [atomselect $molid all]
set seglist [lsort -unique [$segsel get segid]]
$segsel delete
if { [lsearch -exact $seglist ""] >= 0 } {
    puts "WORKER) Empty segids present; applying setseg (chain -> segid)"
    setseg $molid
    set segsel [atomselect $molid all]
    puts "WORKER) segids now: [lsort -unique [$segsel get segid]]"
    $segsel delete
}

set tsv [open [file join $chunkdir "frames.tsv"] "a"]
set wrote_ref 0

foreach frame $PT_FRAMES {

    set local [expr {$frame - $first}]
    set prefix [file join $chunkdir [format "f%06d" $frame]]
    set t0 [clock milliseconds]

    set status "ok"
    if { [catch {
        eval pathways [list \
            "-d" $PT_DONOR \
            "-a" $PT_ACCEPTOR \
            "-b" $PT_BRIDGE \
            "-mol" $molid \
            "-frame" $local \
            "-o" $prefix \
            "-debug" "1" \
            "-rnd" "none"] $PT_PARGS
    } err] } {
        puts "WORKER) ERROR on frame $frame: $err"
        set status "error"
    }

    # Count atoms in the pruned subsystem the plugin actually solved on.
    set natoms 0
    if { [file exists "${prefix}.pdb"] } {
        set fh [open "${prefix}.pdb" "r"]
        while { [gets $fh line] >= 0 } {
            if { [string match "ATOM*" $line] || [string match "HETATM*" $line] } {
                incr natoms
            }
        }
        close $fh
    }

    # Keep one reference copy of the frame-invariant artifacts.
    if { !$wrote_ref && [file exists "${prefix}.pdb"] } {
        file copy -force "${prefix}.pdb" [file join $chunkdir "ref.pdb"]
        if { [file exists "${prefix}.cb"] } {
            file copy -force "${prefix}.cb" [file join $chunkdir "ref.cb"]
        }
        write_da_indices [file join $chunkdir "ref.pdb"] \
                         [file join $chunkdir "ref.da"] \
                         $PT_DONOR $PT_ACCEPTOR
        set wrote_ref 1
    }

    # Discard the bulky per-frame artifacts; keep 1.out and .hb.
    foreach ext {".pdb" ".psf" ".cb" "2.out"} {
        file delete -force "${prefix}${ext}"
    }

    # The plugin leaves one extra molecule per call when headless
    # (the pruned clone at pathways.tcl:829/851).
    while { [llength [molinfo list]] > 1 } {
        mol delete top
    }

    set elapsed [expr {([clock milliseconds] - $t0) / 1000.0}]
    puts $tsv "$frame\t$natoms\t$status\t$elapsed"
    flush $tsv
    puts "WORKER) frame $frame done in ${elapsed}s (status $status)"
}

close $tsv

set fh [open [file join $chunkdir "DONE"] "w"]
puts $fh "ok"
close $fh

quit
