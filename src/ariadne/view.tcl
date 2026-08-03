# ARIADNE 3D pathway viewer.
#
# Loads only the analysed frames, aligns them on donor+acceptor so the
# endpoints stay still, and redraws the pathway whenever the frame changes.
#
# Invoked as: vmd [-dispdev text] -e view.tcl -args <viewparams.tcl>

set viewparams [lindex $argv 0]
if { $viewparams eq "" || ![file exists $viewparams] } {
    puts "VIEW) ERROR: viewparams file not given or missing"
    exit 1
}
source $viewparams

package require pathways

# ---------------------------------------------------------------- resolution

# Turn one path step into exactly one atom index.
#
# paths.csv records segid as X for runs made after the setseg fix and empty for
# earlier ones, while a freshly loaded molecule has neither until setseg runs.
# Matching an empty segid literally would select nothing and draw an invisible
# path, so the clause is omitted when the recorded segid is empty.
proc resolve_atom { molid entry } {
    foreach {resid resname name segid bond} $entry { break }
    set text "resid $resid and resname $resname and name $name"
    if { $segid ne "" } {
        append text " and segid $segid"
    }
    set sel [atomselect $molid $text]
    set n [$sel num]
    set idx [lindex [$sel get index] 0]
    $sel delete
    if { $n != 1 } {
        error "selection '$text' matched $n atoms, expected exactly 1"
    }
    return $idx
}

# ---------------------------------------------------------------- load

proc load_scene { } {
    global VIEW_PRMTOP VIEW_TRAJ VIEW_FRAMES VIEW_FROMEXPORT VIEW_RUNDIR

    if { $VIEW_FROMEXPORT } {
        set dir [file join $VIEW_RUNDIR "view"]
        mol new [file join $dir "view.psf"] type psf waitfor all
        set molid [molinfo top]
        mol addfile [file join $dir "view.dcd"] type dcd waitfor all molid $molid
        puts "VIEW) Loaded exported scene from $dir"
        return $molid
    }

    mol new $VIEW_PRMTOP type parm7 waitfor all
    set molid [molinfo top]
    foreach f $VIEW_FRAMES {
        mol addfile $VIEW_TRAJ type netcdf first $f last $f waitfor all molid $molid
    }
    puts "VIEW) Loaded [molinfo $molid get numframes] frame(s) of [llength $VIEW_FRAMES] requested"

    # Same conditional setseg as the worker, so segids agree on both sides.
    set s [atomselect $molid all]
    set seglist [lsort -unique [$s get segid]]
    $s delete
    if { [lsearch -exact $seglist ""] >= 0 } {
        puts "VIEW) Empty segids present; applying setseg (chain -> segid)"
        setseg $molid
    }
    return $molid
}

# ---------------------------------------------------------------- align

# Fit every frame onto frame 0 using the donor and acceptor termini, so those
# endpoints stay fixed and the route between them is what visibly moves. Done
# once at load; there is no per-frame cost.
#
# The fit is taken over `same residue as` the two selections rather than the
# selections themselves. Single-atom termini are the documented, shipped case
# (examples/fmn_heme_intra_hb_view.in: one FMN N5 and one heme FE), and two
# atoms cannot determine a rotation -- so fitting on the raw selections would
# print "alignment skipped" on every real run and the headline behaviour would
# never happen. A whole residue at each end is rigid enough for a
# well-determined fit while still pinning only the termini. The bare selection,
# and then the warning, remain as fallbacks for exotic cases.
proc align_scene { molid } {
    global VIEW_DONOR VIEW_ACCEPTOR

    set ends "($VIEW_DONOR) or ($VIEW_ACCEPTOR)"
    set text "same residue as ($ends)"
    set ref [atomselect $molid $text frame 0]
    if { [$ref num] < 3 } {
        $ref delete
        set text $ends
        set ref [atomselect $molid $text frame 0]
    }
    if { [$ref num] < 3 } {
        puts "VIEW) WARNING: only [$ref num] atoms match '$text'; alignment skipped"
        $ref delete
        return 0
    }
    set nfit [$ref num]
    set n [molinfo $molid get numframes]
    for { set i 1 } { $i < $n } { incr i } {
        set cur [atomselect $molid $text frame $i]
        set all [atomselect $molid all frame $i]
        $all move [measure fit $cur $ref]
        $cur delete
        $all delete
    }
    $ref delete
    puts "VIEW) Aligned $n frame(s) on $nfit atom(s) of '$text'"
    return 1
}

# ---------------------------------------------------------------- export

proc export_scene { molid } {
    global VIEW_RUNDIR VIEW_BRIDGE VIEW_DONOR VIEW_ACCEPTOR VIEW_FRAMES

    # view.dcd frame i is meaningless without knowing which trajectory frame it
    # came from: --from-export draws bundle frame i with VIEW_PATH for the i-th
    # entry of the current --frames, so a bundle exported at 0:999:10 and
    # reloaded at 0:999:20 would put every cylinder on the wrong coordinates
    # with nothing to show for it. Record the mapping next to the DCD; view.py
    # refuses to reload a bundle whose frames.txt disagrees.
    #
    # If the loaded frame count and VIEW_FRAMES ever disagree the manifest
    # would itself be a lie, so refuse rather than write it. Uses the same
    # SIGKILL as the resolve loop, for the same reason: VMD's exit reports 0.
    set nloaded [molinfo $molid get numframes]
    if { $nloaded != [llength $VIEW_FRAMES] } {
        puts "VIEW) ERROR: loaded $nloaded frame(s) but [llength $VIEW_FRAMES] were\
              requested; refusing to write a bundle whose frames.txt would be wrong"
        flush stdout
        exec kill -9 [pid]
    }

    set dir [file join $VIEW_RUNDIR "view"]
    file mkdir $dir
    set sel [atomselect $molid "($VIEW_DONOR) or ($VIEW_ACCEPTOR) or ($VIEW_BRIDGE)"]
    $sel writepsf [file join $dir "view.psf"]
    $sel writepdb [file join $dir "view.pdb"]
    animate write dcd [file join $dir "view.dcd"] waitfor all sel $sel $molid
    puts "VIEW) Exported [$sel num] atoms x [molinfo $molid get numframes] frames to $dir"
    $sel delete

    set fh [open [file join $dir "frames.txt"] "w"]
    foreach f $VIEW_FRAMES {
        puts $fh $f
    }
    close $fh

    set fh [open [file join $dir "view.vmd"] "w"]
    puts $fh "mol new view.psf type psf waitfor all"
    puts $fh "mol addfile view.dcd type dcd waitfor all"
    close $fh
}

# ---------------------------------------------------------------- main

set molid [load_scene]
align_scene $molid

# Resolve every atom up front so a bad selection fails loudly here rather than
# producing an empty picture later.
#
# A bare Tcl `error` does not reliably yield a non-zero process exit code under
# `vmd -e`: VMD's own "exit" command is a built-in that runs its shutdown
# sequence and always terminates the process with status 0 ("Exiting
# normally"), regardless of the value passed to it -- confirmed by direct
# experiment, not assumption. So on failure this catches the error, prints a
# clear VIEW)-prefixed line, flushes it, and self-delivers SIGKILL: that is
# the one mechanism that reliably produces a non-zero exit code from this VMD
# build, since no Tcl-level exit call can override VMD's own exit handler.
set resolved 0
set resolve_err ""
if { [catch {
    foreach f $VIEW_FRAMES {
        if { ![info exists VIEW_RANKLIST($f)] } { continue }
        foreach r $VIEW_RANKLIST($f) {
            foreach entry $VIEW_PATH($f,$r) {
                resolve_atom $molid $entry
                incr resolved
            }
        }
    }
} resolve_err] } {
    puts "VIEW) ERROR: $resolve_err"
    flush stdout
    exec kill -9 [pid]
}
puts "VIEW) resolved $resolved path atoms"

# INVARIANT: nothing above this line may write a file. The resolve loop above
# ends a failure by SIGKILLing this process, so any output already opened would
# be left half-written with no chance to clean up -- and a truncated view/
# bundle looks loadable. Keep every file write below this point, and if the
# catch above is ever widened, move the writes further down rather than
# relaxing this rule.
if { $VIEW_EXPORT } {
    export_scene $molid
}

# ---------------------------------------------------------------- drawing

# VMD colour ids: blue 0, red 1, grey 2, green 7.
# Colour indices come from viewparams.tcl so they stay configurable from the
# input file. Defaults are warm (orange / light orange / yellow) because VMD
# colours atoms by element, and a blue or red path competes with the nitrogen
# and oxygen of the licorice underneath it.
proc bond_color { bond } {
    global VIEW_COL_COVALENT VIEW_COL_HBOND VIEW_COL_THROUGH_SPACE
    switch -- $bond {
        covalent      { return $VIEW_COL_COVALENT }
        hbond         { return $VIEW_COL_HBOND }
        through_space { return $VIEW_COL_THROUGH_SPACE }
        default       { return 2 }
    }
}

# Atom indices of the dominant path on the first rendered frame. Used for the
# context shell, the clipping plane and the camera, so they all agree.
proc path_indices { molid } {
    global VIEW_FRAMES VIEW_PATH VIEW_RANKLIST
    set f [lindex $VIEW_FRAMES 0]
    if { ![info exists VIEW_RANKLIST($f)] } { return {} }
    if { [lsearch -exact $VIEW_RANKLIST($f) 0] < 0 } { return {} }
    set idxs {}
    foreach entry $VIEW_PATH($f,0) {
        lappend idxs [resolve_atom $molid $entry]
    }
    return $idxs
}

proc setup_reps { molid } {
    global VIEW_DONOR VIEW_ACCEPTOR VIEW_CONTEXT_RADIUS
    mol delrep 0 $molid

    # Context cartoon, restricted to a shell around the pathway. Showing the
    # whole protein buries the path: the transparent ribbons in front of it
    # veil the very thing you are looking at. A radius of 0 drops the cartoon
    # entirely; a large radius restores the whole protein.
    set idxs [path_indices $molid]
    if { $VIEW_CONTEXT_RADIUS > 0 } {
        if { [llength $idxs] } {
            set context "protein and same residue as (within $VIEW_CONTEXT_RADIUS of (index $idxs))"
        } else {
            set context "protein"
        }
        mol representation NewCartoon 0.3 12.0 4.1 0
        mol color ColorID 2
        mol selection $context
        mol material Transparent
        mol addrep $molid
    }

    # Donor and acceptor, always visible.
    mol representation Licorice 0.2 12.0 12.0
    mol color Name
    mol selection "($VIEW_DONOR) or ($VIEW_ACCEPTOR)"
    mol material Opaque
    mol addrep $molid

    # Residues on the current frame's dominant path; selection updated per frame.
    mol representation Licorice 0.15 12.0 12.0
    mol color Name
    mol selection "none"
    mol material Opaque
    mol addrep $molid
    return [expr {[molinfo $molid get numreps] - 1}]
}

# Label a frame that has no pathway, in the scene as well as on stdout: a
# structure with no cylinders and no explanation is indistinguishable from a
# bug, which is exactly the silent-wrong-picture failure this viewer must not
# have. Placed at the donor/acceptor midpoint, which is always on screen.
proc draw_frame_label { molid frame text } {
    global VIEW_DONOR VIEW_ACCEPTOR

    set sel [atomselect $molid "($VIEW_DONOR) or ($VIEW_ACCEPTOR)" frame $frame]
    if { [$sel num] == 0 } {
        $sel delete
        return
    }
    set c [measure center $sel]
    $sel delete
    graphics $molid color 1
    graphics $molid text $c $text size 1.2
}

# Redraw every path for the frame currently displayed.
proc redraw_paths { molid pathrep } {
    global VIEW_FRAMES VIEW_RADIUS VIEW_RANKS VIEW_PATH VIEW_RANKLIST VIEW_STATUS

    graphics $molid delete all

    set i [molinfo $molid get frame]
    set f [lindex $VIEW_FRAMES $i]
    if { ![info exists VIEW_RANKLIST($f)] } {
        # Computed but pathless (status no_path or error): show the structure
        # and say so, rather than presenting an empty scene as if it were one.
        mol modselect $pathrep $molid "none"
        set why "no path"
        if { [info exists VIEW_STATUS($f)] } {
            set why "no path (status $VIEW_STATUS($f))"
        }
        puts "VIEW) frame $f: $why"
        catch { draw_frame_label $molid $i "frame $f: $why" }
        return
    }

    set resids {}
    foreach r $VIEW_RANKLIST($f) {
        if { $VIEW_RANKS eq "0" && $r != 0 } { continue }

        if { $r == 0 } {
            set rad $VIEW_RADIUS
            set mat "Opaque"
        } else {
            set rad [expr {$VIEW_RADIUS / 3.0}]
            set mat "Transparent"
        }

        set prev ""
        foreach entry $VIEW_PATH($f,$r) {
            set idx [resolve_atom $molid $entry]
            set sel [atomselect $molid "index $idx" frame $i]
            set xyz [lindex [$sel get {x y z}] 0]
            $sel delete
            if { $r == 0 } { lappend resids [lindex $entry 0] }

            # Only the previous atom's position is needed: the cylinder's
            # colour and style come from the bond type recorded on the *current*
            # step, which is the bond that brought the path here.
            if { $prev ne "" } {
                set pxyz $prev
                set bond [lindex $entry 4]
                if { $r == 0 } {
                    set col [bond_color $bond]
                } else {
                    set col 2
                }
                # draw_step treats type CB as a solid cylinder and anything
                # else as 0.5A dashes (pathways.tcl:342-361).
                if { $bond eq "covalent" } { set type "CB" } else { set type "TS" }
                draw_step $molid $pxyz $xyz $rad 20 $col $mat $type
            }
            set prev $xyz
        }
    }

    if { [llength $resids] } {
        mol modselect $pathrep $molid "resid [lsort -unique -integer $resids]"
    } else {
        mol modselect $pathrep $molid "none"
    }
    puts "VIEW) frame $f drawn"
}

# Pivot the camera about the dominant path rather than the box centre.
proc center_on_path { molid } {
    global VIEW_FRAMES VIEW_PATH VIEW_RANKLIST

    set f [lindex $VIEW_FRAMES 0]
    if { ![info exists VIEW_RANKLIST($f)] } { return }
    set idxs {}
    foreach entry $VIEW_PATH($f,0) {
        lappend idxs [resolve_atom $molid $entry]
    }
    if { ![llength $idxs] } { return }
    set sel [atomselect $molid "index $idxs" frame 0]
    set c [measure center $sel]
    set mm [measure minmax $sel]
    $sel delete

    # Zoom to the pathway, NOT to the molecule. `display resetview` frames
    # every loaded atom, and 92% of this system is solvent (162,702 of 177,430
    # atoms here), so it leaves the protein and the path as a speck inside the
    # periodic water box -- confirmed on screen before this was changed.
    #
    # VMD's viewport spans roughly 2 world units, so scale ~ 2/extent fits an
    # object exactly; 1.4 leaves a margin around the path.
    set extent [veclength [vecsub [lindex $mm 1] [lindex $mm 0]]]
    if { $extent <= 0.0 } { set extent 10.0 }

    molinfo $molid set center_matrix [list [transoffset [vecinvert $c]]]
    molinfo $molid set rotate_matrix [list [transidentity]]
    molinfo $molid set global_matrix [list [transidentity]]
    scale to [expr {1.4 / $extent}]
    puts [format "VIEW) Centred on the pathway: %.1f A across" $extent]
}

# Hide the protein in FRONT of the pathway, keeping what is behind it as
# context. Without this the transparent ribbons drift across the very thing
# you are trying to look at.
#
# Sign convention, established by rendering both directions and comparing
# against a cartoon-off reference: a clipping plane with normal N keeps the
# geometry on the side N points TOWARDS. So to keep what is behind, the normal
# points away from the camera.
#
# The plane is expressed in molecule coordinates and therefore rotates with the
# molecule. After rotating the view, type `reclip` in the VMD console to
# re-align it with the new viewing direction.
proc clip_front { molid } {
    global VIEW_CLIP_FRONT VIEW_CONTEXT_RADIUS

    if { !$VIEW_CLIP_FRONT || $VIEW_CONTEXT_RADIUS <= 0 } { return 0 }
    set idxs [path_indices $molid]
    if { ![llength $idxs] } { return 0 }

    # Which way is the camera, in molecule coordinates? The viewer looks along
    # -z in view space, so undo the molecule's rotation to find it. With the
    # rotation identity (as set at startup) this is simply {0 0 1}.
    set rot [lindex [molinfo $molid get rotate_matrix] 0]
    if { [catch {set toward [coordtrans [measure inverse $rot] {0 0 1}]}] } {
        set toward {0 0 1}
    }
    set toward [vecnorm $toward]

    set sel [atomselect $molid "index $idxs" frame [molinfo $molid get frame]]
    set centre [measure center $sel]
    set front -1.0e30
    foreach p [$sel get {x y z}] {
        set d [vecdot $p $toward]
        if { $d > $front } { set front $d }
    }
    $sel delete

    # Sit the plane just in front of the frontmost path atom, so no part of the
    # pathway region is cut away.
    set shift [expr {$front - [vecdot $centre $toward] + 1.0}]
    set plane [vecadd $centre [vecscale $toward $shift]]

    mol clipplane center 0 0 $molid $plane
    mol clipplane normal 0 0 $molid [vecinvert $toward]
    mol clipplane status 0 0 $molid 1
    puts "VIEW) Clipped the cartoon in front of the pathway (type 'reclip' after rotating)"
    return 1
}

# Convenience for the VMD console: re-align the clipping plane after rotating.
proc reclip { } {
    clip_front [molinfo top]
}

# Register the frame-change redraw. VMD's `apply` lambda form is the modern
# idiom for `trace add variable ... write`, but its exact callback signature
# (which extra args Tcl appends) has varied across VMD builds. Try that form
# first; if registering it errors (rather than firing it -- registration
# itself can raise on some builds), fall back to a named global proc that
# takes molid and pathrep as ordinary prepended arguments and swallows
# whatever extra args the trace appends, instead of relying on `apply`'s
# argument-passing behaviour.
proc _view_redraw_trace_fallback { molid pathrep args } {
    redraw_paths $molid $pathrep
}

if { [info exists env(DISPLAY)] } {
    set pathrep [setup_reps $molid]
    center_on_path $molid
    clip_front $molid
    redraw_paths $molid $pathrep

    # Fires on slider drags, arrow keys and playback alike.
    global vmd_frame
    if { [catch {
        trace add variable vmd_frame($molid) write \
            [list apply {{molid pathrep args} {redraw_paths $molid $pathrep}} $molid $pathrep]
    } trace_err] } {
        puts "VIEW) WARNING: apply-lambda trace registration failed ($trace_err); using named-proc fallback"
        trace add variable vmd_frame($molid) write \
            [list _view_redraw_trace_fallback $molid $pathrep]
    }
    # Start at the beginning rather than on the last-loaded frame, and make
    # sure the Main window is up: that is where the frame slider and the
    # animation speed control live, which is how you scrub and slow playback.
    menu main on
    animate goto 0
    animate style Loop
    animate speed 0.4

    puts "VIEW) Interactive."
    puts "VIEW)   Scrub frames:  drag the slider in the VMD Main window"
    puts "VIEW)   Play / pause:  the arrow buttons beside it"
    puts "VIEW)   Speed:         the 'speed' slider in VMD Main (set to 0.4 of max)"
    puts "VIEW)   Step one frame: the single-arrow buttons, or +/- keys in the display"
    puts "VIEW)   Quit:          File > Quit, or type 'quit' in this terminal"
} else {
    puts "VIEW) Headless; no representations built"
}

puts "VIEW) READY"
