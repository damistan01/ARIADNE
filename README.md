# ARIADNE

![platform: Linux](https://img.shields.io/badge/platform-Linux-blue)
![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)

> "However blameless the life we lead, the passions and the greed of men could
> bring us to ruin, and there was nothing we could do."
>
> — Jennifer Saint, *Ariadne* (2021)

**A**utomated **R**ate & **I**nteratomic-coupling **A**nalysis of
**D**onor–acceptor **N**onadiabatic **E**lectron-transfer

Batch electron-transfer pathway analysis over molecular dynamics trajectories.

**Linux only.** ARIADNE drives a Linux build of VMD and its Pathways plugin,
and assumes Unix process and filesystem behaviour throughout. On Windows, run
it inside WSL2 — that is how it is developed. See
[Platform support](#platform-support).

**New here?** Install (below), then follow
**[the tutorial](tutorial/README.md)** — a complete ten-frame run with real
output, end to end in a few minutes.

---

## What this does, in plain terms

Electrons tunnel through proteins along preferred routes. The
[Pathways plugin](#credits) for VMD (Balabin, Hu & Beratan) finds those routes
in **one** structure: you give it a donor atom and an acceptor atom, and it
returns the best tunnelling path between them and the electronic coupling
`T_DA` along it.

A single structure is a snapshot, and proteins move. To say anything
statistical — a mean coupling `<T_DA>`, its fluctuation, the coherence
parameter — you need the same calculation on hundreds or thousands of frames of
an MD trajectory. Doing that by hand in the VMD GUI is a day of clicking.

**ARIADNE runs it for you.** One parameter file, one command, every frame, in
parallel, with the results written out as CSVs you can plot, plus a full record
of exactly what was run.

It does **not** reimplement the physics. Bond perception, hydrogen-bond
detection, pruning and the `pathcore` call all happen inside the unmodified
Pathways 1.4 plugin, so a frame processed by ARIADNE and a frame you run by
hand in the VMD GUI traverse identical code. This was verified against a
hand-driven plugin run, whose coupling the batch pipeline reproduces to 1e-12
relative.

Named for the thread that traced a route through the labyrinth — and the branch
it lives on, `clew`, for the ball of thread itself.

### What you need to bring

| You provide | Example |
| --- | --- |
| An AMBER topology (`.prmtop` / `.parm7`) | `system.prmtop` |
| An AMBER NetCDF trajectory (`.nc`) | `system_md.nc` |
| A donor and an acceptor, as VMD atom selections | `resname FMN and name N5` |

### What you get back

A directory containing `couplings.csv` (one row per frame), `paths.csv` (one
row per step of every path), `summary.txt` (the statistics) and `run.json`
(complete provenance). And, optionally, a 3D animation of the pathway moving
through the protein.

---

## Platform support

| Platform | Status |
| --- | --- |
| Linux (x86-64) | Supported. Developed and tested on Ubuntu under WSL2. |
| Windows | Via WSL2 only — install a Linux distribution and follow this README inside it. Not supported natively. |
| macOS | Not supported. Untested, and the plugin layout and `pathcore` binary assumed here are Linux builds. |

Nothing in the Python is deliberately Linux-specific, but the toolchain it
drives is, and no other platform is tested.

---

## Installation

Five steps. Steps 1–3 are the VMD side and are the fiddly part; if you already
run the Pathways plugin from the VMD GUI, skip to step 4.

### 1. Python 3.11 or newer

ARIADNE uses only the Python standard library — nothing to `pip install` to run
it.

```bash
python3 --version        # must print 3.11 or higher
```

If it is older, install a newer Python from your distribution, e.g. on
Debian/Ubuntu:

```bash
sudo apt update && sudo apt install python3.11
```

### 2. VMD

Download the **Linux (x86-64)** build of VMD from
<https://www.ks.uiuc.edu/Research/vmd/> (registration is free) and install it,
typically:

```bash
tar xzf vmd-*.opengl.tar.gz
cd vmd-*
./configure
cd src && sudo make install
```

Check it is on your `PATH`:

```bash
vmd -h        # should print VMD's usage text
```

If you get `vmd: command not found`, add its install directory to `PATH` in
your `~/.bashrc` or `~/.zshrc`:

```bash
export PATH="$PATH:/path/to/vmd/bin"
```

### 3. The Pathways 1.4 plugin and `pathcore`

The Pathways plugin is a separate download from the Beratan group. It ships as
a directory of Tcl (the plugin itself) plus `pathcore`, a compiled binary that
does the actual search.

Put the plugin directory somewhere permanent, for example:

```
/home/you/tools/vmd-2.0.0/plugins/LINUXAMD64/tcl/pathways1.4
```

Then tell VMD where to find it by adding this line to `~/.vmdrc` (create the
file if it does not exist), naming the **parent** directory of `pathways1.4`:

```tcl
lappend auto_path /home/you/tools/vmd-2.0.0/plugins/LINUXAMD64/tcl
```

Put `pathcore` on your `PATH` — a symlink into a directory already on `PATH` is
the simplest way:

```bash
mkdir -p ~/bin
ln -s /path/to/pathways1.4/pathcore ~/bin/pathcore
export PATH="$HOME/bin:$PATH"        # add this to ~/.bashrc too
```

Verify both:

```bash
which pathcore                       # prints a path
echo 'package require pathways' | vmd -dispdev text
```

The second command should load without a `can't find package pathways` error.

### 4. Get ARIADNE

```bash
git clone https://github.com/damistan01/ARIADNE.git
cd ARIADNE
```

There is nothing to build or install. Run `bin/ariadne` directly, or put it on
your `PATH`:

```bash
ln -s "$PWD/bin/ariadne" ~/bin/ariadne
```

### 5. Check the install

```bash
ariadne -i tutorial/pathways.in \
    -p system.prmtop -y system_md.nc -o /tmp/check/ --dry-run
```

`--dry-run` runs every preflight check — VMD, `pathcore`, the plugin, your
selections, the frame range — and prints the plan without computing anything.
If it exits 0, you are ready.

---

## Quickstart

### Step 1 — write an input file

The input file holds both the physics parameters and (later) the display
settings. A minimal one is three lines:

```
donor    = resname FMN and name N5
acceptor = resname FE1 and name FE
bridge   = all
```

`donor`, `acceptor` and `bridge` are ordinary **VMD atom selections** — the
same syntax as the selection box in the VMD GUI. If you are unsure what to
write, open your structure in VMD and try the selection there first; it must
match exactly the atoms you mean.

Everything else has a default. The two you will most likely set:

```
npaths   = 5    # how many ranked paths to return per frame (default 1)
withh    = 1    # see the note below
```

`withh` reads backwards from what you might expect. `withh = 1` means
"hydrogens are explicitly present in the structure", and in that mode the
plugin deliberately writes an **empty** hydrogen-bond list, so no H-bond jumps
are considered at all. `withh = 0` strips hydrogens and runs `measure hbonds`,
making H-bond-mediated steps available. **To enable hydrogen bonds, set
`withh = 0`.**

[`tutorial/pathways.in`](tutorial/pathways.in) is a complete, commented
example you can copy. [`tutorial/README.md`](tutorial/README.md) walks through
a full ten-frame run using it, with real output.

<details>
<summary>Full list of physics keys</summary>

`donor`, `acceptor`, `bridge`, `npaths`, `epsc`, `epsh`, `exph`, `r0h`,
`epsts`, `expts`, `r0ts`, `hcut`, `hang`, `tscut`, `procut`, `withh`, `cda`.

These map one-to-one onto the Pathways plugin's own parameters and keep the
plugin's defaults; consult the plugin documentation for their meaning.
</details>

### Step 2 — run it

```bash
ariadne -i tutorial/pathways.in \
    -p system.prmtop \
    -y system_md.nc \
    -o pathways_WT/ \
    -n 4
```

ARIADNE runs one real frame before fanning out, so a bad selection or a missing
`pathcore` fails in seconds rather than after eight processes have each loaded
a trajectory slice. That probe also prints a measured ETA — read it before
walking away, because a full trajectory can take hours (see
[Measured performance](#measured-performance)).

If a run is interrupted, `--resume` picks it up where it stopped.

### Step 3 — read the output

`pathways_WT/` will contain:

**`summary.txt`** — start here. The statistics over all successful frames:

```
PATHWAY COUPLING STATISTICS

  <T_DA>:        ...      mean coupling
  sigma(T_DA):   ...      its standard deviation across frames
  <T_DA^2>:      ...      mean square coupling (this is what rate theory wants)
  coherence:     ...      <T_DA>^2 / <T_DA^2>, between 0 and 1

FRAMES

  ok: 1000
```

The **coherence parameter** near 1 means the coupling is dominated by a single
persistent path; near 0 means it is dominated by fluctuations, and no one
structure represents the system.

**`couplings.csv`** — one row per frame:

| Column | Meaning |
| --- | --- |
| `frame` | 0-based frame index in the trajectory |
| `time_ps` | simulation time of that frame, from the NetCDF file |
| `t_da` | coupling of the best path in that frame |
| `log10_t_da` | its base-10 logarithm, for plotting |
| `n_paths` | how many paths `pathcore` returned |
| `n_steps` | number of steps in the best path |
| `status` | `ok`, or an error marker |

**`paths.csv`** — one row per step of every path, i.e. the routes themselves:

| Column | Meaning |
| --- | --- |
| `frame` | frame index |
| `path_rank` | 0 is the dominant path, 1 the next best, … |
| `t_da` | coupling of *that* path |
| `step` | position along the path, starting at the donor |
| `atom_index` | 0-based index into the full topology |
| `resid`, `resname`, `atom_name`, `segid` | which atom it is |
| `bond_type` | `covalent`, `hbond` or `through_space` |

**`run.json`** — every input path and checksum, the VMD and plugin versions,
the resolved parameters, the frames requested and completed, and wall time.
Keep it: it is what makes a run reproducible.

### Step 4 — look at it in 3D

```bash
ariadne view pathways_WT/
```

See [Viewing pathways in 3D](#viewing-pathways-in-3d).

---

## Command reference

```
ariadne -i INPUT -p TOPOLOGY -y TRAJECTORY -o OUTDIR [options]
```

| Flag | Meaning |
| --- | --- |
| `-i` | input file with the Pathways parameters |
| `-p` | AMBER topology (parm7) |
| `-y` | AMBER NetCDF trajectory |
| `-o` | output directory (created) |
| `--frames start:stop:stride` | frame selection, 0-based, stop inclusive; overrides the `frames` key in the input file |
| `-n` | parallel workers (default: min(8, cores)) |
| `--water` | keep water in the bridge selection |
| `--resume` | continue an interrupted run, or widen the frame range of a completed one — frames already computed are kept, not recomputed |
| `--keep-raw` | keep the scratch directory |
| `--dry-run` | preflight and plan only |

Exit codes: 0 success, 1 usage or config error, 2 preflight failure,
3 one or more workers failed.

---

## Viewing pathways in 3D

```bash
ariadne view pathways_WT/
```

Opens VMD showing the pathway across the analysed frames, with the donor and
acceptor held stationary so only the route between them moves. Use the frame
slider or press play, exactly as with any trajectory.

This needs a graphical display. Under WSL2, current Windows builds provide one
(WSLg) with no setup; over SSH, use `ssh -X`. With no display at all, use
`--text --export` to build the bundle and view it elsewhere.

| Flag | Meaning |
| --- | --- |
| `--export` | also write a self-contained `view/` bundle (~1.7 MB for 10 frames) |
| `--from-export` | load `view/`, so the topology and the 2.1 GB trajectory need not be present. The render selection must match the one the bundle was exported with; `view/frames.txt` records it and a mismatch is refused |
| `--frames` | override `view_frames` for this invocation |
| `--text` | run headless: validate and export, display nothing |

Presentation is configured in the same input file as the physics:

| Key | Default | Meaning |
| --- | --- | --- |
| `frames` | all | which frames to **compute**; `--frames` overrides it |
| `view_frames` | all computed | which computed frames to render; an open-ended stop (`0:`) means the last computed frame |
| `view_radius` | 0.3 | cylinder radius of the dominant path, in Å |
| `view_ranks` | `all` | `0` for the dominant path only, `all` to draw lower ranks faintly |
| `view_color_covalent` | `orange` | colour of covalent steps |
| `view_color_hbond` | `orange3` | colour of hydrogen-bond steps (light orange) |
| `view_color_through_space` | `yellow` | colour of through-space jumps |
| `view_context_radius` | 12 | Å shell of protein cartoon kept around the path; `0` hides the cartoon entirely |
| `view_clip_front` | `yes` | hide the cartoon that sits in front of the pathway |

The colour defaults are warm on purpose. VMD colours atoms by element, so
nitrogen is blue and oxygen red; drawing the path in those colours makes it
compete with the licorice it runs through. Orange -> light orange -> yellow
reads as one gradient against a grey cartoon and collides with nothing.

Colours accept any VMD colour name, or its index:

```
blue red gray orange yellow tan silver green white pink cyan purple
lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2
cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2
orange3
```

Note that `orange3` is the *lighter* orange (0.96 0.72 0.00); `orange2` is
darker (0.89 0.35 0.00).

### Seeing past the protein

Drawn naively, the transparent cartoon veils the very thing you are looking at:
ribbons in front of the pathway wash it out. Two settings deal with this, and
both are on by default.

`view_context_radius` keeps only the protein within a shell of the path, which
removes ribbons that merely drift across the view from elsewhere in the
structure. `view_clip_front` then clips away whatever cartoon still sits
between the camera and the path, leaving the structure behind it as context.

Measured on the reference run, mean rendered brightness -- a proxy for how much
veil is in the way -- falls from 0.062 with the whole protein to 0.030 with
both settings on, against 0.016 for no cartoon at all.

The clipping plane is defined in molecule coordinates, so it rotates with the
structure. After rotating the view, type `reclip` in the VMD console to
re-align it with the new viewing direction.

Presentation keys are recorded in `run.json` but excluded from the `--resume`
guard, so changing a radius never blocks resuming a long run.

---

## Measured performance

Measured on 1KBI_MD/WT (1000 frames, 177,430 atoms; pruned subsystem 14,728
atoms) with the tutorial input at `npaths = 5`, on 8 cores / 15 GB.

| Workers | Per-frame when concurrent | Slowdown vs. alone | Effective | 1000 frames |
| --- | --- | --- | --- | --- |
| 1 | 73 s | 1.0x | 73 s/frame | ~20.3 h |
| **4** | 144 s | 2.0x | **36.0 s/frame** | **~10.0 h** |
| 8 | 435 s | 6.0x | 54.4 s/frame | ~15.1 h |

**Use `-n 4`, not `-n 8`.** `pathcore` is single-threaded, so this is not CPU
oversubscription -- it is memory/cache contention, and it grows superlinearly.
Going from 4 to 8 workers triples the per-frame cost, which more than cancels
the extra parallelism.

`npaths` is the other cost dial: `npaths = 1` costs 22.4 s/frame on a single
worker against 73 s/frame for `npaths = 5`, so about 3.3x. Use `npaths = 1` if
you only need `<T_DA>` and the coherence parameter.

Stride is a bigger lever than either. Frames in these trajectories are 100 ps
apart, so `--frames 0:999:2` halves the run for nearly the same statistics if
the coupling autocorrelation time is well under 1 ns.

---

## Troubleshooting

| Message | What it means |
| --- | --- |
| `ariadne: vmd is not on PATH` | VMD is not installed, or its `bin` directory is not in `PATH`. See [step 2](#2-vmd). |
| `ariadne: pathcore is not on PATH` | The `pathcore` binary was not found. See [step 3](#3-the-pathways-14-plugin-and-pathcore). |
| `can't find package pathways` | `~/.vmdrc` does not `lappend auto_path` the directory *containing* `pathways1.4`. Point it at the parent directory, not at `pathways1.4` itself. |
| `ariadne: no DISPLAY available; use --text --export instead` | `ariadne view` needs a graphical display. Under WSL2 use a Windows build with WSLg; over SSH use `ssh -X`; otherwise export the bundle and view it on a machine that has one. |
| `line N: unknown key 'foo'` | A typo in the input file. Keys are listed in [Quickstart](#step-1--write-an-input-file) and in the viewer table. |
| `donor must not be empty`, or a selection matching nothing | The VMD selection does not match your topology. Open the structure in VMD and test the selection interactively. |
| `frame range ... exceeds the trajectory` | `--frames` asks for frames the NetCDF file does not have. Frames are 0-based and `stop` is inclusive. |
| Exit code 3 | One or more workers failed. `couplings.csv` marks which frames, and `worker.log` in the scratch directory says why — re-run with `--keep-raw` to preserve it. |

---

## Development

The pytest suite and the design notes are not distributed with the repo: both
pin absolute paths to the reference trajectory and toolchain on the development
machine, so neither is meaningful against anyone else's data as written.

The code is laid out as one module per stage of the pipeline —
`config.py` (parse the input file), `plan.py` (decide the frames and chunks),
`runner.py` (drive the VMD workers), `parse.py` (read what they produced),
`writer.py` (emit the CSVs) — plus `view.py` and `viewconfig.py` for the 3D
viewer, and the two Tcl scripts `worker.tcl` and `view.tcl` that run inside
VMD.

---

## Using and contributing

ARIADNE is MIT-licensed: clone it, change it, and use it for whatever you like,
privately or publicly, without asking anyone. Local modifications are yours to
keep.

If you would like to send changes back, fork the repository and open a pull
request against `clew` — see [CONTRIBUTING.md](CONTRIBUTING.md). Input files for
systems other than the shipped example are especially welcome, and live in
[`contrib/`](contrib/), separate from the code.

Questions go in
[Discussions](https://github.com/damistan01/ARIADNE/discussions); bugs in
[Issues](https://github.com/damistan01/ARIADNE/issues).

---

## Credits

The physics is the **Pathways 1.4** VMD plugin by Ilya A. Balabin, Xiaoqing Hu
and David N. Beratan (Duke University). ARIADNE only orchestrates it; all
coupling calculations happen inside that unmodified plugin. Please cite the
Pathways authors for any scientific use.
