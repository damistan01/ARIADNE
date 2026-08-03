# Tutorial: your first ARIADNE run

A complete worked run, start to finish, on ten frames. It should take a few
minutes. Follow it once with the example as written, then change the donor and
acceptor to your own system.

The example is **flavocytochrome b2**: an electron transferring from the FMN
cofactor to the heme iron. Nothing about the workflow is specific to that
system — only the two selections in [`pathways.in`](pathways.in) are.

You need ARIADNE installed and working; see
[Installation](../README.md#installation) in the main README. You also need
your own topology and trajectory, since neither is distributed here.

---

## 1. Point the input file at your system

Open [`pathways.in`](pathways.in). Every key is commented; the only two lines
you must change are these:

```
donor    = resname FMN and name N5
acceptor = resname FE1 and name FE
```

These are **VMD atom selections**. Each should resolve to a single atom. The
surest way to get them right is to check in VMD before you run anything:

```bash
vmd system.prmtop system_md.nc
```

then in the VMD console:

```tcl
[atomselect top "resname FMN and name N5"] num
```

If that prints `1`, the selection is unambiguous. If it prints `0` your
selection is wrong; if it prints more than `1`, add `and resid N` to pin it to
the copy you mean — otherwise the plugin silently picks whichever gives the
larger coupling, which can even switch between frames.

To find candidate atoms in the first place:

```tcl
lsort -unique [[atomselect top "not water and not protein"] get resname]
```

That lists the non-protein residues — cofactors, metals, ligands — which is
usually where donors and acceptors live.

Note the input file sets `frames = 0:9:1`, i.e. the first ten frames only.
Keep it that way for now.

## 2. Check before you commit

```bash
ariadne -i tutorial/pathways.in \
    -p system.prmtop \
    -y system_md.nc \
    -o tutorial_run/ \
    --dry-run
```

`--dry-run` validates everything — VMD, `pathcore`, the plugin, your
selections, the frame range — and prints the plan without computing anything.
Fix whatever it complains about before going further; the
[troubleshooting table](../README.md#troubleshooting) covers the common
messages.

## 3. Run it

Drop `--dry-run` and add `-n 4`:

```bash
ariadne -i tutorial/pathways.in \
    -p system.prmtop \
    -y system_md.nc \
    -o tutorial_run/ \
    -n 4
```

ARIADNE computes one real frame before starting the other workers, so a bad
selection fails in seconds rather than after four processes have each loaded a
trajectory slice. That probe also prints a measured ETA.

Use `-n 4` even on a bigger machine. `pathcore` is single-threaded, so this is
not CPU oversubscription — it is memory contention, and past four workers it
grows fast enough to cancel out the extra parallelism. See
[Measured performance](../README.md#measured-performance).

## 4. Read the results

```bash
cat tutorial_run/summary.txt
```

For the example system this prints:

```
PATHWAY COUPLING STATISTICS

  <T_DA>:        8.1863155e-06
  sigma(T_DA):   5.747262059449895e-06
  <T_DA^2>:      1.0004678264553251e-10
  coherence:     0.669844243797207

FRAMES

  ok: 10
```

`<T_DA>` is the mean coupling across frames; `<T_DA^2>` is what
non-adiabatic rate theory actually wants. The **coherence parameter** is
`<T_DA>^2 / <T_DA^2>`, between 0 and 1: near 1 means one persistent path
dominates and a single structure represents the system fairly; near 0 means the
coupling is dominated by fluctuation, and no single snapshot is representative.

Ten frames is far too few to trust any of those numbers. It is enough to
confirm the machinery works.

Now the per-frame couplings:

```bash
column -s, -t tutorial_run/couplings.csv
```

```
frame  time_ps  t_da         log10_t_da  n_paths  n_steps  status
0      4100.0   1.95005e-06  -5.7100     1        16       ok
1      4200.0   1.66288e-05  -4.7791     1        12       ok
2      4300.0   9.51076e-06  -5.0218     1        12       ok
3      4400.0   4.58865e-07  -6.3383     1        15       ok
4      4500.0   4.3352e-06   -5.3630     1        15       ok
5      4600.0   1.66912e-05  -4.7775     1        12       ok
6      4700.0   1.01695e-05  -4.9927     1        12       ok
7      4800.0   1.13117e-06  -5.9465     1        13       ok
8      4900.0   1.23852e-05  -4.9071     1        10       ok
9      5000.0   8.60241e-06  -5.0654     1        15       ok
```

(The real file carries full precision in `log10_t_da`; it is truncated here to
keep the table readable.)

One row per frame. `t_da` is the coupling of the best path found in that frame;
`log10_t_da` is there because coupling spans orders of magnitude and is almost
always plotted on a log axis. `status` is `ok` for a frame that computed
cleanly.

Watch how much `t_da` moves between adjacent frames. In these ten frames it
spans a factor of 36 — from 4.6e-07 at frame 3 to 1.7e-05 at frame 5, a hundred
picoseconds apart — and the number of steps in the best path changes from 16 to
10 as the route rearranges. That variation is the entire reason this tool
exists: it is what a single-structure calculation cannot show you, and it is
why the coherence parameter above sits at 0.67 rather than near 1.

And the paths themselves:

```bash
head -6 tutorial_run/paths.csv | column -s, -t
```

```
frame  path_rank  t_da         step  atom_index  resid  resname  atom_name  segid  bond_type
0      0          1.95005e-06  0     7244        927    FMN      N5         X
0      0          1.95005e-06  1     7245        927    FMN      C5A        X      covalent
0      0          1.95005e-06  2     7246        927    FMN      C6         X      covalent
0      0          1.95005e-06  3     7247        927    FMN      C7         X      covalent
0      0          1.95005e-06  4     7248        927    FMN      C7M        X      covalent
```

One row per **step** along a path, in order from the donor. `bond_type` says
how the electron got from the previous atom to this one — `covalent` through a
bond, `hbond` across a hydrogen bond, `through_space` as a jump through
nothing. Step 0 is the donor itself, so its `bond_type` is empty: there is no
preceding atom. Through-space steps are the expensive ones; a path with several is a
poorly coupled path.

`path_rank` is 0 for the dominant path. With `npaths = 1` that is all you get.
Raise `npaths` in the input file to see competing routes.

## 5. Watch it move

```bash
ariadne view tutorial_run/
```

VMD opens with the pathway drawn through the structure, donor and acceptor held
still so that only the route between them moves. Press play, or drag the frame
slider, as with any trajectory.

Colours follow `bond_type`: orange for covalent, lighter orange for hydrogen
bonds, yellow for through-space. Only the protein within 12 Å of the path is
drawn as cartoon, and whatever cartoon still sits between you and the path is
clipped away — without both of those the ribbons wash out the very thing you
are looking at.

If you rotate the view, type `reclip` in the VMD console to re-align the
clipping plane with the new viewing direction.

No graphical display? Use `ariadne view tutorial_run/ --text --export` to write
a self-contained `view/` bundle you can open on a machine that has one.

## 6. Scale up

Once the ten-frame run looks sane, widen the range. Edit `frames` in the input
file, or override it on the command line:

```bash
ariadne -i tutorial/pathways.in \
    -p system.prmtop -y system_md.nc -o tutorial_run/ \
    -n 4 --frames 0:999:1 --resume
```

`--resume` keeps the ten frames you already computed instead of redoing them.
It also protects you if the run is interrupted: re-issue the same command and
it picks up where it stopped.

Before launching something that runs overnight, two things worth deciding:

- **Stride.** `--frames 0:999:2` halves the cost for nearly identical
  statistics, as long as the coupling decorrelates faster than the interval
  between the frames you keep.
- **`npaths`.** Leave it at 1 unless you specifically want to compare competing
  routes. Raising it to 5 costs roughly 3.3x.

Keep `run.json` from any run you intend to report. It records the checksums of
your inputs, the VMD and plugin versions, every resolved parameter and the
frames completed — it is what makes the run reproducible later.
