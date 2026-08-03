# Community contributions

Input files for systems other than the shipped example, contributed by users.
Nothing in this directory is maintained or verified by the ARIADNE maintainer —
each file is the work of whoever submitted it, and works as well as their system
and their care allow.

The tool itself lives in `src/`; nothing here affects how it runs. These are
starting points for people working on similar systems.

## Using one

Point `-i` at it, exactly as with the tutorial input:

```bash
ariadne -i contrib/some_system.in -p your.prmtop -y your.nc -o run/
```

The **selections will not transfer** to your topology unless it is genuinely
the same system with the same residue numbering. Treat the donor and acceptor
lines as something to adapt, not to copy — and check them with `atomselect`
before running, as the [tutorial](../tutorial/README.md#1-point-the-input-file-at-your-system)
describes.

## Adding one

Open a pull request adding a single `.in` file here, named for the system
rather than for you — `azurin_cu_his.in`, not `my_input.in`.

Put a comment header at the top of the file covering:

- **What the system is**, with a PDB code or a citation if there is one.
- **What the donor and acceptor are** chemically, not just as selection strings.
- **Why the selections are written the way they are** — especially if you pinned
  a `resid` to disambiguate between copies of a cofactor, which is the most
  common trap.
- **Anything you had to discover the hard way**: an unusual residue name, a
  numbering offset, a `withh` setting that mattered.
- **Roughly what you got**, so the next person can tell whether their run is
  sane — a mean coupling and the number of frames behind it is plenty.

The commented example in [`tutorial/pathways.in`](../tutorial/pathways.in)
shows the level of detail that is useful.

Please do not commit topologies, trajectories, or any other large binary — the
`.gitignore` blocks the usual extensions deliberately. An input file plus a
note on where the structure came from is what is wanted.
