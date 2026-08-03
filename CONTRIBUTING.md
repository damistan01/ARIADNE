# Contributing

ARIADNE is MIT-licensed. You are free to clone it, change it, and use it for
whatever you like, privately or publicly, without asking anyone. Nothing below
applies to you if that is all you want to do — it only describes how to send
changes back if you would like to.

## Just using it, with local changes

Clone it and edit it. Your changes stay yours; nothing reports back, and there
is no obligation to share anything.

```bash
git clone https://github.com/damistan01/ARIADNE.git
cd ARIADNE
```

If you want your local changes to survive pulling in later updates, commit them
on a branch of your own:

```bash
git checkout -b my-changes
# edit, then
git commit -am "my local changes"
```

Then `git pull --rebase origin clew` brings in upstream work while keeping your
commits on top.

## Sending changes back

Changes arrive as **pull requests**, which appear in the repository's
[Pull requests](https://github.com/damistan01/ARIADNE/pulls) tab, entirely
separate from the main code. Nothing you send is merged unless the maintainer
reviews and accepts it, so there is no way to break anything by trying.

1. **Fork** the repository — the Fork button, top right. You now have your own
   copy that you can push to.
2. **Branch** for your change: `git checkout -b what-it-does`.
3. **Commit and push** to your fork.
4. **Open a pull request** against the `clew` branch. GitHub offers a button
   for this as soon as you have pushed.

Say in the description what the change does and why. If it fixes something
broken, say what you saw go wrong.

Note that the default branch is **`clew`**, not `main` — target that.

## What is especially welcome

**Input files for other systems.** ARIADNE ships one worked example, on
flavocytochrome b2. If you have a working input file for a different
donor–acceptor pair, that is genuinely useful to the next person. These live in
[`contrib/`](contrib/), one file per system, kept separate from the code — see
[`contrib/README.md`](contrib/README.md) for what to include.

**Installation notes for other environments.** The install instructions are
written from one Linux setup. If yours needed different steps — a distribution
where VMD needed extra packages, a cluster module system, a different plugin
layout — that is worth writing down.

**Bug reports.** Open an [issue](https://github.com/damistan01/ARIADNE/issues).
Include the command you ran, what happened, and the contents of `run.json` from
the failed run if there is one — it records versions and every resolved
parameter, which is usually enough to reproduce the problem.

**Questions and discussion** belong in
[Discussions](https://github.com/damistan01/ARIADNE/discussions) rather than
issues.

## A note on scope

ARIADNE deliberately does not reimplement any of the physics — bond perception,
hydrogen-bond detection, pruning and the coupling search all happen inside the
unmodified Pathways 1.4 plugin. Changes that would compute couplings in Python
rather than delegating to the plugin are out of scope, because the guarantee
that a batch frame and a hand-driven GUI frame traverse identical code is the
main thing this tool is for.

Orchestration, output formats, the viewer, documentation, and packaging are all
fair game.

## Code style

Match the surrounding code: standard library only, no runtime dependencies,
type hints on function signatures, and comments that explain *why* rather than
restating the code. Python 3.11+.

The test suite is not distributed with the repository — it hard-codes absolute
paths to a trajectory on the maintainer's machine. If you change behaviour,
describe in the pull request how you checked it.
