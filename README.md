# pp-gen — pseudopotential generation pipeline

Four commands are implemented:

| Command | What it does |
|---|---|
| **`ppgen prepare`** | turns the `input.json` of a run directory into an input file QE `ld1.x` can read directly (such as `U.in`); the scanning workflows generate **a whole batch** at once |
| **`ppgen run`** | runs `ld1.x < U.in > U.out` (equivalent to the reference script `/home/liguozhou/oncv/ppgen/temporary/run.sh`). **Two flows**: without `-i` it runs the existing `U.in` **directly in the current directory** (no generation, no validation); with `-i input.json` it **generates first, then runs** (a scanning workflow generates and runs point by point) |
| **`ppgen plot`** | draws what `plot.json` asks for: radial wavefunctions (`ld1.wfc`) or AE/PS logarithmic derivatives (`ld1.dlog` / `ld1ps.dlog`) |
| **`ppgen report`** | folds one scan into `report.md` (per orbital the best `r_c`, plus the extrema of each channel's poles and zeros at that value and whether they coincide) and draws a "why this `r_c`" figure |

Two workflow shapes:

| Workflow | `prepare` generates | `run` runs |
|---|---|---|
| **`ld1.x`** | one file: `<element>.in` (all-electron, pseudopotential generation, ... all of them; which parameter set is used comes from `parameters.input`) | once |
| `rc_scan` | a batch: `runs/<orbital>/rc_<value>/U.in` (5 orbitals × 11 points = 55 files) | point by point, each point in its own directory |

Six workflows are registered (`pp_gen/inputs/workflows.py`): **`ld1.x`**, **`rc_scan`**,
**`configts`**, **`rlderiv`**, **`rcloc_scan`** and **`energy_scan`**. Each one is
documented in its own section below.

`workflow` only says "run once or scan a batch" — **what exactly runs is decided by the
YAML in `parameters`**. So AE (`iswitch=1`) and the best `r_c` picked after a scan
(`rc_best.yml`, `iswitch=3`) use the **same** `ld1.x` workflow, only with a different
`parameters.input`:

```json
{"workflow": "ld1.x", "parameters": {"input": "ae_radial"}, "qe_path": "..."}   // all-electron
{"workflow": "ld1.x", "parameters": {"input": "rc_best"},   "qe_path": "..."}   // generate a PP at the best r_c
```

> **A note on language.** This README is English, and the package is being translated
> module by module. The `--help` banners, the `prepare` / `run` / `plot` / `report`
> output labels and the generated `report.md` are already English, so the field names in
> the sample blocks below are the ones the tool really prints (verified with `--dry-run`
> against the real run directories). A few module docstrings are still Chinese in the
> source at the time of writing (`diagnostics/selection.py`, `inputs/request.py`,
> `ld1/runner.py`). The test suite is red for the same reason — see "Install and test".

## Usage

```bash
# a single input (all-electron radial wavefunction): run -i generates U.in first, then runs it
ppgen prepare -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/input.json
ppgen run     -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/input.json
ppgen plot    -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/plot.json

# generate only, do not run (to see what would be produced)
ppgen prepare -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/input.json --dry-run

# scanning (rcut per orbital): generate 55 inputs at once and run them one by one
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/rc_scan/input.json
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/rc_scan/input.json --only 5F --rc 1.0   # only one point

# logarithmic derivative reference radius scan (which energy the poles fall on)
ppgen run    -i /home/liguozhou/oncv/uranium-ppgen/best/test_rlderiv/input.json
ppgen report -i /home/liguozhou/oncv/uranium-ppgen/best/test_rlderiv/output.json \
                /home/liguozhou/oncv/uranium-ppgen/best/test_rlderiv/plot_dlog.json

# local potential core radius scan (find the best rcloc)
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/rcloc_scan/input.json
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/rcloc_scan/input.json --rc 2.5   # only one point

# energy scan (find ener for the extra projector with occ=0)
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/energy_scan/input.json
ppgen run -i /home/liguozhou/oncv/uranium-ppgen/energy_scan/input.json --only 5F-2 --ener=-0.6

# you edited one point's U.in by hand and only want to re-run it as is: go there and run
cd /home/liguozhou/oncv/uranium-ppgen/rc_scan/runs/5F/rc_1p2 && ppgen run
```

`prepare`'s output (the products land next to `input.json`, i.e. in the run directory):

```
workflow   : ld1.x
parameters : /home/liguozhou/oncv/pp-gen/src/pp_gen/parameters/input/ae_radial.yml
qe_path    : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin (was a placeholder in input.json,
             filled with the default and written back to input.json)
ld1.x      : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x
generated  : /home/liguozhou/oncv/uranium-ppgen/ae_radial/U.in
```

`run`'s output (really running the local QE 7.4):

```
run dir    : /home/liguozhou/oncv/uranium-ppgen/ae_radial
input      : U.in (266 B, regenerated from the parameters, overwrote an existing file of that name)
qe_path    : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin
ld1.x      : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x
command    : cd /home/liguozhou/oncv/uranium-ppgen/ae_radial && /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x < U.in > U.out
exit code  : 0
seconds    : 0.6 s
stdout     : U.out (5578 B, overwritten)
stderr     : U.err (170 B, overwritten, echoed to the terminal)
products   : ld1.wfc (551902 B, overwritten), vx.pot (181631 B, overwritten)
```

The `input` and `products` lines are deliberately spelled out, because those are the
two things people ask about most:

* whether the input file was **used as is** (no `-i`: "direct run: no input.json read,
  no generation, no validation") or **regenerated from the configuration** (with
  `-i`: "regenerated from the parameters, overwrote an existing file of that name";
  on a dry run it says "dry run: would regenerate and overwrite, the N B file on disk
  is untouched");
* whether each product **overwrote** an existing file or was **created** by this run —
  file names and sizes alone cannot tell those two apart.

### The two halves inside `ld1ps.wfc` (a QE convention)

`ld1.x`'s `ld1ps.wfc` **holds two sets of orbitals in one file** (source
`write_resultsps.f90:181-204`): the first half is the **pseudo** orbitals φ, the second
half is the corresponding **all-electron** orbitals Ψ multiplied by an aligned phase —
so the state labels in the header appear in pairs (for U: `5F 6S 6P 6D 7S | 5F 6S 6P 6D 7S`).
`ppgen plot` honours this convention: the state names it reads are `5F (PS)` / `5F (AE)`,
and when drawing, **the same state shares one colour, solid = AE, dashed = PS** — which
is exactly how Fig. 2 (left panel) of the paper is drawn. As a side effect, duplicate
labels that do not follow the convention are no longer an error; they get a `#2`, `#3`
suffix.

So to draw "pseudo orbital vs all-electron orbital", **listing `ld1ps.wfc` alone in
`data.wfcs` is enough**; adding `ld1.wfc` on top only adds duplicate AE curves (the AE
columns of that file agree with the `(AE)` half — for individual states up to an overall
sign, which is the phase alignment QE performed).

`ppgen run -i input.json` works with relative paths too: paths are resolved to absolute
paths before being displayed (otherwise reports would show `run dir : .` and `cd .`,
which tell you nothing about where the run actually happened).

`plot`'s output:

```
plot.json  : /home/liguozhou/oncv/uranium-ppgen/ae_radial/plot.json
mode       : wavefunction
data       : ld1.wfc (label 'ae_radial', 3427 grid points, r = 1.93e-06–53.08 Å)
             states: 4F, 5D, 5F, 6S, 6P, 6D, 7S
quantity   : rR(r) (inferred from ylabel 'rR(r)')
x unit     : Å (inferred from xlabel 'r(Å)')
axis labels: x='r(Å)', y='rR(r)'
title      : ae_radial
axis ranges: x 0.01–10 Å (log axis), y -2–2
core states: unspecified (all states share one palette; write figure.core to colour by core/valence)
generated  : /home/liguozhou/oncv/uranium-ppgen/ae_radial/ae_radial.png (134863 B)
```

### Switches

| Command | Switch | Effect |
|---|---|---|
| `prepare`, `plot` | `-i/--input PATH` | required; the path of the configuration file (`input.json` / `plot.json`) |
| `run` | `-i/--input PATH` | optional: an `input.json` (render from it first, then run) or a plain `ld1.x` input such as `U.in` (that is a direct run); without it the `U.in` of the current directory is used |
| `report` | `-i/--input JSON [JSON ...]` | required; two files for the scan workflows — `output.json` (run record) and `plot.json` (analysis template) — and a single `output.json` for `configts` |
| all four | `-o/--outdir DIR` | output directory (default: the directory of the configuration file) |
| `prepare` | `-n/--dry-run` | print what would be generated; write no file and do not write back `input.json` |
| `prepare`, `run` | `--only NAME[,NAME]` | handle only these groups (`rc_scan`: an orbital name such as `--only 5F`; `energy_scan`: a card key such as `--only 5F-2`; comma separated) |
| `prepare`, `run` | `--rc VALUE[,VALUE]` | handle only these scan values (`rc_scan`: `rcut` in bohr; `rcloc_scan`: `rcloc` in bohr; comma separated) |
| `prepare`, `run` | `--ener VALUE[,VALUE]` | handle only these `ener` values (Ry, e.g. `--ener=-0.6`; `energy_scan`) |
| `run` | `--qe PATH` | override `qe_path` (QE's `bin` directory or `ld1.x` itself) |
| `run` | `-n/--dry-run` | print the commands that would run; run nothing, write no file |
| `run` | `--timeout SEC` | timeout in seconds (default 3600; `0` means no limit) |
| `plot` | `-n/--dry-run` | read the data and report what would be drawn; write no image |
| `report` | `--minimum` / `--median` / `--average` / `--maximum` | which `\|ΔL\|` statistic ranks the points (mutually exclusive; `--median` is the default) |
| `report` | `--pole` | rank by poles instead: whether the poles coincide (`poles.coincide`) first, then the total pole ΔE (the sum of `pairs[].delta`), and only then the median `\|ΔL\|` |
| `report` | `--scope scanned` / `all` | ranking scope: the channel of the scanned orbital (default) / all channels pooled |
| `report` | `--no-figures` | do not draw per scan point or write `parser.json`, only `report.md` |
| `report` | `--no-decision-figure` | do not draw the "why this `r_c`" figure (`report_<metric>.png`) |
| `report` | `-n/--dry-run` | compute but do not write `report.md` (the body is printed) |
| top level | `--version` | print the version (`ppgen 0.1.0`) |
| top level | `-h/--help` | help for `ppgen` and for every subcommand |

The generated `U.in` **contains no comments at all**, only the namelist itself (the
following corresponds to the current `ae_radial.yml`, i.e. U all-electron with the fine
`dx=0.005` grid):

```
&input
  title       = 'U all-electron',
  atom        = 'U',
  zed         = 92,
  dft         = 'PBE',
  config      = '[Rn] 5f3 6d1 7s2',
  rel         = 1,
  iswitch     = 1,
  max_out_wfc = 7,
  xmin        = -8.0,
  dx          = 0.005,
  rmax        = 100,
/
```

The variables in `U.in` are exactly the ones written in the YAML, no more and no fewer —
so the three grid parameters above are **what you wrote into the YAML**, not something
pp-gen adds. Running `run` after `prepare` makes the product sizes follow the grid:
`dx=0.005` gives 3427 grid points and a 551902 B `ld1.wfc`; `dx=0.008` (the automatic
value of `ld1.x` for Z≥19) gives 2141 points and 344856 B.

## How to write `plot.json`

```json
{
  "mode": "wavefunction",
  "data": {
    "wfcs": ["ld1.wfc"],
    "labels": ["ae_radial"]
  },
  "figure": {
    "xlabel": "r(Å)",
    "ylabel": "rR(r)",
    "title": "ae_radial",
    "xscale": "log",
    "xlim": [0.01, 10.0],
    "ylim": [-2.0, 2.0],
    "fig_name": "ae_radial.png"
  }
}
```

* **`mode`**: the drawing mode; currently `wavefunction` (radial wavefunctions from
  `ld1.wfc`) and `logfunction` (AE/PS logarithmic derivative curves from `ld1.dlog` /
  `ld1ps.dlog`).
* **`data.wfcs`**: which `ld1.wfc` files to draw, **resolved relative to the directory
  holding `plot.json`**. Several files are overlaid on one figure (for comparison).
  **`data.labels`**: one label per entry of `wfcs`, used in the legend; if omitted the
  file name is used (`ld1`).
* **`figure`**: `xlabel` / `ylabel` / `title` / `xlim` / `ylim` / `fig_name` are used
  literally; without `fig_name` the name of `plot.json` is used (`plot.png`), and only
  the file name counts — no directory (use `-o` for the output directory). Optional:
  `xscale` / `yscale` (`linear`/`log`), `figsize` (inches), `dpi` (default 180),
  `quantity`, `x_unit`, `core`, `csv_name`, `parser_name`.
* **Unknown keys are an error** (typing `fig_name` as `figname` should not be silently
  ignored).

### Logarithmic x axis (`"xscale": "log"`)

There is no 0 on a logarithmic axis, so the lower `xlim` bound must be positive. Writing
it wrong is not fatal, but it is **fixed and explained**:

| What you write | What is used | What the report says |
|---|---|---|
| `"xlim": [0.0, 10.0]` | `0.01`–`10` Å | `xlim lower bound 0 is invalid on a log axis (a log axis has no 0 and no negatives), replaced by 0.01 Å` |
| no `xlim` | `0.01`–`10` Å | `the x axis is logarithmic but xlim was not given, so the reference figure's 0.01–10 Å is used; write figure.xlim to change it` |
| `"xlim": [0.001, 5.0]` | as written | — |

(**The y axis is different**: the x axis has a conventional lower bound of 0.01 Å taken
from the reference figure, the y axis has no such convention, so an unusable `ylim` on a
logarithmic y axis is "dropped + taken automatically from the data + a warning" — no
guessed lower bound.)

The "same physical radius" bound is converted through the display unit: with
`x_unit: "bohr"`, 0.01 Å = 0.0189 bohr. The default range is the one of the reference
figure `fig1_ae_radial_wfc_u_logx_lineary.png` (0.01–10 Å) — without it, autoscaling
would spread the whole grid, 1.9e-6 Å to 53 Å, over the axis and squeeze seven orders of
magnitude together.

Also: grid points with `r ≤ 0` cannot be drawn on a log axis and are skipped (a note is
printed when points really were skipped); the y axis stays linear (only the x axis
changes) — write `"yscale": "log"` to make the y axis logarithmic too.

### Which quantity, and which x unit

The explicit `figure.quantity` / `figure.x_unit` win; when they are absent the choice is
inferred from the axis labels, and **the inference is reported together with its source**:

| `ylabel` written as | Quantity drawn |
|---|---|
| `rR(r)`, `u` | `rR(r)` (signed, so nodes are visible) |
| `R(r)` | `R(r) = rR/r` (the origin counts as 0, where R itself diverges) |
| `4πr²R²(r)`, `rho` | radial probability density `(rR)²` |

| `xlabel` written as | x axis unit |
|---|---|
| `r(Å)`, `r (ang)` | Å (`ld1.wfc` stores bohr × 0.529177210903) |
| `r (bohr)`, `r (a.u.)` | bohr (no conversion) |

When nothing is recognised, the defaults apply: draw `rR(r)` and use Å on the x axis,
and the report says that these are the defaults.

### Are the `ld1.wfc` columns R(r) or rR(r)?

**They are `rR(r)`** (that is, `u(r) = r·R(r)`). Checked against the QE 7.4 source:

1. Nothing is rescaled where the file is written: `write_results.f90` does
   `psiaux(:,counter)=psi(:,1,i)` → `write_files.f90`'s `write_wfcfile`, which is
   `WRITE(19,'(38f20.12)') grid%r(n), (wfc(n,ns),…)` (header `# r <labels…>`, 12 decimals
   per row — matching the actual file).
2. U uses `rel=1`, so `scf.f90` goes through `CALL lschps(…, psi(1,1,n), nstop)`, and the
   output argument of `lschps.f90` is literally **called** `u`; its normalization is
   `sn = Σ grid%dx·r(i)·u(i)²` — on a logarithmic grid `dr = r·dx`, i.e. **∫u²dr = 1**,
   which is precisely the normalization of `u = rR`.
3. Measured on this `ld1.wfc`: all seven states give `∫ψ²dr = 1.000000`, while
   `∫ψ²r²dr` ranges from 0.35 to 18.7, and `∫ψ²r dr` / `∫ψ²r²dr` equal the `<r>` / `<r2>`
   that `ld1.x` itself prints (7S: 4.0641 / 18.6879).

`ppgen plot` therefore treats these columns as `rR(r)` directly; to get `R(r)` or the
probability density, write `figure.quantity` (see the table above) — no manual conversion.

### How figures are drawn

* **Colour by state**: if `figure.core` is given (e.g. `["4F", "5D"]`) the reference
  script's "core = blue family / valence = warm colours" scheme is used; otherwise all
  states share one qualitative palette, and the figure **does not claim** which state is
  core and which is valence — that split is element dependent and should not be decided
  by the code.
* **Line style by file**: when several `wfcs` are overlaid, the first file is solid, the
  second dashed, the third dash-dot, the fourth dotted; a single file is all solid. The
  legend then splits into two blocks: `states` (colour) and `files` (line style).
* States are drawn in `n`/`l` order (`7S 6D 5F 6P 4F 5D 6S` → `4F 5D 5F 6S 6P 6D 7S`),
  and the line colours are assigned in that order.
* A linear y axis gets a `y=0` reference line; a logarithmic y axis plots `|quantity|`
  and the axis label becomes `|rR(r)|` automatically (nodes show up as downward spikes).
* A note is added only when an entire curve falls outside `ylim`, or when the grid is
  shorter than the upper `xlim`.
* **Text inside the figure is English**: the default matplotlib fonts contain no CJK
  glyphs, so Chinese text would be drawn as boxes (writing some produces a warning).
  Latin-extended characters such as `Å` are fine.

Drawing needs matplotlib (`pip install -e '.[plot]'`); `cli` imports it lazily, so
`prepare` / `run` work fine without matplotlib.

### The logarithmic derivative figure (`"mode": "logfunction"`)

When scanning `r_c`, what really matters is the **logarithmic derivative curve**: at a
fixed radius `r_d` (the `rlderiv` of `U.in`) one solves the radial equation for
`L_l(E) = u_l'(r_d,E)/u_l(r_d,E)`; `ld1.x` writes the AE side to `ld1.dlog` and the PS
side to `ld1ps.dlog`, one row per energy as `energy  L_{l=0}  L_{l=1} …` (Ry).

```json
{
  "mode": "logfunction",
  "data": {"wfcs": ["ld1.dlog", "ld1ps.dlog"], "labels": ["AE", "PS"]},
  "figure": {"xlabel": "Energy(Ry)", "ylabel": "Log. der. (arbs. units)",
             "yscale": "log", "xlim": [-5, 5], "fig_name": "logfunction.png"}
}
```

* **`data.wfcs`**: any number of `dlog` files (**the first counts as AE, the second as
  PS** — the pole comparison happens between those two). `labels` are the names used in
  the legend.
* **Colour by channel `l`, line style by file** (file 1 solid, file 2 dashed), so the
  AE/PS pair of the same channel is easy to match by eye.
* **A logarithmic y axis plots `|L|`** (the axis label becomes `|L|` automatically): `L`
  diverges at a pole (in this data `|L|` spans 3e-4 to 1e4, seven orders of magnitude),
  so a linear axis would have to clip.
* **A `ylim` that is unusable on a log y axis is dropped with a warning** rather than
  guessed: the sensible lower bound depends on the data, so guessing would only mislead.
  `ylim: [-1.5, 2.0]` (non-positive lower bound) therefore becomes "automatic from the
  data" plus a warning; a valid but narrow `ylim` (two positive numbers) is used as
  given, with a note about how many points were clipped. **Without a `ylim` it is
  automatic**, and the report's `axis ranges` line marks each axis as "given" or
  "automatic from the data".
* **No extrema CSV is written** (that table holds the extrema of `rR`, which has nothing
  to do with the logarithmic derivative): this mode does not write one by default, and a
  `figure.csv_name` is ignored with a note.
* No reference line is drawn at 0 on a linear y axis — for `L`, 0 is an asymptote, not a
  reference value.

### `data.orbitals`: draw only these orbitals

```json
"data": {"wfcs": ["ld1.wfc", "ld1ps.wfc"], "labels": ["AE_radial", "PS_radial"],
         "orbitals": ["6S", "7S", "6P", "5F", "6D"]}
```

Everything is drawn when it is absent. The rules:

* matching is by **state name**, case insensitive: `6S` selects `6S`, `6S (PS)` and
  `6S (AE)` at once (both halves of the paired headers of `ld1ps.wfc` count as the same
  state);
* filtering happens "after reading the files, before states are compared" — so the
  **figure, the extrema CSV and the report's "states:" entry** list only the selected
  orbitals, and the cross-file "states differ" note is no longer disturbed by states that
  are not drawn (filter out `4F`/`5D`, which have no pseudo orbital in `ld1.x`, and that
  note disappears);
* a name found in none of the files → a **warning**, and it is ignored; when nothing
  matches at all → an error (no empty figure is drawn).

### The `analysis` section: pole tolerance and error statistics

```json
"analysis": {
  "pole_threshold": 0.02,
  "mark_poles": true,
  "error": {"median": true, "average": true, "maximum": true}
}
```

**Whether vertical lines appear on the figure is explicitly up to you**: `mark_poles` /
`mark_zeros` default to **not marking** when absent (the same rule as `error` below —
writing the block means following it, missing keys count as `false`; defaults apply only
when the whole block is absent).

* **`mark_poles`** (absent = **do not mark**): marks the poles **on the figure** — one
  vertical line per pole, **coloured by channel** (the colour of the curve) and with
  **line style distinguishing the source** (AE dotted / PS dash-dot); the pair that could
  not be matched (the PS pole **and its AE partner**) is drawn thicker and more solid,
  with **`ΔE`** written vertically at the top (e.g. `Δ=0.180 Ry`) — so "the pair that did
  not line up" shows up as two thick vertical lines. The legend gains the two entries
  `AE: pole / zero` and `PS: pole / zero`. Absent (or `false`) means not a single
  vertical line — **the report still gives everything, it is just not drawn on the
  figure**.

  > The lines were once too faint (`lw=1.0 / alpha=0.35` dotted, almost invisible among
  > eight curves — measured, that column was only 15% darker than its neighbours). Now a
  > matched pole is `lw=1.4 / alpha=0.85` and an unmatched one is `lw=2.0 / alpha=1.0`,
  > measured 1.5–4× darker over the same columns.
* **`mark_zeros`** (absent means **follow `mark_poles`**): also marks the **zeros of
  `L(E)`** (`u'(r_d,E) = 0`, which on a logarithmic y axis are the places where `|L|`
  spikes downwards to 0) — drawn as **short ticks along the bottom**, so **the height
  tells them apart**: a line spanning the figure is a pole, a bottom tick is a zero;
  colour still follows the channel and line style still distinguishes the source (AE
  dotted, PS dash-dot). The two legend entries then read `AE: pole / zero` and
  `PS: pole / zero`. So writing only `mark_poles: true` marks both poles and zeros; write
  `"mark_zeros": false` for poles only, or `"mark_poles": false, "mark_zeros": true` for
  zero ticks only.

  > The criterion is complementary to the pole one, with the same threshold `amp=20`: **a
  > sign change with a very large `|L|` on one side is a pole** (`u=0`), **a sign change
  > with moderate `|L|` on both sides is a zero** (`u'=0`); zero positions are
  > interpolated linearly to `L=0`. The report also gives the number of zeros per channel
  > (`L zeros AE/PS 2/2`).
* **`pole_threshold`** (Ry, default `0.15`, same as the reference script): a PS pole that
  is no further than this from the **nearest** AE pole on the x axis (energy) counts as
  **the same pole**; anything further is "unmatched", and the report lists its **x axis
  distance `ΔE`** to the nearest AE pole one by one (otherwise it shows the pole's x axis
  distance). When the window holds no AE pole at all, the report says "(no AE pole in the
  window)". Poles with `ΔE > 1 Ry` are also flagged as **suspected ghosts** (the
  reference script's criterion), as are "poles the AE side has and the PS side does not" —
  and the other direction is visible too.
* **`error`**: which **`|L_PS − L_AE|`** statistics to show (inside the `xlim` window,
  **per channel**): **`minimum` / `median` / `average` / `maximum`**. **Writing `error`
  decides the set**: keys that are not listed count as `false` — `{"median": true}` shows
  only the median, no mean or maximum appears; all four set to `false` (or an empty
  block) prints not a single number. **Only omitting the whole block** means all four.
  **This is the only set of statistics** — pole position differences get no additional
  mean/median style summary (the per-pair differences are in `parser.json`'s `pairs`).
* Both only **take effect in `mode="logfunction"`** (writing them in another mode gives a
  warning and they are ignored). The statistics go into the command line report and
  **into the figure legend** (the legend uses the English abbreviations `med/avg/max` —
  matplotlib's default fonts cannot draw Chinese and would render boxes).

Besides the axis ranges, the command line report gives the **pole diagnostics** (the
point of this mode):

```
analysis   : pole_threshold=0.02 Ry, poles marked on the figure, error: minimum, median, average, maximum
             (the energy tolerance for "the same pole"; the error is |L_PS − L_AE|)
poles      : the first file counts as AE, the second as PS
           : l=0  AE/PS 2/2, matched 2, unmatched 0 (of which 0 suspected ghosts); L zeros AE/PS 2/2
                  |ΔL| min 3.666e-09, median 0.00041067, average 0.9187, max 263.82
           : l=1  AE/PS 2/2, matched 1, unmatched 1 (of which 0 suspected ghosts)
                  |ΔL| min 5.6659e-05, median 0.0074168, average 10.884, max 1981
             · unmatched pole 4.490->4.310  Δ=0.180
           : l=2  AE/PS 1/1, matched 1, unmatched 0 (of which 0 suspected ghosts)
                  |ΔL| min 1.3406e-07, median 0.0028996, average 3.2466, max 767.33
           : l=3  AE/PS 2/2, matched 1, unmatched 1 (of which 0 suspected ghosts)
                  |ΔL| min 0.00028888, median 0.019203, average 16.197, max 4846.7
             · unmatched pole 4.070->4.390  Δ=0.320
conclusion : no channel has extra poles
```

The criterion is exactly the one of `/home/liguozhou/oncv/ppgen/script/rc_scan_transferability.py` (compared
channel by channel): a pole is a sign change of `L` between adjacent points with `|L| > 20`
on one side; PS poles are split by their distance to AE poles into `matched` (≤0.15 Ry),
`shifted` (≤1 Ry — almost unavoidable in the high energy region, **not** a ghost) and
`extra` (no match within 1 Ry — **suspected ghost**). **The AE side has poles of its
own**, so ghosts are judged by what is *extra*, not by counting.

### `parser.json`: the machine readable analysis result

In logfunction mode `plot` also writes a **`parser.json`** next to the image (rename with
`figure.parser_name`; write `""` to skip it), packaging this analysis in a form that is
easy to **parse uniformly after scanning 55 points**:

```json
{
  "plot_json": ".../plot.json",
  "generated_at": "2026-09-24T15:46:03+08:00",
  "mode": "logfunction",
  "files": [{"label": "AE", "path": ".../ld1.dlog", "n_points": 501,
             "channels": [0,1,2,3], "energy_min": -5.0, "energy_max": 5.0},
            {"label": "PS", "path": ".../ld1ps.dlog", ...}],
  "window": {"energy_min": -5.0, "energy_max": 5.0},
  "settings": {"pole_threshold": 0.02, "mark_poles": true, "mark_zeros": true,
               "error": {"minimum": true, "median": true, "average": true,
                         "maximum": true}},
  "channels": [
    {"l": 1,
     "poles": {"ae": {"n": 2, "energies": [-1.45, 4.49], "min": -1.45, "max": 4.49},
               "ps": {"n": 2, "energies": [-1.45, 4.31], "min": -1.45, "max": 4.31},
               "coincide": false,
               "pairs": [{"ae": -1.45, "ps": -1.45, "delta": 0.0},
                         {"ae": 4.49, "ps": 4.31, "delta": 0.18}]},
     "zeros": {"ae": {"n": 2, "energies": [-1.904947, 1.347777], ...},
               "ps": {...}, "coincide": false, "pairs": [...]},
     "error": {"minimum": 5.67e-05, "median": 0.0074168,
               "average": 10.884, "maximum": 1981.0}}
  ],
  "figure": ".../logfunction.png"
}
```

* **`settings`** records the parameters of `plot.json` that decided the result
  (the window, `pole_threshold`, whether poles/zeros are marked, which error statistics
  are wanted), so this one file says how the numbers were computed.
* **`poles` / `zeros`** have the same structure: the extreme positions of AE and PS
  (`min`/`max`, i.e. the two ends "largest/smallest"), **whether they coincide**
  (`coincide`: same count and every one matches within `pole_threshold`; both sides empty
  also counts as coinciding), and the per-pair differences (`delta` in `pairs`, where
  `ae: null` means the window holds no AE feature to match).
  **No additional mean/median style summary** (removed on request).
* **`error`** is the minimum / median / average / maximum of `|L_PS − L_AE|` for that
  channel.
* A dry run (`-n`) writes nothing, but the report still gives the path; writing
  `parser_name` in a non-`logfunction` mode is ignored with a note.

### The extrema CSV

Besides the PNG, `plot` writes an **extrema table** (every local maximum/minimum of
`rR(r)`), by default named `<fig_name stem>_extrema.csv` — so `ae_radial.png` gives
`ae_radial_extrema.csv`, with the same name and format as the copy in your reference run:

```
state,n_nodes,idx,r_bohr,r_angstrom,rR
4F,0,1,2.8599281895e-01,1.5134088227e-01,1.7585306514e+00
5D,2,3,7.9708928745e-01,4.2180148597e-01,1.2443807668e+00
5D,2,2,2.7340841474e-01,1.4468150235e-01,-6.9094497348e-01
5D,2,1,8.4856933398e-02,4.4904355341e-01,4.8458155693e-01
…
```

* **Columns**: only `state, n_nodes, idx, r_bohr, r_angstrom, rR` (the `class` / `kind` /
  `abs_rR` of the reference file are dropped). Every number uses `%.10e`, as in the
  reference file.
* **Ordering**: by state first (`n`/`l`, the same order as the figure legend and the
  colours: `4F 5D 5F 6S 6P 6D 7S`), then **inside a state by decreasing `|rR|`**.
* **`idx` is still the index within the state in increasing `r`** and is not renumbered
  by the sort — so after sorting, `idx` skips (the `5D` rows above are 3, 2, 1). That is
  deliberate: it shows at a glance which extremum this is.
* The **node and extremum** algorithms follow `/home/liguozhou/oncv/ppgen/script/ae_radial_wfc.py`: nodes are
  sign changes (linearly interpolated to the zero crossing), extrema below `< 1e-8` are
  filtered out as numerical jitter — `ld1.wfc` stores only 12 decimals, and without that
  filter the tails would sprout a crowd of spurious extrema.
* **Renaming / switching off**: `figure.csv_name` (e.g. `"peaks.csv"`) renames it;
  writing `""` writes no CSV at all. With several `wfcs` there is **one table per file**
  (`<stem>_<label>.csv`) — states of the same name from different files cannot be told
  apart in a single table.
* The state order and the colour assignment use the same `n`/`l` ordering, so the CSV row
  order matches the figure legend one to one.

## The two flows of `ppgen run`

```bash
ppgen run -i .../ae_radial/input.json      # (1) render from input.json first, then run
cd .../runs/5F/rc_0p9 && ppgen run         # (2) run the U.in of the current directory as is
ppgen run -i U.in                          # same as (2) (a plain ld1.x input was given)
```

**(1) with `-i input.json`: render first, then run.** Every invocation renders the inputs
from the configuration again (`<element>.in`; a scanning workflow renders one file per
point, `runs/<orbital>/rc_<value>/U.in`), **overwriting same-named files**, and then runs
`ld1.x` on each one. Hand-edited `U.in` files are therefore overwritten — keep them by
using flow (2), or by editing the corresponding parameter in `input.json` so that what is
rendered is what you want.

**(2) without `-i` (or with `-i U.in`): run directly in the current directory.** That is
"run this directory's input once": **no `input.json` is read, no input is rendered and
nothing is validated** (a bad input is reported by `ld1.x` itself). The input file is
found as `U.in` → the only `*.in` in the directory → otherwise an error asking you to name
it; `ld1.x` comes from `PPGEN_QE_PATH` (default
`/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin`) or from `--qe PATH`. `--only` / `--rc` /
`--ener` are rejected in direct mode (those belong to the `input.json` flow);
`--dry-run` / `--timeout` work as usual.

Both modes write `output.json` next to the input (in a direct run it records
`"direct": true` and an empty `"workflow"`), and the report's input line reads
"direct run: no input.json read, nothing rendered, nothing validated".

## `ppgen report`: folding a scan into a report

```bash
ppgen report -i output.json plot.json --median   # pick the r_c with the smallest |ΔL| median
ppgen report -i output.json plot.json --pole     # pick the r_c by pole coincidence + pole ΔE
```

Both inputs are passed through `-i` (order does not matter, **they are identified by
content**: an `output.json` has `tasks`, a `plot.json` has `mode`/`figure`):

* **`output.json`** — the run record written by `ppgen run`: which scan points exist, in
  which directory, for which orbital, with which `r_c`, and whether each one succeeded;
* **`plot.json`** — the analysis template: which two files to compare
  (`ld1.dlog` / `ld1ps.dlog`), the energy window (`xlim`), the pole tolerance and the
  switches (`analysis`).

It **recomputes the same analysis for every scan directory**, then picks the `r_c` with
the smallest ranking statistic per orbital, and writes **`report.md`** next to
`output.json` (`-o` changes the directory, `-n` prints instead of writing), together with
a **`report_<metric>.png`** decision figure (`--no-decision-figure` turns that one off).

Along the way, **every scan point that ran successfully** also produces, inside its own
directory (`--no-figures` turns this off):

| Product | Description |
|---|---|
| `logfunction.png` | exactly what `ppgen plot` draws (pole lines / zero ticks, `ΔE` annotations, legend with `min/med/avg/max`) |
| `parser.json` | the machine readable per-channel result (pole/zero extrema and coincidence, `pairs`, `|ΔL|` statistics) |
| `plot.json` | an analysis template written **only when the directory has none** — an existing one is never overwritten (you may have edited it) |

So a directory such as `runs/5F/rc_0p9` carries its own figure and analysis result right
after `report`, with no need to hand-write a `plot.json` or call `ppgen plot` per point.
A point whose figure fails records one note and does not affect the rest of the report.

### The decision figure: `report_<metric>.png`

The figure next to `report.md` answers **"why this `r_c`"**: one subplot per orbital,
**the x axis is always the scan value of every point**, the y axis is **the quantity that
decided the ranking**, and the selected point is marked with a star and a dashed vertical
line:

| Command | y axis (main) | light grey on the right (comparison) |
|---|---|---|
| `--pole` | **pole ΔE sum** (the sum of `pairs[].delta`), dashed line = `pole_threshold` | median of `\|ΔL\|` (the secondary criterion when ΔE ties) |
| `--median` / `--average` / `--maximum` / `--minimum` | **that statistic of `\|ΔL\|`** (log axis) | pole ΔE (dashed line = `pole_threshold`) |

A few things to read from the figure:

* **filled marker = poles coincide, open marker = they differ** — so with `--pole` you
  see at a glance which points are inside the tolerance (`6S` at 1.2/1.3/1.7/1.8 are
  filled, `5F` is open everywhere);
* **`∞`**: that channel has no AE pole in the window to match, and it is noted at the top
  of the subplot;
* failed points get no curve; their `r_c` is recorded in the corner of the subplot
  (`failed r_c: 0.6, 0.7, 1`);
* on the `--median` figure you can see directly that the point with the smallest median
  already has its poles outside the tolerance line — which is exactly why `--pole` gives a
  different answer.

(These products are **not written into `report.md`** — look at the directory when you need
them; the report only lists "points that failed to plot".)

```
run record : /home/liguozhou/oncv/uranium-ppgen/rc_scan/output.json
template   : /home/liguozhou/oncv/uranium-ppgen/rc_scan/plot.json
window     : -5 ... 5 Ry
pole tol   : 0.15 Ry
ranking    : median of |ΔL| (scope: scanned channels)
best rcut:
  5F   l=3  r_c = 1.3 bohr   |ΔL| median = 0.018965   rc_1p3
  6S   l=0  r_c = 1.2 bohr   |ΔL| median = 4.8892e-05   rc_1p2
  6P   l=1  r_c = 1.8 bohr   |ΔL| median = 0.0063025   rc_1p8
  6D   l=2  r_c = 2.4 bohr   |ΔL| median = 0.0028996   rc_2p4
  7S   l=0  r_c = 1.4 bohr   |ΔL| median = 6.2735e-05   rc_1p4
usable pts : 45/55
per point  : 45 directories each get logfunction.png and parser.json (e.g. .../runs/5F/rc_1p3)
decision   : .../report_median.png (x: rcut, y: median of |ΔL|; the light grey on the right is the pole ΔE)
generated  : .../report.md
```

* **The ranking statistic**: `--median` (default) / `--average` / `--maximum` /
  `--minimum` — the point whose chosen statistic of `|L_PS − L_AE|` is smallest wins.
* **`--pole`: pick `r_c` by "do the poles line up"** (never by which `|ΔL|` is smaller).
  The order is:

  1. **whether the poles coincide** (`poles.coincide` in `parser.json`) — coinciding
     points come first;
  2. when both coincide (or both differ), the **pole ΔE sum** — the sum of
     `poles.pairs[].delta` (e.g. 0.32 for `l=3` of 5F, 0.14 for `l=1` of 6P), **smaller
     is better**; pairs for which the window holds no AE pole at all (`delta: null`) have
     no ΔE to compare, count as "structurally inconsistent", are recorded as **∞** and
     ranked last within the same group;
  3. only when the ΔE also ties, the **median of `|ΔL|`** (the secondary criterion, always
     the median).

  In other words: "poles lining up wins first, then which one lines up most accurately, and
  only if all that ties does the curve difference matter". `--scope` applies here too:
  `scanned` looks only at the channel of the scanned orbital, `all` adds up the ΔE of
  every channel. The report lists the pole ΔE **per channel** (bold = not coinciding), so
  it is obvious which channel fails to match:

  ```
  ranking    : pole coincidence -> pole ΔE total -> median of |ΔL| (scope: scanned channels)
  best rcut:
    5F   l=3  r_c = 1.2 bohr   poles differ  ΔE total 0.32  |ΔL| median = 0.019056   rc_1p2
    6S   l=0  r_c = 1.2 bohr   poles coincide  ΔE total 0  |ΔL| median = 4.8892e-05   rc_1p2
    6P   l=1  r_c = 1.8 bohr   poles differ  ΔE total 0.14  |ΔL| median = 0.0063025   rc_1p8
    6D   l=2  r_c = 2.4 bohr   poles coincide  ΔE total 0.02  |ΔL| median = 0.0028996   rc_2p4
    7S   l=0  r_c = 1.4 bohr   poles coincide  ΔE total 0  |ΔL| median = 6.2735e-05   rc_1p4
  ```

  Measured, this criterion **changes the answer**: `--median` picks `r_c = 1.3` for 5F
  (median 0.018965, the smallest), while `--pole` picks `r_c = 1.2` — because the `l=3`
  pole ΔE is 0.34 Ry at 1.3 and 0.32 Ry at 1.2; only among the points with the same ΔE
  (0.9 / 1.1 / 1.2) does the median decide.
* **The ranking scope `--scope`** (this one matters a lot): only one orbital is scanned at
  a time, so **the other** channels stay the same throughout a scan and pooling them in
  dilutes the difference — measured, the pooled median of the 5F run differs only in the
  fourth significant digit between points (0.0043246 vs 0.0043247), while the median of
  the **scanned channel** (`l=3`) moves from 0.0210 to 0.0190, which does separate good
  from bad. Hence the default `--scope scanned`: rank by the **channel of the scanned
  orbital** (the trailing letter of the orbital label → `l`: s=0 p=1 d=2 f=3 …);
  `--scope all` pools every channel.
* Besides the per-point detail table (with the statistics of every channel and the failure
  reasons), `report.md` has a **"Poles and zeros at the best parameter"** section: for the
  `r_c` selected per orbital it lists, channel by channel, the **extreme positions** of the
  poles (`L` diverges) and of the zeros (`L=0`) — `min`/`max` plus every position — and
  **whether they coincide** (same count and every one matching within `pole_threshold`),
  with the per-pair differences. Skipped directories (failed runs, missing
  `ld1.dlog`/`ld1ps.dlog`) get their own section with the reason.

  Pole positions are located **by interpolation** on the energy grid (`deld` of `&input`),
  and two poles count as the same when they agree within `pole tolerance + grid step`
  (0.02 + 0.02 = 0.04 Ry): a pole falls between two grid points, so interpolation is the
  only option and its resolution is one cell, which is why the criterion keeps one cell of
  slack. `ΔE` reports the true interpolated difference, not a multiple of the grid step.

## `configts`: test one pseudopotential against configuration after configuration (transferability, Tables 3 and 5 of the paper)

```json
{"workflow": "configts",
 "parameters": {"input": "rc_scan", "configts": "configts"},
 "qe_path": "/…/bin"}
```

```yaml
# parameters/configts.yml -- holds nothing but "which configurations to test"
configts:
  - '5f3 6s2 6p6 6d1 7s2'      # reference configuration (the one used for generation)
  - '5f3 6s2 6p6 6d0 7s2'
  - '5f3 6s2 6p6 6d2 7s1'
```

```yaml
# parameters/input/test_configts.yml -- the numbers live here (the only source of truth)
input:
  iswitch: 2                   # test an existing PP (not generate one)
  # …config/dft/rel/nld/rlderiv/eminld/emaxld/deld as usual
test:
  nconf: 1
  'configts(1)': '5f3 6s2 6p6 6d1 7s2'
  file_pseudo: '/abs/path/U.pbe-n-nc.UPF'   # the file under test (fixed)
```

**Two modes** (decided by the parameter file):

1. **A single run (recommended, simplest)**: the parameter file writes `nconf = N` and
   `configts(1..N)` directly → `ppgen run -i input.json` runs **`ld1.x` once** (about 5 s
   for 5 configurations), and QE prints `dEtot_ae` / `dEtot_ps` / `Delta E` per
   configuration (writing `ld1<n>*.dlog`). The auxiliary group `parameters.configts` may
   be omitted (and is ignored when given).
2. **Scan one configuration at a time**: the parameter file lists only one `configts(1)`
   (`nconf = 1`) and the list moves into the auxiliary group `parameters.configts` → one
   task per configuration, landing in `runs/<configuration>/`, where each directory's
   `U.in` **changes only `configts(1)` of `&test`**.

In both modes `iswitch`, `file_pseudo`, `nconf` and the presence of `&inputp` or cards
follow the parameter file (the workflow does not second-guess you). Both modes give exactly
the same δ (measured difference <1 µRy).
* **The parameter file must use `iswitch=2`** (`iswitch=2` together with a card is
  rejected): `iswitch=3` means "generate the PP + test it", and the in-memory PP is
  unscreened with the generation configuration, so as soon as the test configuration
  changes the PS total energy is off (measured on the same configuration: a separate run
  per configuration with `iswitch=3` gives −183.196 Ry, while a single `nconf>1` run or
  `iswitch=2` with a fixed `file_pseudo` gives −181.60 Ry — 1.6 Ry apart) ⇒ ΔE_PS and δ
  become meaningless. Writing `iswitch=3` makes `prepare` / `run` **warn** about this
  (without changing any data).
* `ppgen report -i output.json`: **one input only** (no `plot.json` needed) and it produces
  `report.md`: a transferability table (one row per configuration with
  `E_AE / E_PS / ΔE_AE / ΔE_PS / δ`, where δ = ΔE_AE − ΔE_PS, using the same sign
  convention as QE's `dEtot`) plus a per-configuration AE/PS level table. Its warnings
  come automatically: configurations whose valence electron count differs from the
  reference, an absurdly large `|ΔE_AE|` (the AE side fell onto another solution), and a
  reference configuration missing from the task list (in which case the first successful
  configuration becomes the reference).
* `--only` selects a few by index, directory tag or configuration string; `--rc` is
  meaningless in this workflow (rejected with an error).

## `rc_scan`: scan `rcut` per orbital

The `input.json` uses **two** parameter groups: one for the `ld1.x` namelists, one for the
scan range:

```json
{
  "workflow": "rc_scan",
  "parameters": {"input": "rc_scan", "rc_scan_range": "rc_scan_range"},
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
}
```

* `parameters/input/rc_scan.yml`: `&input` / `&inputp` / the **channel cards** / `&test`,
  where `orbitals` holds the 8-column cards per orbital
  (`label n l occ ener rcut rcutus jj`) and `n_orbitals` the number of cards (a mismatch is
  an error).
* `parameters/rc_scan_range.yml`: the offsets, plus the optional `update` (see below).

```yaml
rc_range: [-0.5, 0.5, 0.1]     # [start, stop, step], both ends included
update: true                   # after an orbital is scanned, set its rcut to the best value
```

**Scan rules** (one orbital moves at a time):

* **The scan centre is the `rcut` of that orbital in the `orbitals` card** (5F → 1.0641
  bohr), with no rounding whatsoever — the `input` parameter is the single source of
  truth, so no separate "centre value" has to be written.
* For every offset of orbital X, `card rcut + offset` (5F → 0.5641, 0.6641 … 1.5641)
  produces one input under `runs/<orbital>/rc_<value>/U.in`, e.g. `runs/5F/rc_1p0641/U.in`
  (**absolute value**, not the offset; the decimal places are written out so the name is
  the true `rcut` of that point).
* That orbital's `rcut` **and** `rcutus` are set to the same value (``ld1.x``'s
  ``read_psconfig`` requires `rcut <= rcutus`); **all other orbitals keep** the exact value
  from the `input` parameter (6S stays at 1.3527), so only one variable moves at a time.
* One directory per point is mandatory: every point's `file_pseudopw` is called
  `U.pbe-n-nc.UPF`, so they would overwrite each other in a single directory.

### Updating `rcut` after the scan (`update`)

With **`update: true`** in `rc_scan_range.yml`, `ppgen run -i input.json` becomes
**serial, orbital by orbital**:

1. render and run every scan point of this orbital with the current parameters;
2. use the resulting `ld1.dlog` / `ld1ps.dlog` to pick the best point — the criterion is
   **exactly** that of `ppgen report --pole`: **do the poles coincide → pole ΔE sum →
   median `|ΔL|`** (looking only at the channel of this orbital, i.e. `--scope scanned`);
3. write that `rcut` **in place** into the parameter file that supplies the `orbitals`
   cards (`parameters/input/rc_scan.yml`), then render and scan the **next** orbital with
   the new value.

So only one variable moves at a time, and **an orbital that has been optimised stays at its
best value from then on**: after 6S is scanned down to 1.2, the inputs of the later
7S/6P/6D/5F runs use 1.2 for 6S. After one pass, the `rcut`s in that parameter file are all
at their current optimum — running `rc_scan` again re-scans finely around them (iterating
to convergence).

The write-back touches **only the two numbers on that line** (`rcut`/`rcutus`);
indentation, quotes, trailing comments, the other columns and the rest of the file are kept
verbatim (these YAML files are hand written, and reflowing the whole file would destroy the
comment layout).

```
update     : pole coincidence -> pole ΔE sum -> median |ΔL| (pole tolerance 0.15 Ry);
             energy window from eminld/emaxld of &input; groups run serially: a group's
             scan values are replaced by its best values (written back to the parameter
             file) before the next group is scanned
…
update rcut: 6S -> 1.2 bohr (was 1.3)
           : best point .../runs/6S/rc_1p2, poles coincide, ΔE total 0 Ry, |ΔL| median 4.8892e-05
           : written back to /home/liguozhou/oncv/uranium-ppgen/rc_scan/parameters/input/rc_scan.yml (later groups use it)
```

* **The criterion can be tuned** (absent means the `true` variant):

  ```yaml
  update: true                                    # same as {"criterion": "pole"}
  update: {"criterion": "pole", "pole_threshold": 0.02}
  update: {"criterion": "median"}                 # rank by the |ΔL| median only
  update: {"criterion": "pole", "window": [-5, 5]}  # default: eminld/emaxld of &input
  ```

  * `criterion`: `pole` (default) / `minimum` / `median` / `average` / `maximum`;
  * `pole_threshold` (default `0.15` Ry, the same default as in `report`) decides what
    counts as "coinciding";
  * `window` defaults to `eminld` / `emaxld` of `&input` (the same window as `xlim` of
    `plot.json`).
* **`ld1.dlog` / `ld1ps.dlog` are required** (i.e. `iswitch=3` plus `rlderiv`, which this
  parameter set has). Points missing those two files (failed runs, or filtered out by
  `--rc`) do not take part in the comparison; **if not a single point can be picked, no
  update happens** and the report and the warnings say why — no value is guessed and
  written back.
* The update record goes into **`updates` of `output.json`** (new value / previous value /
  criterion / whether the poles coincide / ΔE / best point directory / file written back),
  so you can always look up "which scan decided this `rcut`".
* A dry run (`-n`) neither runs nor updates, and the report says so with a warning; the
  direct run (no `-i`) is unrelated to all of this.

`prepare`'s report:

```
scan       : 5 orbitals × 11 points = 55 scan points (centre from the rcut of the
             orbitals card + rc_range offsets…)
           : 5F  rc = 0.5641, 0.6641, 0.7641, 0.8641, 0.9641, 1.0641, 1.1641, 1.2641, 1.3641, 1.4641, 1.5641   ->  runs/5F/rc_*/
output dir : /home/liguozhou/oncv/uranium-ppgen/rc_scan/runs
generated  : 55 U.in files (e.g. .../rc_scan/runs/6P/rc_1p6133/U.in)
```

`run` goes point by point (`ld1.x` is single threaded, so this is sequential), one line per
point plus a summary at the end; **a failure in the middle does not stop anything** — the
point of a scan is to finish the points and then choose:

```
run dir    : /home/liguozhou/oncv/uranium-ppgen/rc_scan
targets    : 55 scan points
input      : rendered 55 U.in files from the parameters (one per scan point, overwriting same-named files)
qe_path    : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin
ld1.x      : /home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x
stderr     : each point writes its own U.err (not echoed point by point)
per point  :
  [  1/55] 6S rc=0.8527  ok  1.6 s  U.pbe-n-nc.UPF, ld1.dlog, ld1.pwe, ld1.test, ld1.wfc, ld1ps.dlog, ...
  [  2/55] 6S rc=0.9527  ok  1.6 s  U.pbe-n-nc.UPF, ld1.dlog, ...
  …
summary    : succeeded 55/55, total time 84.2 s
record     : .../rc_scan/output.json (exit code / exit reason / run directory of every task)
```

Exit codes: 0 when everything succeeded; with **a single target** the exit code of `ld1.x`
is passed through unchanged; with several targets it is 1 when anything failed (the failing
targets are listed on the `failed` line, and the actual error lines are printed as well). A
timeout **counts as a failure of that one point** and is recorded; it does not abort the
batch.

After the run, an **`output.json`** is written next to `input.json` (a dry run writes
nothing), recording **exit code / exit reason / run directory** per task, plus the input
file name, the elapsed time and the artifacts needed to reproduce it:

```json
{
  "input": "/home/liguozhou/oncv/uranium-ppgen/rc_scan/input.json",
  "direct": false,
  "workflow": "rc_scan",
  "run_dir": "/home/liguozhou/oncv/uranium-ppgen/rc_scan",
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin",
  "ld1": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin/ld1.x",
  "started_at": "2026-09-24T14:29:56+08:00",
  "finished_at": "2026-09-24T14:30:00+08:00",
  "summary": {"total": 3, "succeeded": 2, "failed": 1, "seconds": 3.975, "exit_code": 1},
  "tasks": [
    {
      "label": "5F rc=0.6", "orbital": "5F", "scanned": "rcut", "rcut": 0.6,
      "config": null,
      "directory": ".../runs/5F/rc_0p6",
      "input": "U.in",
      "exit_code": 1,
      "exit_reason": "Error in routine run_pseudo (1): / Errors in PS-KS equation",
      "seconds": 0.858,
      "artifacts": ["CRASH", "input_tmp.in", "ld1.dlog", "ld1.wfc", "vx.pot"]
    },
    { "label": "5F rc=0.9", "...": "...", "exit_code": 0, "exit_reason": "finished normally" }
  ],
  "updates": [
    {
      "orbital": "6S", "scanned": "rcut", "rcut": 1.2, "previous": 1.3,
      "criterion": "pole", "pole_threshold": 0.15,
      "coincide": true, "delta": 0.0, "stat": 4.8892e-05,
      "file": ".../parameters/input/rc_scan.yml",
      "written": true, "directory": ".../runs/6S/rc_1p2", "note": ""
    }
  ]
}
```

`updates` has content only with `update: true`: the best value picked after each group was
scanned, the value before and after, the criterion with its numbers, which parameter file
was written, and the best point's directory (see "Updating `rcut` after the scan").

The possible values of `exit_reason`: `finished normally` (exit code 0), QE's
`Error in routine …` lines (non-zero exit code), `still running after N s, terminated`
(timeout, with `exit_code` null), and `input file U.in not found` (that point had no input
and was skipped).

**Some points failing is what a scan is meant to discover.** Scanning a small `rcut` makes
QE's `run_pseudo` report

```
Error in routine run_pseudo (1):
Errors in PS-KS equation
```

(that point has no UPF; the AE half did finish). This is not a tool failure: `ppgen run`
still runs every point, prints the failure list and the error lines, and exits non-zero.
To re-run the failing points, select them with `--only` / `--rc` (measured, the same input
fails reproducibly — it is not random).

In batch mode each point's stderr goes only into its own `U.err` and is **not echoed** to
the terminal point by point (55 points would mean 55 blocks of OpenMPI / `MPI_ABORT`
noise); only a single target is echoed as before.

A scan has filter switches (accepted by `prepare` and by `run`, comma separated and
combinable):

```bash
ppgen run -i .../rc_scan/input.json --only 5F                    # only the 11 points of 5F
ppgen run -i .../rc_scan/input.json --rc 1.0641                  # only the point with rcut=1.0641
ppgen prepare -i .../rc_scan/input.json --only 5F --rc 1.0641    # render only that one point
```

`--rc` takes the **true `rcut`** (card centre + offset), not the offset; when unsure which
values exist, run `--dry-run` first — the report lists every `rcut` each orbital scans.

## `rlderiv`: scan the logarithmic derivative reference radius of `&input`

`rlderiv` decides "at which energy the poles appear" (a pole is the wavefunction crossing
zero at `r=rlderiv`), so changing it can move the poles out of the energy range you care
about — and **it does not change the pseudopotential** (see below). It is not a column of a
card but a scalar of `&input`, using the same namelist scanning machinery as `rcloc`:

```json
{
  "workflow": "rlderiv",
  "parameters": {"input": "energy_scan", "rlderiv": "rlderiv"},
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
}
```

* `parameters.input`: **the PP setup under test** (`&input` / `&inputp` / channel cards /
  `&test`; sections unrelated to this workflow stay as they are), whose `rlderiv` is the
  scan centre;
* `parameters.rlderiv`: the range, plus the optional `update`:

```json
{"rlderiv": [4.0, 8.0, 1.0]}
```

**How a range is interpreted** (this also applies to `rcloc`):

| Form | Meaning |
|---|---|
| `[4.0, 8.0, 1.0]` | both ends positive **and** the centre (the current `rlderiv` of `&input`) falls inside → **absolute** values: 4, 5, 6, 7, 8 |
| `[-1.0, 3.0, 1.0]` | the centre is not inside (or a value is non-positive) → **offsets from the centre** (4, 5, 6, 7, 8 when the centre is 5.0) |
| `{"offset": [1.0, 3.0, 1.0]}` | explicitly offsets |
| `{"absolute": [6.0, 8.0, 1.0]}` | explicitly absolute values |

One directory per value, `runs/rlderiv/rlderiv_<value>/U.in`, where **only `rlderiv` of
`&input`** is changed to that value; `ppgen run -i input.json` runs them one by one:

```
scan       : rlderiv of &input × 5 points = 5 scan points (rlderiv of &input, range from
             ['rlderiv'], one directory per value)
           : rlderiv  = 4, 5, 6, 7, 8   ->  runs/rlderiv/rlderiv_*/
```

Afterwards, in the run directory, `ppgen report -i output.json plot_dlog.json`: **every
point directory gets its own `logfunction.png` and `parser.json`** (using the window,
tolerance and line styles of your `plot_dlog.json`), and `report.md` is written — one row
per `rlderiv`, **one column per `l` giving the median `|ΔL|`**, and at the end the
pole/zero comparison for the selected `rlderiv`:

```
| rlderiv (a.u.) | status | l=0 median | l=1 median | l=2 median | l=3 median | pooled median |
| 4.0 | ok | 0.00021245 | 8.8633e-05 | 4.5817e-05 | 1.7052e-05 | 9.035e-05 |
| 5.0 | ok | 0.00011071 | 2.6424e-05 | 1.6769e-05 | 7.3193e-06 | 2.9083e-05 |
| 6.0 | ok | 8.6703e-05 | 5.6425e-05 | 6.108e-06 | 1.6319e-06 | 1.7644e-05 |
| 7.0 **←** | ok | 5.2452e-05 | 1.8828e-05 | 8.26e-06 | 9.9612e-06 | 1.3269e-05 |
| 8.0 | ok | 5.1653e-05 | 2.5881e-05 | 2.63e-05 | 2.3941e-05 | 2.7342e-05 |
```

Because `rlderiv` affects the whole pseudopotential (every channel), the ranking and the
report automatically pool **all channels** (equivalent to `--scope all`) and say so.

With `update: true` the best `rlderiv` is written back into `&input` of the file given as
`parameters.input` (only the number on that line; comments and other keys are untouched).

> Unlike `rcut` / `ener` / `rcloc`, **`rlderiv` does not change the pseudopotential itself**
> (with the same cards, `rlderiv=5` and `rlderiv=6` produce byte-identical UPF data blocks,
> and the same `ld1.test`); it only changes *where* the logarithmic derivative is looked
> at. So picking a `rlderiv` is a matter of picking a **nicer viewpoint** from "which
> parameter set is most self-consistent", not of changing the physics of the PP.

## `rcloc_scan`: scan `rcloc` of `&inputp` (the core radius of the local potential)

`rcloc` is the core radius of the **local potential** (with `lloc=-2` the local potential is
the one of the `l=2` channel): it belongs to no card but is a scalar of `&inputp`, so this
workflow scans a different axis:

```json
{
  "workflow": "rcloc_scan",
  "parameters": {"input": "rcloc_scan", "rcloc_scan_range": "rcloc_scan_range"},
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
}
```

```yaml
rcloc_range: [-0.5, 0.5, 0.25]   # [start, stop, step], both ends included
update: {"criterion": "pole", "pole_threshold": 0.02}
```

**Scan rules**:

* **The scan centre is the current `rcloc` of `&inputp`** (missing → an error: rather than
  guessing the default of `ld1.x`, you are asked to write a centre explicitly);
  `rcloc_range` gives **offsets from it**, the same convention as `rc_range` /
  `energy_range`;
* every point lands in `runs/rcloc/rcloc_<value>/U.in`, e.g. `runs/rcloc/rcloc_2p5/U.in`;
* only `rcloc` of `&inputp` changes per point; `lloc`, `rclow` (if written) and every other
  section stay untouched;
* it affects **the whole pseudopotential** (all channels) — there is no "scanned orbital /
  card" — so both the criterion and the report use **all channels** (all poles coincide →
  the ΔE sum of every pole → the pooled median `|ΔL|` of every channel), and `ppgen report`
  treats `--scope` as `all` automatically and says so.

```
scan       : rcloc of &inputp × 5 points = 5 scan points (centre from the current rcloc of
             &inputp + rcloc_range offsets, one directory per value)
           : rcloc  = 2, 2.25, 2.5, 2.75, 3   ->  runs/rcloc/rcloc_*/
```

### Updating `rcloc` after the scan (`update`)

The same criterion as the other two scan workflows (the serial iteration simply has a single
group here: `rcloc` itself):

```
update rcloc: rcloc -> 2.2 bohr (was 2.5)
           : best point .../runs/rcloc/rcloc_2p2, poles coincide, ΔE total 0 Ry, |ΔL| median 4.8892e-05
           : written back to .../parameters/input/rcloc_scan.yml (later groups use it)
```

The write-back replaces only **that number** on the `rcloc:` line of `&inputp`: indentation,
trailing comments, the other keys of the section and the rest of the file stay verbatim;
when the key is missing it is inserted at the first line of that section.

Filtering: `--rc` means "which `rcloc` values" in this kind of workflow (`--rc 2.5` runs
that single point); `--only` is useless here (there is only one group).

## `energy_scan`: scan `ener` for the extra projection operator with `occ=0`

Multiple projectors (two cards for the same `l`) raise the question of which reference
energy the extra projector should use: the first card carries the valence occupancy, with
`ener` written as `0.00` (`ld1.x` takes the one-electron energy), and **the second card must
have `occ = 0`** with an `ener` you supply — that `ener` is what `energy_scan` looks for.

```json
{
  "workflow": "energy_scan",
  "parameters": {"input": "energy_scan", "energy_scan_range": "energy_scan_range"},
  "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
}
```

* `parameters/input/energy_scan.yml`: the same structure as `rc_scan` (`&input` / `&inputp`
  / channel cards / `&test`), with at least one card having `occ == 0` and `ener != 0`;
* `parameters/energy_scan_range.yml`: the offsets, plus the optional `update`:

```yaml
energy_range: [-0.2, 0.2, 0.05]   # [start, stop, step], both ends included
update: true                      # after a card is scanned, set its ener to the best value
```

**Scan rules** (one card moves at a time):

* the scan covers **every card with `occupancy == 0`** (`occ` is column 4 of a card) —
  occupied cards are "the state being described", they are not scanned;
* **the scan centre is the current `ener` of that card**, and `energy_range` gives
  **offsets from it** (the same convention as `rc_range`): with `-0.60` on the 5F card and
  the range `[-0.2, 0.2, 0.05]`, the scan walks `-0.80, -0.75 … -0.40`;
* every point lands in `runs/<card key>/ener_<value>/U.in`, e.g. `runs/5F-2/ener_m0p6/U.in`.
  **Card key**: the label when it is unique among the cards (`5F`), and otherwise with an
  index — the parameter set above has two `5F` cards, so the scanned one is called `5F-2`;
* only **that one card's** `ener` changes; the other `occ=0` card and every occupied card
  (`ener`, `rcut`) keep the exact values from the parameters — one variable at a time;
* an `occ=0` card with `ener` written as `0.0` is skipped with a warning (there `ld1.x`
  uses the one-electron energy, so scanning it is pointless).

```
scan       : 2 occ=0 cards × 9 points = 18 scan points (centre from the current ener of
             each occ=0 card + energy_range offsets, one directory per card)
           : 6D-2  ener = -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15   ->  runs/6D-2/ener_*/
           : 5F-2  ener = -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4   ->  runs/5F-2/ener_*/
output dir : .../energy_scan/runs
```

### Updating `ener` after the scan (`update`)

Exactly the same machinery as `update` in `rc_scan` (criterion, tolerance, window and the
serial iteration per group are all the same); the only difference is which column is
changed and which card is written back:

1. render and run all points of **this card**;
2. pick the best point by **poles coincide → pole ΔE sum → median `|ΔL|`** (the channel
   comes from the `l` written on the card, not from guessing it off the label, so a key
   like `5F-2` is still judged correctly);
3. write that `ener` in place into **that card** (located by card index, so several cards
   with the same label cannot be mixed up), then render and scan the **next** card.

```yaml
update: true                                       # same as {"criterion": "pole"}
update: {"criterion": "pole", "pole_threshold": 0.02, "window": [-3, 3]}
```

```
update ener: 5F-2 -> -0.60 Ry (was -0.55)
           : best point .../runs/5F-2/ener_m0p6, poles coincide, ΔE total 0 Ry, |ΔL| median 4.8892e-05
           : written back to .../parameters/input/energy_scan.yml (later groups use it)
```

> **Do not open the window too wide.** Without `window` the `eminld`/`emaxld` of `&input` is
> used; when that is `-5…10` (15 Ry), a crowd of poles in the high energy region that have
> nothing to do with PP transferability defeats the "poles coincide" level, the ranking
> falls back to the ΔE sum, and the point picked can be meaningless (measured: 6P was once
> picked at a point with `ΔE 3.88 Ry / |ΔL| 1.42`). So when `update` is on and the window
> is wider than 6 Ry, `run` warns and suggests writing `window: [-3, 3]`.

Every task in `output.json` records **`"scanned": "ener"`** and **`"ener": <value>`**
(`rc_scan` records `"scanned": "rcut"` and `"rcut"`; older records without `scanned` are
read as `rcut`), and the `updates` entries use the same key names:

```json
{"label": "5F-2 ener=-0.6", "orbital": "5F-2", "scanned": "ener", "ener": -0.6,
 "directory": ".../runs/5F-2/ener_m0p6", "exit_code": 0, "exit_reason": "finished normally"}
```

The header, the units and the x axis of the decision figure in `ppgen report` follow that
axis as well (`ener (Ry)` / "why this ener"), never `r_c (bohr)`.

Filter switches: `--only` takes **card keys** (`--only 5F-2`, or a bare label to select all
scannable cards with that label at once) and `--ener` takes **true `ener` values**
(`--ener=-0.6`):

```bash
ppgen run -i .../energy_scan/input.json --only 5F-2            # scan only that card
ppgen run -i .../energy_scan/input.json --only 5F-2 --ener=-0.6  # run only that one point
```

## A few conventions of `run`

* **Giving `-i` re-renders the input.** That flow means "render the files from
  `input.json` and run them", and it overwrites same-named `.in` files every time; to run a
  hand-edited file as it is, go to its directory and use `ppgen run` (or
  `ppgen run -i U.in`) — the two flows each do one job.
* **The working directory is the run directory.** The products of `ld1.x` (`ld1.wfc`,
  `vx.pot`, `CRASH`, …) land next to `input.json` and do not pollute the directory you are
  in.
* **stdout → `U.out`, stderr → `U.err`** (and echoed to the terminal). The OpenMPI
  `librdmacm.so.1` noise still shows up in the terminal, but it is also kept in `U.err`.
* **Failures must be visible.** With a non-zero exit code, the `Error in routine …` lines
  of `U.out` and a few trailing lines are put right in front of you, together with the list
  of files that appeared in this run; the exit code is returned to the shell unchanged.
* **Axis ranges**: the report marks every axis as "given" or "automatic from the data"
  (`x -5-5 Ry (given), y 0.0003074-1.069e+04 (log axis) (automatic from the data)`), so
  after deleting `xlim`/`ylim` you still know which range was actually drawn.
* **The artifact list is computed**: the directory is compared before and after the run,
  listing the files that appeared or were rewritten with their "overwritten / created"
  marker, so you never have to guess what this run produced or whether it overwrote the
  previous one.

## How the three inputs line up

```
run dir/input.json                     pp-gen/src/pp_gen/parameters/
------------------                     -----------------------------
workflow   : "ld1.x"          -----+
parameters : {"input":             |   (the workflow registry: ae_radial needs the "input" group only)
              "ae_radial"}  -----+-->  input/ae_radial.yml --+
qe_path    : "xx"                  |                          |
                                   +-> pp_gen/inputs/workflows.py    +--> &input section
                                                                     |
                                                          run dir/U.in (<element>.in, no comments)
```

* **`workflow`**: a workflow name registered in `pp_gen/inputs/workflows.py`. Six are
  registered today: **`ld1.x`** (one input at a time, needing the `input` group only),
  **`rc_scan`** (`input` + `rc_scan_range`), **`configts`** (`input` + `configts`),
  **`rlderiv`** (`input` + `rlderiv`), **`rcloc_scan`** (`input` + `rcloc_scan_range`) and
  **`energy_scan`** (`input` + `energy_scan_range`). A misspelled workflow name, or a
  parameter group the workflow does not use, is an error right away — so that no section
  `ld1.x` should not read ever gets rendered.
* **`parameters`**: `group name -> YAML name`. `{"input": "ae_radial"}` means
  `parameters/input/ae_radial.yml` (the lookup order is
  `parameters/<group>/<name>.yml`, then `parameters/<name>.yml`, and finally a recursive
  search under the parameter root that succeeds only when it finds exactly one such file).
  The top level keys of the YAML are namelist section names, and the keys inside a section
  are the variables.
* **The file name** is decided by the `atom` of `&input` (when only `zed` is written, the
  element is recovered from the atomic number — the pslibrary inputs are written that way),
  so `atom: 'U'` → `U.in`.

## `qe_path`

`qe_path` points at QE's **`bin` directory** (recommended: `ld1.x` and `pw.x` both live
there), or directly at the `ld1.x` executable itself.

When `input.json` holds a placeholder (`"xx"`, empty, `TODO`…), `prepare` fills in the
default **and writes it back to `input.json`** (other keys stay as they are); a real path
already there is respected and left alone.

Defaults and overrides:

| | Value |
|---|---|
| Default QE path | `/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin` (`DEFAULT_QE_PATH` in `pp_gen/core/settings.py`) |
| Environment overrides | `PPGEN_QE_PATH` (QE path), `PPGEN_PARAMETERS` (root of the parameter YAMLs) |

The `ld1.x` under the default path on this machine is confirmed to exist and to be
executable; when it is not found, `prepare` only warns and does not stop (rendering input
files does not need QE).

## Install and test

```bash
python3 -m pip install -e /home/liguozhou/oncv/pp-gen     # installs the ppgen command
python3 -m pip install -e '.[plot]'                       # matplotlib is additionally needed for plotting
cd /home/liguozhou/oncv/pp-gen && python3 -m pytest -q     # 489 passed
python3 -m pytest -q -m "not slow"                        # 486 passed, 3 deselected, without touching real QE
```

It also runs without installing: `PYTHONPATH=src python3 -m pp_gen prepare -i .../input.json`.

`pytest -q` collects 489 tests and all of them pass on this machine; three carry the `slow`
marker because they really run the local `ld1.x`. Four tests that compare against the
reference artifacts of the earlier analysis (`.../uranium-ppgen/rc_scan/runs/5F/rc_1p1`, a
reference `ld1.wfc`, a reference scan product) skip themselves when those files are absent.

`run`'s tests mostly use a **fake `ld1.x`** (a bash script that reads stdin, writes stdout
and drops one file in the cwd), so exit codes, stderr and the working directory can be
asserted exactly; `plot`'s tests inspect the state of the matplotlib `Axes` (xlim, axis
labels, number of lines, legend, line styles/colours) instead of comparing PNG pixels. Only
the `slow` tests actually invoke the local QE.

## Verification

The generated `U.in` was really run through the local `ld1.x` (QE 7.4) — in practice
`ppgen run -i .../ae_radial/input.json`.

**Byte-for-byte identical to the reference run**: when `U.in` uses the same grid as the
reference file `U.ae.iswitch1.in` (`xmin=-8.0` / `dx=0.008` / `rmax=100`, i.e. the
automatic values of `ld1.x` for Z≥19), the products `ld1.wfc` and `vx.pot` are **byte for
byte the same** as the copies in `/home/liguozhou/oncv/ppgen/script/runs/ae_radial/`, and `U.out` differs only
in the timestamp, the available memory and the text of the title line.

> `ae_radial.yml` later wrote `dx=0.005` explicitly (a finer grid), so a run today gives an
> `ld1.wfc` of 551902 B / 3427 points instead of 344856 B / 2141 points. **That is the grid
> changing, not a failure to overwrite** — whether something was overwritten is what the
> "overwritten / created" marker in the `run` report and the file mtime are for.

The failure path was verified with real QE too: putting a variable name that does not exist
into `U.in` makes `ld1.x` exit with code 1, and `ppgen run` reports
`Error in routine ld1_readin (5010): reading input namelist` plus the tail of `U.out`, and
lists the `CRASH` and `input_tmp.in` files QE left behind.

The results of `plot` were compared item by item against the earlier analysis
(`ld1.wfc` → peak table):

* the largest `|rR|` of all seven states agrees **exactly** with
  `ae_radial_max_per_orbital.csv` (4F 1.7585, 5D 1.2444, 5F 0.8692, 6S 1.0017, 6P 0.8810,
  6D 0.5929, 7S 0.5720), which shows that reading the file and converting the units are
  correct;
* the peak positions (Å) differ from it by ≤ 0.004 Å (the difference comes from the
  `dx=0.005` grid versus the reference `dx=0.008`);
* after drawing, the image was checked pixel by pixel: all seven line colours of the seven
  states really are on the canvas (`ae_radial.png` is 1260×864 = 7.0×4.8 in @180 dpi), so no
  curve disappeared because of a wrong range;
* **the extrema table was compared row by row with the reference file**: on that reference
  `ld1.wfc`, the 28 extrema rows this tool computes are **exactly the same set of rows** as
  in `ae_radial_extrema.csv`, the largest deviation of `r_bohr`/`r_angstrom`/`rR` is 4e-11
  (the output precision of 10 decimals), and `n_nodes` and the `max`/`min` types match one
  by one — so the node/extremum algorithm is the same one used by the earlier analysis.
  With your finer `dx=0.005` grid, the row count, `n_nodes` and the types still correspond
  one to one, and the positions differ by ≤ 7e-3 bohr.
* **`rc_scan` was really run for a full pass**: all 55 scan points produced a
  `U.pbe-n-nc.UPF` (plus `ld1.dlog` / `ld1.pwe` / `ld1.wfc` / `ld1ps.*` / `vx.pot`), and the
  log ends with `End of Bessel function test`. Card rendering was also checked against
  `read_psconfig` of QE's `el_config.f90`: with `rel<2` only 7 columns are read (`rcutus`
  directly follows `rcut`, and `rcut <= rcutus` is required), the 8th column is the optional
  `isws`/`jj`, and extra columns are dropped by the list-directed read at the end of the
  record — so writing 8 columns the way this tool does is safe (and matches the reference
  `U.gen.in`), and at rendering time it already catches the "label does not match `l`,
  `n < l+1`, occupancy out of range, `rcut > rcutus`" classes of error for QE.
* **the pole diagnostics of the logarithmic derivative figure agree channel by channel with
  the reference script**: `pole_report` of `/home/liguozhou/oncv/ppgen/script/rc_scan_logderiv.py` was loaded unchanged and
  compared on the same point for all 4 channels against `n_ae/n_ps/matched/shifted/extra`
  for three `r_c` values (0.8, 1.1, 1.4) — **all identical**.
* **`rcloc_scan` was run once with real QE** (5 points scanning `rcloc` of `&inputp`: 2.0,
  2.25, 2.5, 2.75, 3.0, 1.5 s each, 5/5 succeeded): in `runs/rcloc/rcloc_2p5/U.in` only
  `rcloc` changed (`lloc`, the cards and `&test` were untouched); the criterion used **all
  channels**, `update` picked **2.25** (poles coincide, ΔE total 0.06 Ry, pooled `|ΔL|`
  median 1.4e-3) and wrote it back to the `&inputp` line; `report.md`'s header is
  `best rcloc (bohr)`, its channel column reads `l=all`, and it explains automatically that
  all channels were pooled.
* **a hole in the `|ΔE|` level was fixed**: `pairs` is built from the **PS** side, so extra
  AE poles used to be missing from it entirely — as a result a point where "PS has no pole
  in the window at all" got a `ΔE` sum of `0` (looking best) and was picked by the ranking.
  Now **a different count on the two sides is recorded as `inf`** (shown as "unmatched" in
  the report). The `rcloc=2.5 / 3.0` case above is exactly that situation: the smallest
  `|ΔL|` median (1.2e-3 / 4.6e-3) but mismatched AE/PS pole counts, and `--pole` no longer
  ranks them first.
* **`energy_scan` was run once with real QE** (9 points scanning the `ener` of the `5F`
  `occ=0` card, `-0.80…-0.40`, 1.9 s per point, 9/9 succeeded): the rendered `U.in` changed
  only that card's `ener` (the fifth, occupied `5F` card still reads `0.00` and the other
  `occ=0` card `6D` still reads `-0.05`); `update: true` wrote the picked best value back
  into the **7th** card (the other card with the same label was untouched); `output.json`
  records `"scanned": "ener"` / `"ener"`, `report.md`'s header is `best ener (Ry)`, the
  conclusion line reads `ener = -0.60 Ry`, and the x axis of the decision figure is
  `ener (Ry)` too.
* the logarithmic x axis (`"xscale": "log"`) was checked separately:
  `ax.get_xscale() == "log"`, major ticks land on the four decimal positions
  `0.01 / 0.1 / 1 / 10`, the y axis stays linear `[-2, 2]`, the curves span 85% of the canvas
  width, and **every column** of that span has a curve passing through it (so it is not just
  a few points being drawn).

## Module map

```
pp-gen/
├── pyproject.toml               # console script: ppgen = pp_gen.cli:main
├── src/pp_gen/
│   ├── __init__.py              # package docstring and __version__
│   ├── __main__.py              # `python -m pp_gen` entry point (no console script needed)
│   ├── core/                    # shared foundations
│   │   ├── __init__.py
│   │   ├── settings.py          # global settings: QE path, parameter root, env overrides, Bohr radius, timeout
│   │   ├── errors.py            # the domain errors the CLI turns into a readable message and exit code 2
│   │   ├── elements.py          # zed -> element symbol, used only to decide the input file name
│   │   └── cards.py             # field layout of the orbital cards (columns, card-line pattern, labels)
│   ├── inputs/                  # the three inputs -> ld1.x input text
│   │   ├── __init__.py
│   │   ├── request.py           # input.json: read (order preserving), validate, write back
│   │   ├── parameter_files.py   # locate, load and merge the parameter YAMLs
│   │   ├── workflows.py         # workflow registry: name -> required parameter groups, scan axis
│   │   ├── ld1_input_writer.py  # render parameter sections into ld1.x input text (comment free)
│   │   └── generation.py        # the prepare use case: input.json -> <element>.in, or one per scan point
│   ├── scans/                   # the scan axes and their points
│   │   ├── __init__.py
│   │   ├── settings.py          # names of the scan layer: runs/, update, the auxiliary range keys
│   │   ├── ranges.py            # a range specification -> the concrete values
│   │   ├── axis.py              # which quantity a scan walks over: headings, units, directory tags
│   │   ├── points.py            # the scan points themselves (one per rendered input, card keys)
│   │   └── apply.py             # apply one scan point to the parameter sections (pure, no file access)
│   ├── ld1/                     # running ld1.x and reading what it wrote
│   │   ├── __init__.py
│   │   ├── runner.py            # run ld1.x (direct or point by point): exit code, timing, artifacts
│   │   ├── wavefunctions.py     # read ld1.wfc and derive rR / R / the probability density
│   │   ├── log_derivatives.py   # read ld1.dlog / ld1ps.dlog, pole criterion, AE/PS classification
│   │   ├── extrema.py           # nodes and extrema of the radial wavefunctions + the extrema CSV
│   │   └── records.py           # read back output.json tasks and the ld1.test energy blocks
│   ├── diagnostics/             # choosing a value from a finished scan
│   │   ├── __init__.py
│   │   └── selection.py         # `update: true`: pick the best value of a group and decide what to write back
│   ├── plotting/                # everything that draws
│   │   ├── __init__.py
│   │   ├── spec.py              # plot.json: read, validate, and infer (ylabel -> quantity, xlabel -> unit)
│   │   ├── dispatch.py          # the plot entry point: pick the figure body from `mode`
│   │   ├── builder.py           # drawing only: palette, line styles, legend, saving
│   │   ├── wavefunction.py      # the radial wavefunction figure (rR(r), R(r) or the density)
│   │   ├── log_derivatives.py   # the AE/PS logarithmic derivative figure and its parser.json
│   │   ├── decision.py          # report's "why this scan value" figure
│   │   ├── checks.py            # warnings and side outputs collected while drawing (limits, CSV)
│   │   ├── orbitals.py          # data.orbitals filtering and the (PS)/(AE) variants of a state
│   │   └── series.py            # data holders the figure bodies return and the report consumes
│   ├── reporting/               # scan results -> report.md
│   │   ├── __init__.py
│   │   ├── metrics.py           # the numbers behind the report: one entry per scan point, ranking keys
│   │   ├── markdown.py          # render the report body as markdown (scan and configts variants)
│   │   └── report.py            # the entry points that assemble a report from a run directory
│   ├── writing/                 # write a chosen value back into the parameters
│   │   ├── __init__.py
│   │   ├── cards.py             # rewrite rcut/rcutus/ener inside the orbital cards, in place
│   │   └── namelists.py         # rewrite a namelist scalar (rcloc / rlderiv) in place
│   ├── cli/                     # the command line
│   │   ├── __init__.py          # `main`: build the parser, dispatch, report errors and exit codes
│   │   ├── arguments.py         # the argparse tree: subcommands, options and their defaults
│   │   ├── console.py           # small formatting helpers shared by the command modules
│   │   ├── filters.py           # parse --only / --rc / --ener into lists
│   │   └── commands/
│   │       ├── prepare.py       # the prepare command: render inputs and print what was written
│   │       ├── run.py           # the run command: direct run or generated run, plus updates
│   │       ├── plot.py          # the plot command: draw the figure plot.json describes
│   │       └── report.py        # the report command: summarise a scan into report.md
│   └── parameters/              # parameter YAMLs: input/, rc_scan_range.yml, rlderiv.yml, rcloc_scan_range.yml …
└── test/                        # pytest: rendering, run, plot, report, error paths, CLI
```

Command line → module: `ppgen prepare` → `inputs/generation.py` (`cli/commands/prepare.py`),
`ppgen run` → `ld1/runner.py` (`cli/commands/run.py`), `ppgen plot` →
`plotting/dispatch.py` (`cli/commands/plot.py`), `ppgen report` → `reporting/report.py`
(`cli/commands/report.py`). The scan value chosen by `diagnostics/selection.py` is written
back through `writing/cards.py` (a card column) or `writing/namelists.py` (a namelist
scalar).
