# Running the Phase 1 stability test on CENAPAD-SP (Lovelace)

## The machine layout, which is not what the guide's wording suggests

`cenapad.unicamp.br` answers as a host named **`frontend`** and is a
*connection server only*: `module` and `qstat` refuse to run there. The work
machine is reached with `ssh lovelace` and is where you edit, compile and
submit. The guide calls that the "Lovelace frontend", which reads backwards
from the hostnames.

```
your laptop  --ssh -p 31459-->  frontend  --ssh lovelace-->  lovelace
                                (auth only)                  (build, qsub)
```

## Environment, read off the machine

| | |
|---|---|
| system gcc | 8.5.0 — **too old**, Castro needs C++20 |
| toolchain used | `module load openmpi/5.0.6-gcc-12.2.0` |
| scheduler | PBS (`qsub`, `qstat -u $USER`, `qdel`) |
| `git` | present, and has outbound internet (`git ls-remote` works) |
| `$HOME` | CephFS, 212 TB free |
| `/work` | CephFS, 41 TB free |

Queue state when this was written — the reason the jobs target `parexp`:

| queue | walltime | nodes | running | **queued** |
|---|---|---|---|---|
| `par128` | 168 h | 1 (128 cores) | 35 | **91** |
| `paralela` | 72 h | 10 | 6 | 22 |
| `parexp` | 24 h | 11 (48 cores) | 4 | **0** |
| `testes` | 1 h | 1 | 0 | **0** |

`par128` has 2.7x the cores per node and a 91-deep queue. For a run estimated
at one to two hours, `parexp` starting now beats `par128` starting tomorrow.

## Procedure

From the laptop, in the repository root:

```bash
scp -P 31459 -r cluster/cenapad rcrlima@cenapad.unicamp.br:~/
scp -P 31459 models/phase1_scenario.txt models/phase1_control.txt \
    rcrlima@cenapad.unicamp.br:~/
```

Then on lovelace:

```bash
ssh -p 31459 rcrlima@cenapad.unicamp.br
ssh lovelace

bash ~/cenapad/bootstrap.sh ~/wd-mag        # clones Castro 26.07, pins the
                                            # submodules, re-applies the two
                                            # core patches, installs the problem
cp ~/phase1_*.txt ~/wd-mag/Castro/Exec/science/wd_scf_stability/
cd ~/wd-mag/Castro/Exec/science/wd_scf_stability
bash ~/cenapad/build.sh .

qsub ~/cenapad/job_ic_check.pbs             # smoke test first, queue testes
qstat -u $USER
# read wdscf96ic.out: it must show div B ~ 1e-16 and the density peak on target

qsub ~/cenapad/job_stability.pbs            # production, queue parexp
```

## The core patches are not optional

`castro_core.patch` touches two files in Castro's own `Source/`:

- `Source/gravity/Gravity.cpp` — without it, the half-shift geometry that puts
  a cell *centre* at r = 0 makes Castro divide by zero and abort with NaN
  density on the first step.
- `Source/driver/Castro_io.cpp` — without it the `ztwd` EOS does not build
  (`network_rp has not been declared`), because `ztwd` has no runtime
  parameters and so never pulls in `extern_parameters.H` transitively.

`bootstrap.sh` applies them and is idempotent. A fresh `git clone` or
`git submodule update` loses them.

## What the run is for

See `docs/PLAN_ULTRAMASSIVE_SNIA.md`, Phase 1. Short version: the 2 M☉
configuration is toroidal-dominated (E_tor/E_pol ~ 3e6 at the scenario's 1e9 G
surface dipole), which is the regime most exposed to the Tayler m = 1
instability. If it rearranges on an Alfven time (1.26 s) it can never sit
still long enough for magnetic braking to act over Myr, and the whole SN Ia
scenario needs a different field. The run damps for 1.68 s to let the star
settle onto the Cartesian mesh, then switches damping OFF and measures for ten
Alfven times.
