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
| toolchain used | `module load openmpi/5.0.8-gcc-15.2.0` |
| why not GCC 12.2 | Castro 26.07 includes `<format>`, which libstdc++ only has from GCC 13 |
| `python` | absent; only `python3`. `build.sh` shims it, Castro's scripts need `python` |
| home directories | frontend and lovelace are SEPARATE filesystems -- scp needs two hops |
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

The repository is public and lovelace's `git` has outbound internet, so the
source does not need transferring -- clone it there. Only the model files move
by `scp`, because they are gitignored derived data (15 MB each).

On lovelace:

```bash
git clone https://github.com/RafaelCRdeLima/WD_MAG.git
bash WD_MAG/cluster/cenapad/bootstrap.sh ~/wd-mag
```

From the laptop, the only transfer -- note the absolute path, the repository
is at ~/wd-magnetizada and ~/Codes/WD_MAG merely contains it:

```bash
scp -P 31459 /home/rafael/wd-magnetizada/models/phase1_scenario.txt \
             /home/rafael/wd-magnetizada/models/phase1_control.txt \
    rcrlima@cenapad.unicamp.br:/home/lovelace/proj/proj503/rcrlima/wd-mag/Castro/Exec/science/wd_scf_stability/
```

**Check first that this lands somewhere lovelace can see it.** The two machines
report different Ceph monitors and different mount subpaths -- frontend has
`/home` (29 TB) from 192.168.193.x mounted at `:/`, lovelace has
`/home/lovelace` (759 TB) from 172.27.254.x mounted at `:/lovelace/home`. They
may not be the same filesystem. Test by writing a file on lovelace and looking
for it from the frontend. If they are separate, copy from lovelace in a second
hop rather than pushing from the laptop.

Then, back on lovelace:

```bash
cp ~/phase1_*.txt ~/wd-mag/Castro/Exec/science/wd_scf_stability/
cd ~/wd-mag/Castro/Exec/science/wd_scf_stability
bash ~/wd-mag/../WD_MAG/cluster/cenapad/build.sh .

qsub ~/WD_MAG/cluster/cenapad/job_ic_check.pbs     # smoke test, queue testes
qstat -u $USER
# read wdscf96ic.out: density peak on target, max|div B| h/|B| ~ 1e-16

qsub ~/WD_MAG/cluster/cenapad/job_stability.pbs    # production, queue parexp
```

## Watching a run without typing qstat repeatedly

Three ways, roughly in order of usefulness:

```bash
bash ~/WD_MAG/cluster/cenapad/watch_job.sh    # state, steps, rho_c drift
tail -f run_prod96.log                        # the raw output, live
qsub -m abe -M you@example.com job_stability.pbs   # mail on begin/end/abort
```

The job redirects into `run_prod96.log` in the working directory precisely so
it can be followed live. PBS only delivers its own `-o` file when the job
ENDS, which makes a multi-hour run invisible until it is over.

`-m abe` is set at submit time rather than in the script, so no address is
committed to a public repository.

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
