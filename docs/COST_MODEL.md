# Cost model and calibration

With Courant factor $C$ the time step is $\Delta t = C/r$ (Meep, 2D), so a run of $T$ time units
on an $s_x \times s_y$ cell at resolution $r$ costs

$$W = N_x N_y N_t = (s_x r)(s_y r)\,\frac{T}{\Delta t} = \frac{s_x s_y\,T\,r^{3}}{C}\ \text{pixel-steps},$$

which is $2\,s_x s_y T r^3$ for $C = 0.5$. When the mirror sector is used, half the cell is computed
and $W$ is halved. $T$ is the source duration $2\,{\rm cutoff}/f_{\rm width}$ plus the measurement
time: the harmonic-inversion signal $t$ (or the second pass of the two-pass rule), the spectrum
ceiling $\min(t_{\max}, 6 Q_{\rm est}/\pi f_{\rm cen})$, the field ceiling $\min(t_{\max}, 4 Q_{\rm est}/\pi f_{\rm res})$,
or 600 time units for the cavity-less reference.

The reference computation of the paper ($30a \times 29a$, $r = 24$, $t = 17\,750$, mirror sector)
is $W = 2.16\times10^{11}$ pixel-steps: one reference-run unit (RRU). It took 677 s on one core of
an AMD EPYC 7742 node that ran 32 such tasks at once, i.e. $3.19\times10^{8}$ pixel-steps per
second per task, the built-in default.

Memory per task is estimated as 64 bytes per pixel (six field arrays and two permittivity
arrays) plus the DFT buffers of the flux planes, plus 300 MB for the interpreter; the runner caps
the number of workers so that the sum fits the available memory.

`sprat calibrate` runs the reference geometry with a short signal at resolutions 12 and 24 and
writes `calibration.json` (throughput of the resolution-24 sample, host, cores, date) into the
records directory; `plan`, `run` and `reproduce` read it. The estimate is an upper bound on
efficiency: small cells run at a lower pixel-step rate (the 693 coarse-sweep cells of the
systematic runs of 2026, $17a \times 15a$ at resolution 20, ran at $3.6\times10^{7}$ pixel-steps/s),
and the rate drops when many tasks share the memory bandwidth of one socket (the 14 spectrum runs
with the guide continued into the absorber, 14 concurrent tasks of four cores each on one node,
ran at $5.6\times10^{7}$ to $1.1\times10^{8}$ pixel-steps/s per task).

The wall time of a batch with $P$ workers is $\max(\max_i T_i,\ \sum_i T_i / P)$ for longest-first
scheduling; `sprat plan` prints it for several $P$.
