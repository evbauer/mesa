# TASK070 regeneration

Run from the workspace root only on the inherited compute node (`rzwhippet25`, job `5549798`). All numerical stages enforce placement and stop at 00:25 PST on 3 November 2026. Use one numerical-library thread per process:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
code=mesa/eos/eosFreeEOS_builder/biquintic/freeeos_helm_map
python3 -B "$code/run.py" init
python3 -B "$code/run.py" prepare
python3 -B "$code/run.py" prepare_routes
python3 -B "$code/run.py" route_T
python3 -B "$code/run.py" route_d
python3 -B "$code/run.py" kernel
python3 -B "$code/run.py" kernel_d
python3 -B "$code/run.py" freeeos
python3 -B "$code/analyze.py"
python3 -B "$code/validate.py"
python3 -B "$code/diagnostics.py"
python3 -B "$code/display.py"
python3 -B "$code/norm_audit.py"
python3 -B "$code/plot.py" --stem freeeos_helm_native_faces_corrected
python3 -B -m unittest discover -s "$code" -p 'test_*.py' -v
```

These exact output paths already exist. The stages deliberately refuse to overwrite raw streams, audits, and plot products. To regenerate an EOS survey, use a fresh experiment directory and adjust the single `OUTPUT` constant in a copied scoped code directory; preserve experiment 070 as evidence. For a plot-only redraw of the existing arrays, run `python3 -B "$code/plot.py" --stem NEW_NAME` from the workspace root. It reads saved arrays and `display_bounds_corrected.json`, never either EOS. The corrected numerical analysis consumes `routing_T_corrected/`, `routing_d_corrected/`, and the linear-coordinate `helm_{T,d}_kernel.tsv` files. Its deliverable is `helm_T_traces_corrected.npz` and `helm_summary_corrected.json`; `helm_d_traces.npz` is analytically zero on valid external density faces. `display.py` replays the original 060 plotted points and margins, then crops the saved HELM arrays; `norm_audit.py` reads archived 048 traces without source evaluations.

The initial `routing_T/` and `routing_d/` outputs are preserved failed probes: linear coordinates were supplied to a logarithmic probe, producing status `-1003`. The first `helm_T_traces.npz` and `helm_summary.json` are preserved failed analyses from swapped mixed-derivative slots. Never use these for the figure. `freeeos_helm_native_faces_corrected.{png,pdf}` is the proposed Figure 4. The earlier `_delivery`, `_final`, and unqualified figure files are presentation drafts: `_delivery` used the full 060 manifest bounds rather than the old plot's autoscaled bounds. The figure's own presentation metadata names and hashes its plotted arrays; `in_display_stats_corrected.json` distinguishes cropped populations from the broad survey.

The external HELM selection is measured at each new face midpoint with the original-library MESA probe and exact scaled-solar mixture. It is a diagnostic selection, not an assertion about default-MESA HELM routing. Raw `probe.tsv` and `helm_*_kernel.tsv` files, coordinates, one-sided traces, masks, summaries, checks and hashes remain under experiment 070. No MESA source or installed cache is changed; all new runtime caches are private.

Compile the editable task PDF from `freeeos_biquintic_progress/progress_reports/` with `pdflatex -interaction=nonstopmode -halt-on-error 2026-11-02_freeeos_helm_map.tex` (twice for stable references). After the final report and command log are frozen, `python3 -B "$code/provenance.py"` saves their SHA-256 hashes in `provenance_final.json`; that command also refuses to overwrite an existing manifest.
