## Python codes
List of Python codes.

### GPU-enabled LDOS scripts
The GPU variants are `LDOS_kwant_SOC_supercell_gpu.py` and
`LDOS_kwant_nSOC_supercell_gpu.py`. Run them with the `--use-gpu` flag to route
the hydrogen-coupling block math through CuPy (Kwant remains CPU-bound).

#### nSOC (spin up/down)
```bash
python3.10 LDOS_kwant_nSOC_supercell_gpu.py \
  $INPF_TB_UP $OUTPUT_PY_DATA_UP \
  $INPF_TB_DN $OUTPUT_PY_DATA_DN \
  --use-gpu >> $mainf/$OUTF
```

#### SOC
```bash
python3.10 LDOS_kwant_SOC_supercell_gpu.py \
  $INPF_TB $OUTPUT_PY_DATA \
  --use-gpu >> $mainf/$OUTF
```
