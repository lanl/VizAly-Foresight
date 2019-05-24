NYX_reamde.md

## gimletAnalysis.py
Add missing attributes to the compressed hdf5 files and run the gimlet analysis sim_stats.ex

## gimletPowerSpectrum.py
Creates a ratio of different metrics used

A sample input is specified in gimletInput.json


## Questions:
/native_fields:
baryon_density
dark_matter_density
temperature
vx
vy
vz

gimlet (sim_stats.ex) metrics:
 _rhodm_ps3d.txt - dark_matter_density
  _temp_ps3d.txt - temperature
  _rhob_ps3d.txt - baryon_density
  _rhom_ps3d.txt - ????
_velmag_ps3d.txt - vx, vy, vz
    _vz_ps3d.txt - vz


### What is _rhom_ps3d.txt"


Apps in gimlet:
1. lya_all_axes
2. lya_all_axes_rhoT
3. master_pk
4. p1d_all_axes
5. process_gadget_lya
6. sim_stats

### Which app to use?