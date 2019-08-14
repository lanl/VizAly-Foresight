Scripts to load modules on supercomputers:
1. Darwin @ LANL: VizAly-CBench.bash.darwin 
2. Cooley @ ANL : VizAly-CBench.bash.cooley
3. Cori @ NERSC: VizAly-CBench.bash.cori
4. Darwin @ LANL with GPU : VizAly-CBench.bash.darwin.gpu

To allocarte a GPU node on Darwin interactively:
salloc -p shared-gpu

It also has scripts to run CBench on supercomputers, e.g  runJob_cori_NYX_test.sh is an example batch script for darwin

A working [SPACK](https://github.com/spack/spack) script is included, to automatically compile on local environments, supercomputers, e.g.
