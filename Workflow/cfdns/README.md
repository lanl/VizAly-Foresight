## Converstion data Steps:
1. hdf5 to binary;
	python3 dataConversion/hdf5_to_raw.py original_data scalar1 scalar2 ... min_timestep max_timestep output folder
	e.g. python3 dataConversion/hdf5_to_raw.py /bigData/Turbulence/scalarHIT_fields100.h5 one two vx vy vz 0 50  data/conversion3 

2. create foresight input
	python3 dataConversion/createForesightInput.py /home/pascal/projects/VizAly-Foresight/inputs/turbulence start_timestep numTimesteps gda-turbulence_template.json
	e.g. python3 dataConversion/createForesightInput.py /home/pascal/projects/VizAly-Foresight/inputs/turbulence 0 50 gda-turbulence_template.json

3. Run foresight
	python3 dataConversion/runCbench.py

4. convert back to hdf5
	python3 dataConversion/raw_to_hdf5.py input_file dims_x dims_y dims_z data_folder min_timestep max_timestep scalar1 scalar2 scalar3 scalar4 scalar5 prefix
	e.g. python3 dataConversion/raw_to_hdf5.py raw50.h5 128 128 128 data/gdaOriginalFiles 0 50 one two vx vy vz ts_


## Analytics step
1. Run inference engine
	# python3 turbulence_analytics/inference.py inputData checkpointspath outputPath
	# python3 turbulence_analytics/inference.py /home/pascal/projects/VizAly-Foresight/Analysis/pat/turbulence/raw50.h5 turbulence_analytics/checkpoints /home/pascal/projects/VizAly-Foresight/Analysis/pat/turbulence/data/conversion2

2. Plot the results
	# python3 turbulence_analytics/plot_all.py data_path output_prefix output_path
	# python3 turbulence_analytics/plot_all.py /home/pascal/projects/VizAly-Foresight/Analysis/pat/turbulence/data/test1 original original_plot



 

1. hdf5 to binary;
python3 dataConversion/hdf5_to_raw.py original_data scalar1 scalar2 ... min_timestep max_timestep output folder
python3 dataConversion/hdf5_to_raw.py /bigData/Turbulence/scalarHIT_fields100.h5 one two vx vy vz 0 50  data/turbulence


2. convert back to hdf
python3 dataConversion/raw_to_hdf5.py output_file dims_x dims_y dims_z data_folder min_timestep max_timestep scalar1 scalar2 scalar3 scalar4 scalar5 prefix
e.g. python3 dataConversion/raw_to_hdf5.py raw50.h5 128 128 128 data/turbulence 0 50 one two vx vy vz ts_



python3 turbulence_analytics/debugInference.py /home/pascal/projects/VizAly-Foresight/Analysis/pat/turbulence/raw50.h5 /bigData/Turbulence/scalarHIT_fields100.h5


Friday Oct 25 @ 10:00
python dataConversion/hdf5_to_raw2.py /projects/ml_compression/data/scalarHIT_fields100.h5 one two vx vy vz 0 50  temp
python dataConversion/raw_to_hdf52.py testAll50.h5 128 128 128 temp 0 50 one two vx vy vz turbulence_ts_
python dataConversion/compareHDF5.py testAll50.h5 ../../../../data/scalarHIT_fields100.h5



## Running Script
ssh to Darwin and then:
salloc -N 1 -p volta-x86
cd /projects/ml_compression/VizAly-Foresight/Analysis/pat/turbulence
source /projects/ml_compression/VizAly-Foresight/scripts/VizAly-CBench.bash.darwin-autoencoder
mkdir diagnosticsCAE
python turbulence_analytics/debugInference.py /projects/ml_compression/data/scalarHIT_fields100.h5 /projects/ml_compression/data/raw50.h5
python turbulence_analytics/plot_all.py diagnosticsCAE test1 autoencoder_Oct_23


Also, running "python turbulence/dataConversion/compareHDF5.py /projects/ml_compression/data/raw50.h5 /projects/ml_compression/data/scalarHIT_fields100.h5" will indicate that these two files are identical.

Steps that I used to create raw50.h5 if you want to recreate the data:
python turbulence/dataConversion/hdf5_to_raw.py /projects/ml_compression/data/scalarHIT_fields100.h5 one two vx vy vz 0 50 /projects/ml_compression/data/conversion
python turbulence/dataConversion/raw_to_hdf5.py raw50.h5 128 128 128 /projects/ml_compression/data/conversion 0 50 one two vx vy vz ts_



# Steps:
python3 /usr/projects/ml_compression/VizAly-Foresight-CBench/Analysis/pat/turbulence/dataConversion/hdf5_to_raw2.py /usr/projects/ml_compression/data/original/scalarHIT_fields100.h5 one two vx vy vz 0 50 data/convesion3