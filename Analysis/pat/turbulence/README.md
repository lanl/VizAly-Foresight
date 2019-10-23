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