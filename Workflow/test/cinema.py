#! /usr/bin/env python
"""
pat_nyx fun the foresight pipeline for nyx analysis. In other words, it does:
- run cbench
- run the analysis tool
- create a cinema databse of the result

To run:
python pat_nyx.py --input-file <absolute path to input file> <--submit>

--submit: causes the pipeline to be launched. Without it, you can see the scripts generated for debugging.

e.g.
python -m tests.workflow --input-file
"""

import argparse
import csv
import os

import Workflow.draw.utils as utils
import Workflow.draw.cinema as cnm

class NYXCinema(cnm.CinemaCreator):

	def prepare_cinema(self):
		# Open CSV file
		# if "metrics-csv" in self.json_data["input"]:
		# 	metrics_csv = self.json_data["input"]["metrics-csv"]
		# else:
		# 	metrics_csv = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/cbench/" + self.json_data['cbench']['output']['metrics-file'] + ".csv"

		
  
		# Copy to folder cinema 
		csv_filename = self.json_data['data-reduction']['cbench-output']['metrics-file']
		cbench_csv_path = self.json_data['project-home'] + self.json_data['wflow-path'] + "/reduction/cbench/" + csv_filename + ".csv"
		csv_dst = self.json_data['project-home'] + self.json_data['wflow-path'] + "/cinema/data.csv"
  
		print("cp " + cbench_csv_path + " " + csv_dst)
		utils.run_command("cp " + cbench_csv_path + " " + csv_dst)
  
		output_file_name = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/cinema/data.csv"


		all = []
		#reader = futils.open_csv_file(metrics_csv)
		with open(csv_dst,'r') as csvinput:
			reader = csv.reader(csvinput)

			# Modify Cinema files
			row = next(reader)
			row.append('FILE_SimStats_Pk')
			row.append('FILE_lya_all_axes_x_Pk')
			row.append('FILE_lya_all_axes_y_Pk')
			row.append('FILE_lya_all_axes_z_Pk')
			all.append(row)

			values = ["sim_stats_rhob.png", "sim_stats_rhodm.png", "sim_stats_temp.png", "sim_stats_velmag.png", "sim_stats_velmag.png", "sim_stats_vz.png"]
			count = 0
			for row in reader:
				row.append(values[count])
				row.append("lya_all_axes_x.png")
				row.append("lya_all_axes_y.png")
				row.append("lya_all_axes_z.png")
				all.append(row)

				count = count + 1
				if (count == 6):
					count = 0

			
			utils.write_csv(output_file_name, all)

# Parse Input
parser = argparse.ArgumentParser()
parser.add_argument("--input-file")
opts = parser.parse_args()


# Create Cinema DB
cinema = NYXCinema( opts.input_file )
cinema.create_cinema()