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
import os

import Workflow.draw.utils as utils
import Workflow.draw.cinema as cnm

class NYXCinema(cnm.CinemaCreator):
	def prepare_cinema():
		# copy compression param from reduction pre-compressed to csv

		# Copy to folder cinema 
		cbench_csv_path = self.json_data['project-home'] + self.json_data['wflow-path'] + "/reduction/cbench/" + self.json_data['data-reduction']['cbench-output']['metrics'] + ".csv"
		csv_dst = self.json_data['project-home'] + self.json_data['wflow-path'] + "/cinema/"
  
		utils.run_command("cp " + cbench_csv_path + " " + csv_dst)


# Parse Input
parser = argparse.ArgumentParser()
parser.add_argument("--input-file")
opts = parser.parse_args()


# Create Cinema DB
cinema = NYXCinema( opts.input_file )
cinema.create_cinema()