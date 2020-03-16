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

from draw.cinema import cinema as cnm

class NYXCinema(workflow.CinemaCreator):
	def prepare_cinema():
		# copy compression param from reduction pre-compressed to csv
  
		# Copy to folder cinema 
  
		pass