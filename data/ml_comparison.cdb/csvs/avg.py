import csv

def average_csv(file_prefix, num_files, output_filename):
	filename = file_prefix + "0.csv"

	master_csv = []
	num_rows = 0
	num_cols_to_skip = 2


	# Read first file into a list of lists
	with open(filename) as csv_file:
		reader = csv.reader(csv_file)
		for row in reader:
			num_cols = len(row)
			num_rows = num_rows + 1
			master_csv.append(row)


	# Read other files and add
	for i in range(num_files-1):
		filename = file_prefix + str(i+1) + ".csv"

		with open(filename) as csv_file:

			row_count = 0
			reader = csv.reader(csv_file)
			for row in reader:

				#skip header
				if row_count == 0:		
					row_count = row_count + 1		
					continue


				col_count = 0
				for c in range(num_cols):
					# skip the first two cols
					if c < num_cols_to_skip:
						col_count = col_count + 1
						continue

					# Do addition
					if row[c] != 'inf':
						master_csv[row_count][col_count] = str( float(row[c]) + float(master_csv[row_count][col_count]) )

					col_count = col_count + 1

				row_count = row_count + 1

	# Write to CSV file
	with open("sum.csv", 'w', newline='\n') as file:
		writer = csv.writer(file)
		for row in range(num_rows):
			writer.writerow(master_csv[row])


	# Do average
	for r in range(num_rows):
		if r == 0:
			continue

		for c in range(num_cols):
			# skip the first two cols
			if c < num_cols_to_skip:
				col_count = col_count + 1
				continue

			# Do addition
			if master_csv[r][c] != 'inf':
				master_csv[r][c] = " " + str( float(master_csv[r][c])/num_files )

			col_count = col_count + 1


	# Write to CSV file
	with open(output_filename, 'w', newline='\n') as file:
		writer = csv.writer(file)
		for row in range(num_rows):
			writer.writerow(master_csv[row])
            
average_csv("metrics", 50, "data.csv")