#! /usr/bin/env python
import hashlib
import os.path
from os import path

md5_hashes = {
	"metrics_HACC_Travis_.csv" : "d424f552b91c360e79e7f9bba44fcc6c",
	"HACC_Travis_/SZ_HACC_pos_vel___m000.full.mpicosmo.50" : "90a620e32a09202389cfdbd24ead6ec0",
	"metrics_NYX_Travis_.csv": "382c19f8b6cf0b39ab5c3cfcef29b2de",
	"NYX_Travis_/SZ___z255_32.h5" : "db4a04885122ff9957f60c6102325528"
}


def computeMD5Hash(file):
	if ( not path.exists(file) ):
		print file + " does not exit!!!!!"
		return ""

	with open(file, 'rb') as afile:
    		buf = afile.read()
    		hasher.update(buf)
		return hasher.hexdigest()

if __name__ == "__main__":

	all_passed = True

	for item in md5_hashes:
		computed_hash = computeMD5Hash(item)

		if (md5_hashes[item] != computed_hash):
			print "   Error in " + item + "\n"
			all_passed = False
		else:
			print "   " + item + " is identical\n"
		

	if all_passed == True:
		print "All tests passed!"
