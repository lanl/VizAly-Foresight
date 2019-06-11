#! /usr/bin/env python
import hashlib
import os.path
from os import path

md5_hashes = {
	"metrics_HACC_Travis_.csv" : "75763e77c6880ab6c49bd1779d6aaa19",
	"HACC_Travis_/BLOSC_HACC___m000.full.mpicosmo.50" : "5c31f4c4639788ae1d69ffeb36c859a4",
	"HACC_Travis_/SZ_HACC_pos_vel___m000.full.mpicosmo.50" : "f8e1dd7a7aab0dca39440b80f3b733e9",
	"metrics_NYX_Travis_.csv": "193eac46691fd68c73844b83118554df"
}

def computeMD5Hash(file):
	if ( not path.exists(file) ):
		print file + " does not exit!!!!!"
		return ""
	return hashlib.md5(file).hexdigest()

if __name__ == "__main__":

	all_passed = True

	for item in md5_hashes:
		computed_hash = computeMD5Hash(item)

		if (md5_hashes[item] != computed_hash):
			print "   Error in " + md5_hashes[item] + "\n"
			all_passed = False
		

	if all_passed == True:
		print "All tests passed!"
