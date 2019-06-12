#! /usr/bin/env python
import hashlib
import sys
import os.path
from os import path

 
hasher = hashlib.md5()
with open(sys.argv[1], 'rb') as afile:
    buf = afile.read()
    hasher.update(buf)

print(hasher.hexdigest())

