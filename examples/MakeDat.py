#!/usr/bin/env python

# makes dat text files with some of the information contained in the SkyNet h5
# output files in the given directory

import os
import sys
import itertools
from SkyNet import *
import numpy as np
import multiprocessing

num_done = multiprocessing.Value('i', 0)

def H5ToDat(args):
  global num_done
  fname, total_size = args

  if os.path.exists(fname + ".dat"):
    print(fname + ": skipped")
    return

  NetworkOutput.MakeDatFile(fname)

  num_done.value = num_done.value + 1
  print("[%i/%i] done %s" % (num_done.value, total_size, fname))


if __name__ == '__main__':
  if (len(sys.argv) != 2):
    print("Usage: %s <directory with .h5 files>" % sys.argv[0])
    sys.exit(1)

  output_dir = sys.argv[1]

  num_cores = multiprocessing.cpu_count()
  print("Running with %i worker threads" % num_cores)

  pool = multiprocessing.Pool(num_cores, maxtasksperchild=1)

  big_files = []
  size_limit = 1024 * 1024 * 1024 # 1 GB

  files = []
  num_done.value = 0

  for file in os.listdir(output_dir):
    if (not file.endswith(".h5")):
      continue

    if (file.startswith("FAILED")):
      continue

    # check file size
    path = output_dir + "/" + file
    if (os.path.getsize(path) > size_limit):
      print("%s is big, defer until later" % path)
      big_files.append(path)
    else:
      files.append(path)

  size = len(files)
  args = [(p, size) for p in files]
  result = pool.map_async(H5ToDat, args)
  result.wait()

  # now do big files serially because they require lots of memory
  num_done.value = 0
  print("now doing big files individually")
  for f in big_files:
    H5ToDat((f, len(big_files)))
