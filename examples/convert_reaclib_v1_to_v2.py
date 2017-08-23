#!/usr/bin/env python

import sys

def get_line(fin):
  l = fin.readline()

  # fix Windows line endings
  l = l.replace("\r\n", "\n")

  # fix exponential format with D instead of E
  l = l.replace("D+", "E+")
  l = l.replace("D-", "E-")
  l = l.replace("d+", "e+")
  l = l.replace("d-", "e-")

  return l

if len(sys.argv) != 2:
  print("Usage: %s <old REACLIB file>" % sys.argv[0])
  sys.exit(1)

fin = open(sys.argv[1])

while (True):
  l = get_line(fin)

  if (l == ""):
    # reached end of file
    break

  if ((len(l.strip()) == 1) or (len(l.strip()) == 2)):
    chapter = int(l)
    # skip next two lines (they should be empty)
    fin.readline()
    fin.readline()
  else:
    sys.stdout.write("%i\n" % chapter)
    sys.stdout.write(l[0:-1].ljust(74) + "\n")

    l = get_line(fin)
    sys.stdout.write(l[0:-1].ljust(74) + "\n")

    l = get_line(fin)
    sys.stdout.write(l[0:-1].ljust(74) + "\n")

fin.close()
