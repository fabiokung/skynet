#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
  print("Usage: %s <old REACLIB file>" % sys.argv[0])
  sys.exit(1)

fin = open(sys.argv[1])

while (True):
  l = fin.readline()
  l = l.replace("D+", "E+")
  l = l.replace("D-", "E-")

  if (l == ""):
    # reached end of file
    break

  if ((len(l) == 2) or (len(l) == 3)):
    chapter = int(l)
    # skip next two lines (they should be empty)
    fin.readline()
    fin.readline()
  else:
    sys.stdout.write("%i\n" % chapter)
    sys.stdout.write(l[0:-1].ljust(74) + "\n")

    l = fin.readline()
    l = l.replace("D+", "E+")
    l = l.replace("D-", "E-")
    sys.stdout.write(l[0:-1].ljust(74) + "\n")

    l = fin.readline()
    l = l.replace("D+", "E+")
    l = l.replace("D-", "E-")
    sys.stdout.write(l[0:-1].ljust(74) + "\n")

fin.close()
