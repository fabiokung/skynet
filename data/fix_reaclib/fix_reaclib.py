#!/usr/bin/env python

from SkyNet import *
import numpy as np
import sys

old_reaclib = sys.argv[1]
fixed_reaclib = old_reaclib + "_fixed"

# get nuclear data
nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
    + "/data/webnucleo_nuc_v2.0.xml")

# get the list of bad nuclides
bad_beta_minus_file = open("bad_beta_minus", "r")
bad_beta_minus = bad_beta_minus_file.readlines()
bad_beta_minus_file.close()

bad_beta_plus_file = open("bad_beta_plus", "r")
bad_beta_plus = bad_beta_plus_file.readlines()
bad_beta_plus_file.close()

# remove newlines
for i in range(len(bad_beta_minus)):
  bad_beta_minus[i] = bad_beta_minus[i].strip()
for i in range(len(bad_beta_plus)):
  bad_beta_plus[i] = bad_beta_plus[i].strip()

# get the mo03 half-lives for beta minus
mo03_Z, mo03_N, mo03_t12, mo03_p0, mo03_p1, mo03_p2, mo03_p3 = \
    np.loadtxt("mo03_beta_minus.dat", unpack=True)

# get the bqa+ half-lives for beta plus
bqap_Z, bqap_N, bqap_t12 = np.loadtxt("bqa+_beta_plus.dat", unpack=True)

# open input and output files
fin = open(old_reaclib, "r")
fout = open(fixed_reaclib, "w")

num_changed_beta_minus = 0
num_changed_beta_plus = 0

# loop over reaclib entries
while(True):
  line = fin.readline()

  if (len(line) == 0):
    # end of file, we're done
    break

  # read lines that make up entry
  line0 = fin.readline()
  line1 = fin.readline()
  line2 = fin.readline()

  # beta decays are in chapter 1
  if ((line == "1\n") or (line == "2\n")):
    # this is an entry we potentially have to change
    if ((line0[43:48] == 'wc12w')
        and (line1[13:-1].lower().strip() == "0.000000e+00 0.000000e+00 0.000000e+00")
        and (line2[0:-1].lower().strip() == "0.000000e+00 0.000000e+00 0.000000e+00")):
      # this is a weak reaction with a constant rate, i.e. what we want
      if ((line0[5:10].strip() in bad_beta_minus)
          or ((line0[5:10].strip() in bad_beta_plus) and (line == "1\n"))):
        # this is a bad beta minus decay

        # get the reactant and product information
        numN = 0
        rIdx = nuclib.NuclideIdsVsNames()[line0[5:10].strip()]
        pIdx = nuclib.NuclideIdsVsNames()[line0[10:15].strip()]

        if (line == "2\n"):
          # one of the products is a neutron
          numN = 1
          if (line0[10:15].strip() == "n"):
            pIdx = nuclib.NuclideIdsVsNames()[line0[15:20].strip()]
          elif (line0[15:20].strip() != "n"):
            print "bad beta minus decay with two products but no neutron: %s" % line0[0:-1]
            sys.exit(1)

        rZ = nuclib.Zs()[rIdx]
        rA = nuclib.As()[rIdx]
        pZ = nuclib.Zs()[pIdx]
        pA = nuclib.As()[pIdx]

        new_rate = 0
        new_label = ""

        if (line0[5:10].strip() in bad_beta_minus):
          if ((rA != pA + numN) or (rZ != pZ - 1)):
            print "bad beta minus decay is not a beta minus decay: %s" % line0[0:-1]
            sys.exit(1)

          new_label = "mo03"
          mo03_idx = (mo03_Z == rZ) & (mo03_N == (rA - rZ))

          # total beta minus half-life
          t12 = mo03_t12[mo03_idx]

          # total beta minus decay rate
          total_rate = np.log(2.0) / t12

          # beta minus decay rate without any emitted neutrons
          new_rate = total_rate * mo03_p0[mo03_idx]
          if (line == "2\n"):
            # beta minus decay rate with 1 neutron emitted
            new_rate = total_rate * mo03_p1[mo03_idx]

          num_changed_beta_minus = num_changed_beta_minus + 1

        elif (line0[5:10].strip() in bad_beta_plus):
          if ((rA != pA + numN) or (rZ != pZ + 1)):
            print "bad beta plus decay is not a beta plus decay: %s" % line0[0:-1]
            sys.exit(1)

          new_label = "bqa+"
          bqap_idx = (bqap_Z == rZ) & (bqap_N == (rA - rZ))
          t12 = bqap_t12[bqap_idx]
          new_rate = np.log(2.0) / t12

          num_changed_beta_plus = num_changed_beta_plus + 1


        # change reaction label
        new_line0 = line0[:43] + new_label + line0[47:]
        line0 = new_line0

        # change reaction rate
        coefficient = np.log(new_rate)
        new_line1 = ("% 13.6e" % coefficient) + line1[13:]
        line1 = new_line1

        #print "%s: '% 13.6e'" % (line0[5:10], coefficient)

  # write the potentially modified entry
  fout.write(line)
  fout.write(line0)
  fout.write(line1)
  fout.write(line2)

# close input and output files
fin.close()
fout.close()

print "changed %i beta minus and %i beta plus entries" % (num_changed_beta_minus, num_changed_beta_plus)
