import math
import subprocess

## INPUTS ##
little_name="mtb"
input_file = "bands.dat"
#
nbnd = 9
exclude_first = 0
exclude_last = 0

## OUTPUTS ##
bands_file = "bands_parsed_{}.dat".format(little_name)
kpts_file = "kpt_parsed_{}.dat".format(little_name)
merged_file = "output_{}.dat".format(little_name)



## ACTUAL CODE ##

# 1) Extract k-coordinates and list of KS eigenvalues from input_file
kpts = []
eigs = []
with open(input_file, "r") as f:
    for line in f:
        parts = line.split()
        if not parts:
            continue

        values = [float(x) for x in parts]
	
        k = values[:3]
        ev = values[3:3+nbnd]

        kpts.append(k)
        eigs.append(ev)


# 2) Calculates the relative k-distance in reciprocal space + normalization
d = [0.0]

for i in range(1, len(kpts)):
    dx = kpts[i][0] - kpts[i-1][0]
    dy = kpts[i][1] - kpts[i-1][1]
    dz = kpts[i][2] - kpts[i-1][2]
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
    d.append(d[-1] + dist)

for i in range(1, len(kpts)):
    d[i] = d[i] / d[-1]


# 4) Writes data in kpts_file
with open(kpts_file, "w") as f:
    for row in kpts:
        line = "".join(" %12.6f" % x for x in row)
        f.write(line + "\n")


# 3) Writes data in bands_file (for Gnuplot plotting):  Norm_k | Band 1 | Band 2 | ...
with open(bands_file, "w") as f:
    for di, ev in zip(d, eigs):
        line = "%.8f" % di
        for v in ev:
            line += " %12.6f" % v
        f.write(line + "\n")


# 4) Writes data in merged_file (if needed):  | kx | ky | kz | Norm_k | Band 1 | Band 2 | ...

cmd2 = (r'''#!/bin/bash

awk -v Nbands="{nbnd}" '
NR==FNR {{ a[NR] = $1 " " $2 " " $3; next }}
{{
    split(a[FNR], k, " ")
    printf "%9.4f%9.4f%9.4f    ", k[1], k[2], k[3]

    for(i=1;i<=1;i++)
        printf "%8.6f%s", $i, (i<NF ? " " : "")

    for(i=2+{exf};i<=Nbands+1-{exl};i++)
        printf "%11.6f%s", $i, (i<NF ? " " : "")

    print ""
}}
' {out1} {out2} > {out3}
''').format(nbnd=nbnd,exf=exclude_first,exl=exclude_last,out1=kpts_file,out2=bands_file,out3=merged_file)

subprocess.call(cmd2, shell=True, executable="/bin/bash")


# 5) All done.
print("Output data written in:   {} , {}".format(bands_file, merged_file))
