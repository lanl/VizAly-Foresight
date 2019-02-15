#! /usr/bin/env python

input_file = "old.csv"
vinput_file = "v.csv"
output_file = "data.csv"

with open(input_file, "r") as fp:
    lines = fp.readlines()

with open(vinput_file, "r") as fp:
    vlines = fp.readlines()

for i, line in enumerate(lines):

    line = line.rstrip("\n")

    if i == 0:
        line += ",FILE_1,FILE_2,FILE_3\n"
    else:
        anx = ",,,\n"
        for vline in vlines:
            oname = line.split(",")[0]
            vname = vline.split(",")[0]
            name = oname.split("_")[0]
            args = "__" + "__".join(oname.split("__")[1:])
            if vname == name + args:
                anx = "," + vline.rstrip("\n")[len(vname) + 1:] + "\n"
        line += anx

    lines[i] = line

with open(output_file, "w") as fp:
    for line in lines:
        fp.write(line)
