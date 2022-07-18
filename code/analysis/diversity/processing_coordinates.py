#!/usr/bin/python

import sys

diversity = open(sys.argv[1], "r")

header = []
pi_dict = {}
for line in diversity:
    if line.startswith("Scaffold"):
         header.append(line.rstrip("\n"))
    else:
        line = line.strip().split()
        key = line[0]+":"+line[1]+":"+line[2]
        value = line[3:]
        if not key in pi_dict:
            pi_dict[key] = [value]
        else:
            pi_dict[key].append(value)

# print(pi_dict)

print("\t".join(header))
for key in pi_dict:
    if key.split(":")[0] == "superscaffold36":
        if len(pi_dict[key]) == 2:
            if pi_dict[key][0][1] == "NA" and pi_dict[key][1][1] == "NA":
                print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+"\t".join(pi_dict[key][0]))
            elif pi_dict[key][0][1] != "NA" and pi_dict[key][1][1] == "NA":
                print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+"\t".join(pi_dict[key][0]))
            elif pi_dict[key][0][1] == "NA" and pi_dict[key][1][1] != "NA":
                print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+"\t".join(pi_dict[key][1]))
            elif pi_dict[key][0][1] != "NA" and pi_dict[key][1][1] != "NA":
                pi = str(float(pi_dict[key][0][1]) + float(pi_dict[key][1][1]))
                theta = str(float(pi_dict[key][0][2]) + float(pi_dict[key][1][2]))
                Td = str(float(pi_dict[key][0][3]) + float(pi_dict[key][1][3]))
                pi_resamp_mean = str(float(pi_dict[key][0][5]) + float(pi_dict[key][1][5]))
                pi_resamp_low = str(float(pi_dict[key][0][6]) + float(pi_dict[key][1][6]))
                pi_resamp_high = str(float(pi_dict[key][0][7]) + float(pi_dict[key][1][7]))

                print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+pi_dict[key][0][0]+"\t"+pi+"\t"+theta+"\t"+Td+"\t"+pi_dict[key][0][4]+"\t"+pi_resamp_mean+"\t"+pi_resamp_low+"\t"+pi_resamp_high)
            else:
                raise Exception("Error!")
    else:
        if len(pi_dict[key]) == 1:
            print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+"\t".join(pi_dict[key][0]))
        else:
            raise Exception("Error!")

