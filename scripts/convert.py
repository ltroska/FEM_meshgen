import sys

if len(sys.argv) != 3:
    print "input output"
    sys.exit()

outfile = open(sys.argv[2], "w")

outfile.write("{")
for i, line in enumerate(open(sys.argv[1])):
    if i > 0:
        outfile.write(",")
    splitline = line.split()
    outfile.write("{")
    outfile.write(splitline[0]+","+splitline[2])
    outfile.write("}")

outfile.write("}")
outfile.close()

