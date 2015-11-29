import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 8:
    print "input1 input2 input3 label1 label2 label3 title"
    sys.exit()

x = []
y = []
for line in open(sys.argv[1]):
    splitline = line.split()
    x.append(float(splitline[0]))
    y.append(float(splitline[2]))

plt.plot(x, y, '.r-', label=sys.argv[4])

x = []
y = []
for line in open(sys.argv[2]):
    splitline = line.split()
    x.append(float(splitline[0]))
    y.append(float(splitline[2]))

plt.plot(x, y, 'xb-', label=sys.argv[5])

x = []
y = []
for line in open(sys.argv[3]):
    splitline = line.split()
    x.append(float(splitline[0]))
    y.append(float(splitline[2]))

plt.plot(x, y, 'og-', label=sys.argv[6])

ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.title(sys.argv[7])
plt.xlabel("#Knoten")
plt.ylabel("Fehler (Energienorm)")


plt.legend(loc='upper right')
plt.show()



