import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 12:
    print "input1 input2 input3 input4 input5 label1 label2 label3 label4 label5 title"
    sys.exit()

x = []
y = []
for line in open(sys.argv[1]):
    splitline = line.split()
    x.append(float(splitline[0]))
    y.append(float(splitline[3]))

plt.plot(x, y, '.r-', label=sys.argv[6], linewidth=1.5)

x = []
y = []
for line in open(sys.argv[2]):
    splitline = line.split()
    x.append(float(splitline[0]))
    y.append(float(splitline[3]))

plt.plot(x, y, 'xb-', label=sys.argv[7], linewidth=1.5)

x = []
y = []
for line in open(sys.argv[3]):
    splitline = line.split()
    x.append(float(splitline[0]))
    y.append(float(splitline[3]))

plt.plot(x, y, 'og-', label=sys.argv[8], linewidth=1.5)

x = []
y = []
for line in open(sys.argv[4]):
    splitline = line.split()
    x.append(float(splitline[0]))
    y.append(float(splitline[3]))

plt.plot(x, y, 's-', color="purple", label=sys.argv[9], linewidth=1.5)

x = []
y = []
for line in open(sys.argv[5]):
    splitline = line.split()
    x.append(float(splitline[0]))
    y.append(float(splitline[3]))

plt.plot(x, y, 'p-', color="black", label=sys.argv[10], linewidth=1.5)
plt.grid(b=True, which='major', color='black', linestyle='-')

ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.title(sys.argv[11])
plt.xlabel("#Knoten")
plt.ylabel("Fehler (Energienorm)")


plt.legend(loc='best')
plt.savefig('test.png', bbox_inches='tight')



