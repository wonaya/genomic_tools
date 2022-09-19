import matplotlib.pyplot as plt
import matplotlib.dates as  mdates
import math

an = []
cs = []
sn = []
for i in range(0,101):
  a = 2 * math.pi * i / 100.0
  an.append(1290000000 + (i/100.0)*0030000000)
  cs.append(math.cos(a))
  sn.append(math.sin(a))

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(an, cs)
ax.plot(an, sn)
plt.show()
