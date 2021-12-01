import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as s

with open('MeteoWindSpeed.txt', 'r') as f:
    data = f.readlines()


wind_speed = np.arange(0, 10.2, 0.3)

nbre_cases = np.zeros(len(wind_speed))
for line in data:
    speed = float(line)
    for i in range(len(wind_speed)):
        if speed < wind_speed[i]:
            nbre_cases[i-1] += 1
            break


def weib(s,n,a):

    return (a / n) * (s / n)**(a - 1) * np.exp(-(s / n)**a)


plt.rc('axes', linewidth=2)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

ax = plt.gca()

#2 ways that does not work...

#popt, pcov = curve_fit(weib, wind_speed, nbre_cases)
#plt.plot(wind_speed, weib(wind_speed, *popt), 'g-', label='Weibull distribution')

new_data = []
for line in data:
    new_data.append(float(line))
#(loc, scale) = s.exponweib.fit(nbre_cases)
#plt.plot(new_data, s.exponweib.pdf(new_data, *s.exponweib.fit(new_data, 1, 1, scale=0.2, loc=0)))

plt.plot(wind_speed, nbre_cases, 'o')


for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
    tick.label1.set_fontweight('bold')

plt.xlabel('Wind speed [m/s]', fontsize=16, fontweight='bold')
plt.ylabel('Number of cases [/]', fontsize=16, fontweight='bold')

plt.legend()
plt.show()
