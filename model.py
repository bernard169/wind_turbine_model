import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# invariants
R = 63  # Radius
B = 3   # Number of blades

# conditions initiales
Uinf = np.arange(3, 26, 1)   # Mean Wind Speed [m/s]
beta = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 7, 8.5, 10.5, 12, 13.5, 15, 16, 17, 18, 20, 21, 22, 23])    # Local element pitch angle [°]
rpm = np.array([7, 7.2, 7.5, 8, 8.5, 9.1, 10.2, 11.5, 11.9, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1, 12.1])  # Angular velocity [rpm]
rho = 1.225 # p/(287.1*T) # Air density [kg/m^3]


def chord(r):
    if r < 14.05:
        c = 0.1143*r + 3.2143
    elif r < 58.9:
        c = -0.0597*r + 5.6584
    else :
        c = -0.244*r + 16.459
    return c


dr = 1
r = np.arange(1.5, 63.1, dr)
Q = 0
P = np.zeros(len(Uinf))

a8 = np.zeros(len(r))
aprime8 = np.zeros(len(r))
a22 = np.zeros(len(r))
aprime22 = np.zeros(len(r))

data = []
with open('DU25_A17.dat', 'r') as f : 
    data = f.readlines()

for h in range(len(Uinf)):

    omega = 2 * np.pi / 60 * rpm[h]
    for n in range(len(r)):
        i = -1
        c = chord(r[n])
        alpha = []
        Cl = []
        Cd = []
        for line in data:
            i += 1
            if i < 60:
                continue

            alpha.append(float(line.split()[0])/360*2*np.pi)   # Local angle of attack
            Cl.append(float(line.split()[1]))  # Lift coefficient
            Cd.append(float(line.split()[2]))  # Drag coefficient
        i = -1
        mi = 1000

        for line in data:
            print(line)
            i += 1
            if i < 80:
                for j in range(0, 1000):
                    alpha1 = alpha[i] + j * (alpha[i + 1] - alpha[i]) / 1000
                    Cl1 = Cl[i] + j * (Cl[i + 1] - Cl[i]) / 1000
                    Cd1 = Cd[i] + j * (Cd[i + 1] - Cd[i]) / 1000
                    phi = alpha1 + beta[h] / 360 * 2 * np.pi # Local flow angle
                    cst = 4 * np.pi * r[n]
                    dst = 1 / 2 * B * c / ((np.sin(phi)) ** 2) * (Cl1 * np.cos(phi) + Cd1 * np.sin(phi))
                    a = dst / (cst + dst)  # Axial Induction Factor
                    aprime = Uinf[h] * (1 - a) / (np.tan(phi) * omega * r[n]) - 1  # rotational induction factor
                    Vres = (1 - a) * Uinf[h] / np.sin(phi)
                    dQ1 = B / 2 * rho * Vres * Vres * (Cl1 * np.sin(phi) - Cd1 * np.cos(phi)) * c * r[n] * dr
                    dQ2 = 4 * np.pi * r[n] ** 3 * rho * Uinf[h] * (1 - a) * aprime * omega * dr
                    # print(dQ1)
                    if abs(dQ1 - dQ2) < mi and dQ1 > 0 and dQ2 > 0 and a < 0.5:
                        alphamin = alpha1
                        amin = a
                        mi = abs(dQ1 - dQ2)
                        dQ = dQ1
                        #print(dQ)
        if h == 7:
            a8[n] = a
            aprime8[n] = aprime
        if h == 21:
            a22[n] = a
            aprime22[n] = aprime

        # print(aphamin/2*360/np.pi)
        #print(mi)
        #print(dQ)
        Q += dQ

    #print(Q)
    P[h] = Q*omega
    print(P[h])
    
print(a8)
print(aprime8)
print(a22)
print(aprime22)

print(P)

"""
Plot Power curve
"""
plt.figure(1)

plt.rc('axes', linewidth=2)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

ax = plt.gca()

f = interp1d(Uinf, P)
xnew = np.linspace(3, 25, endpoint=True)

plt.plot(Uinf, P, 'o', xnew, f(xnew), '-')

#plt.scatter(wind_speeds, P, color="orange")

plt.xlim(0, 26)
plt.ylim(0, 6*1e6)

plt.grid(True, which='both')

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
    tick.label1.set_fontweight('bold')

plt.xlabel('Wind Speed [$m/s$]', fontsize=16, fontweight='bold')
plt.ylabel('Power [$W$]', fontsize=16, fontweight='bold')

#plt.legend()
plt.show()

"""
Plot Axial Induction Factors
"""
plt.figure(2)

plt.rc('axes', linewidth=2)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

ax = plt.gca()

#f = interp1d(Uinf, P)
#xnew = np.linspace(3, 25, endpoint=True)

plt.plot(r, a8, 'o')
plt.plot(r, a22, 'o')

#plt.scatter(wind_speeds, P, color="orange")

#plt.xlim(0, 26)
#plt.ylim(0, 6*1e6)

plt.grid(True, which='both')

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
    tick.label1.set_fontweight('bold')

plt.xlabel('Blade Radius [m]', fontsize=16, fontweight='bold')
plt.ylabel('Axial Induction Factor [/]', fontsize=16, fontweight='bold')

#plt.legend()
plt.show()

"""
Rotational Induction Factor
"""
plt.figure(3)

plt.rc('axes', linewidth=2)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

ax = plt.gca()

#f = interp1d(Uinf, P)
#xnew = np.linspace(3, 25, endpoint=True)

plt.plot(r, aprime8, 'o')
plt.plot(r, aprime22, 'o')

#plt.scatter(wind_speeds, P, color="orange")

#plt.xlim(0, 26)
#plt.ylim(0, 6*1e6)

plt.grid(True, which='both')

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
    tick.label1.set_fontweight('bold')

plt.xlabel('Blade Radius [m]', fontsize=16, fontweight='bold')
plt.ylabel('Rotational Induction Factor [/]', fontsize=16, fontweight='bold')

#plt.legend()
plt.show()