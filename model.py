import numpy as np

# invariants
R = 63 # Radius
B = 3 # Number of blades

# conditions initiales
Uinf = 11.4 # Mean Wind Speed [m/s]
beta = 0 # Local element pitch angle [°]
rpm = 12.1  # Angular velocity [rpm]
rho = 1.225 # p/(287.1*T) # Air density [kg/m^3]

omega = 2*np.pi/60*rpm

c = 3.5 # Chord length , ça doit varier normalement, j'ai pris une valeur au hazard
dr = 1
r = np.arange(1.5,63.01,dr)
for n in range(len(r)):
    fichier = open('DU25_A17.dat', 'r')
    i = -1
    for ligne in fichier :
        i += 1
        if i < 61:
            continue
        alpha = float(ligne.split()[0])   # Local angle of attack
        Cl = float(ligne.split()[1])      # Lift coefficient
        Cd = float(ligne.split()[2])      # Drag coefficient
        alpha = alpha/360*2*np.pi
        print(r[n])
        phi = alpha + beta         # Local flow angle
        cst = 4*np.pi*r[n]
        dst = 1/2*B*c/((np.sin(phi))**2)*(Cl*np.cos(phi)+Cd*np.sin(phi))
        a = dst/(cst+dst)   # Axial Induction Factor
        aprime = Uinf*(1-a)/(np.tan(phi)*omega*r[n]) - 1   # rotational induction factor
        Vres = (1-a)**2 *Uinf*Uinf/((np.sin(phi))**2)     # Velocity seen by the blades
        dQ1 = B/2 *rho*Vres*Vres*(Cl*np.sin(phi)-Cd*np.cos(phi))*c*r[n]*dr
        dQ2 = 4*np.pi*r[n]**3*rho*Uinf*(1-a)*aprime*omega*dr
        print(dQ1-dQ2)
