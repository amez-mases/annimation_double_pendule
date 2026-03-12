import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# constantes
g = 9.81
L1 = 1.0
L2 = 1.0
m1 = 1.0
m2 = 1.0


def double_pendule(y): #les dérivées

    theta1, omega1, theta2, omega2 = y

    delta = theta2 - theta1

    denominateur1 = (m1 + m2)*L1 - m2*L1*np.cos(delta)**2 # équations du mouvement
    denominateur2 = (L2/L1)*denominateur1

    domega1 = (
        m2*L1*omega1**2*np.sin(delta)*np.cos(delta)
        + m2*g*np.sin(theta2)*np.cos(delta)
        + m2*L2*omega2**2*np.sin(delta)
        - (m1 + m2)*g*np.sin(theta1)
    ) / denominateur1 # accélération angulaire du premier bras

    domega2 = (
        -m2*L2*omega2**2*np.sin(delta)*np.cos(delta)
        + (m1 + m2)*g*np.sin(theta1)*np.cos(delta)
        - (m1 + m2)*L1*omega1**2*np.sin(delta)
        - (m1 + m2)*g*np.sin(theta2)
    ) / denominateur2 #accélération angulaire du deuxième bras

    return np.array([omega1, domega1, omega2, domega2])


def simulation_double_pendule(dt=0.01, tmax=20):

    n = int(tmax/dt)

    thetas1 = np.zeros((n,4))
    thetas2 = np.zeros((n,4))

    # conditions initiales (presque identiques)
    thetas1[0] = [np.pi/2, 0, np.pi/2, 0]
    thetas2[0] = [np.pi/2 + 0.01, 0, np.pi/2, 0]

    # intégration RK4
    for i in range(n-1):

        k1 = double_pendule(thetas1[i])
        k2 = double_pendule(thetas1[i] + dt*k1/2)
        k3 = double_pendule(thetas1[i] + dt*k2/2)
        k4 = double_pendule(thetas1[i] + dt*k3)

        thetas1[i+1] = thetas1[i] + dt*(k1 + 2*k2 + 2*k3 + k4)/6

        k1 = double_pendule(thetas2[i])
        k2 = double_pendule(thetas2[i] + dt*k1/2)
        k3 = double_pendule(thetas2[i] + dt*k2/2)
        k4 = double_pendule(thetas2[i] + dt*k3)

        thetas2[i+1] = thetas2[i] + dt*(k1 + 2*k2 + 2*k3 + k4)/6

    return thetas1, thetas2


def positions(thetas):

    theta1 = thetas[:,0]
    theta2 = thetas[:,2]

    x1 = L1*np.sin(theta1)
    y1 = -L1*np.cos(theta1)

    x2 = x1 + L2*np.sin(theta2)
    y2 = y1 - L2*np.cos(theta2)

    return x1,y1,x2,y2


# simulation
thetaA, thetaB = simulation_double_pendule()

x1a,y1a,x2a,y2a = positions(thetaA)
x1b,y1b,x2b,y2b = positions(thetaB)

pas = len(x2a)

# figure
fig,ax = plt.subplots()
ax.set_xlim(-2.2,2.2)
ax.set_ylim(-2.2,2.2)
ax.set_aspect("equal")
ax.set_title("Chaos double pendule")

ligne1, = ax.plot([],[],"o-",color="blue")
ligne2, = ax.plot([],[],"o-",color="red")

trace1, = ax.plot([],[],color="blue",linewidth=1)
trace2, = ax.plot([],[],color="red",linewidth=1)

trajectoire1x,trajectoire1y = [],[]
trajectoire2x,trajectoire2y = [],[]


def mise_a_jour(i):

    ligne1.set_data([0,x1a[i],x2a[i]],[0,y1a[i],y2a[i]])
    ligne2.set_data([0,x1b[i],x2b[i]],[0,y1b[i],y2b[i]])

    trajectoire1x.append(x2a[i])
    trajectoire1y.append(y2a[i])

    trajectoire2x.append(x2b[i])
    trajectoire2y.append(y2b[i])

    trace1.set_data(trajectoire1x,trajectoire1y)
    trace2.set_data(trajectoire2x,trajectoire2y)

    return ligne1,ligne2,trace1,trace2


ani = FuncAnimation(fig,mise_a_jour,frames=pas,interval=20)

plt.show()