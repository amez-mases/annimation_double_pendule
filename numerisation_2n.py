import plotly.graph_objects as go
import numpy as np

def animer_double_pendule_2D(thetas, L1=1, L2=1, step=1):

    theta1 = thetas[:,0][::step]
    theta2 = thetas[:,2][::step]

    # Positions des masses
    x1 = L1*np.sin(theta1)
    y1 = -L1*np.cos(theta1)

    x2 = x1 + L2*np.sin(theta2)
    y2 = y1 - L2*np.cos(theta2)

    fig = go.Figure(
        data=[
            # Tiges
            go.Scatter(
                x=[0, x1[0], x2[0]],
                y=[0, y1[0], y2[0]],
                mode="lines",
                line=dict(width=4, color="black")
            ),

            # Masse 2 (extrémité libre)
            go.Scatter(
                x=[x2[0]],
                y=[y2[0]],
                mode="markers",
                marker=dict(size=10, color="red")
            ),

            # Trajectoire
            go.Scatter(
                x=[x2[0]],
                y=[y2[0]],
                mode="lines",
                line=dict(width=2, color="blue")
            )
        ]
    )

    frames = []

    for k in range(1, len(x1)):
        frames.append(
            go.Frame(
                data=[
                    # Tiges
                    go.Scatter(
                        x=[0, x1[k], x2[k]],
                        y=[0, y1[k], y2[k]]
                    ),

                    # extrémité libre
                    go.Scatter(
                        x=[x2[k]],
                        y=[y2[k]]
                    ),

                    # trajectoire
                    go.Scatter(
                        x=x2[:k],
                        y=y2[:k]
                    )
                ]
            )
        )

    fig.frames = frames

    fig.update_layout(
        xaxis=dict(range=[-2.2, 2.2], zeroline=False),
        yaxis=dict(range=[-2.2, 2.2], zeroline=False, scaleanchor="x", scaleratio=1),
        width=650,
        height=650,
        title="Double Pendule",
        title_x=0.5,
        updatemenus=[dict(
            type="buttons",
            buttons=[dict(
                label="Play",
                method="animate",
                args=[None,
                      {"frame": {"duration": 20, "redraw": True},
                       "fromcurrent": True}]
            )]
        )]
    )

    fig.show()


# Paramètres physiques
g = 9.81
L1 = 1.0
L2 = 1.0
m1 = 1.0
m2 = 1.0

def double_pendule(thetas):
    theta1, omega1, theta2, omega2 = thetas

    delta = theta2 - theta1

    den1 = (m1 + m2)*L1 - m2*L1*np.cos(delta)**2
    den2 = (L2/L1)*den1

    domega1 = (
        m2*L1*omega1**2*np.sin(delta)*np.cos(delta)
        + m2*g*np.sin(theta2)*np.cos(delta)
        + m2*L2*omega2**2*np.sin(delta)
        - (m1 + m2)*g*np.sin(theta1)
    ) / den1

    domega2 = (
        -m2*L2*omega2**2*np.sin(delta)*np.cos(delta)
        + (m1 + m2)*g*np.sin(theta1)*np.cos(delta)
        - (m1 + m2)*L1*omega1**2*np.sin(delta)
        - (m1 + m2)*g*np.sin(theta2)
    ) / den2

    return np.array([omega1, domega1, omega2, domega2])

def simulation_double_pendule(dt=0.01, tmax=20):

    n = int(tmax/dt)
    thetas = np.zeros((n,4))

    # Conditions initiales
    thetas[0] = [np.pi/2, 0, np.pi/2 + 0.01, 0]

    # Intégration RK4
    for i in range(n-1):
        k1 = double_pendule(thetas[i])
        k2 = double_pendule(thetas[i] + dt*k1/2)
        k3 = double_pendule(thetas[i] + dt*k2/2)
        k4 = double_pendule(thetas[i] + dt*k3)
        thetas[i+1] = thetas[i] + dt*(k1 + 2*k2 + 2*k3 + k4)/6

    return thetas

thetas = simulation_double_pendule()
animer_double_pendule_2D(thetas, step=2)

