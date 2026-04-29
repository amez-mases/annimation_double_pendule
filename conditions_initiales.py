import numpy as np
import plotly.graph_objects as go

# constantes
g = 9.81
L1 = 1.0
L2 = 1.0
m1 = 1.0
m2 = 1.0


def double_pendule(y):  # les dérivées

    theta1, omega1, theta2, omega2 = y
    delta = theta2 - theta1

    denominateur1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    denominateur2 = (L2 / L1) * denominateur1

    domega1 = (
        m2 * L1 * omega1**2 * np.sin(delta) * np.cos(delta)
        + m2 * g * np.sin(theta2) * np.cos(delta)
        + m2 * L2 * omega2**2 * np.sin(delta)
        - (m1 + m2) * g * np.sin(theta1)
    ) / denominateur1

    domega2 = (
        -m2 * L2 * omega2**2 * np.sin(delta) * np.cos(delta)
        + (m1 + m2) * g * np.sin(theta1) * np.cos(delta)
        - (m1 + m2) * L1 * omega1**2 * np.sin(delta)
        - (m1 + m2) * g * np.sin(theta2)
    ) / denominateur2

    return np.array([omega1, domega1, omega2, domega2])


def simulation_double_pendule(dt=0.001, tmax=20):

    n = int(tmax / dt)

    thetas1 = np.zeros((n, 4))
    thetas2 = np.zeros((n, 4))

    # conditions initiales (presque identiques)
    thetas1[0] = [np.pi/2, 0, np.pi/2, 0]
    thetas2[0] = [np.pi/2 + 0.0001, 0, np.pi/2, 0]

    # intégration RK4
    for i in range(n - 1):

        k1 = double_pendule(thetas1[i])
        k2 = double_pendule(thetas1[i] + dt * k1 / 2)
        k3 = double_pendule(thetas1[i] + dt * k2 / 2)
        k4 = double_pendule(thetas1[i] + dt * k3)
        thetas1[i + 1] = thetas1[i] + dt * (k1 + 2*k2 + 2*k3 + k4) / 6

        k1 = double_pendule(thetas2[i])
        k2 = double_pendule(thetas2[i] + dt * k1 / 2)
        k3 = double_pendule(thetas2[i] + dt * k2 / 2)
        k4 = double_pendule(thetas2[i] + dt * k3)
        thetas2[i + 1] = thetas2[i] + dt * (k1 + 2*k2 + 2*k3 + k4) / 6

    return thetas1, thetas2


def positions(thetas):

    theta1 = thetas[:, 0]
    theta2 = thetas[:, 2]

    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)

    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)

    return x1, y1, x2, y2


def animer_double_pendule(x1a, y1a, x2a, y2a, x1b, y1b, x2b, y2b, step=1):

    # réduction du nombre de points pour accélérer l'animation
    x1a = x1a[::step]
    y1a = y1a[::step]
    x2a = x2a[::step]
    y2a = y2a[::step]

    x1b = x1b[::step]
    y1b = y1b[::step]
    x2b = x2b[::step]
    y2b = y2b[::step]

    xm = min(np.min(x2a), np.min(x2b), np.min(x1a), np.min(x1b)) - 0.3
    xM = max(np.max(x2a), np.max(x2b), np.max(x1a), np.max(x1b)) + 0.3
    ym = min(np.min(y2a), np.min(y2b), np.min(y1a), np.min(y1b)) - 0.3
    yM = max(np.max(y2a), np.max(y2b), np.max(y1a), np.max(y1b)) + 0.3

    n = len(x2a)

    fig = go.Figure(
        data=[
            # bras du pendule A
            go.Scatter(
                x=[0, x1a[0], x2a[0]],
                y=[0, y1a[0], y2a[0]],
                mode="lines+markers",
                line=dict(width=3, color="blue"),
                marker=dict(size=8, color="blue"),
                name="Pendule A"
            ),

            # bras du pendule B
            go.Scatter(
                x=[0, x1b[0], x2b[0]],
                y=[0, y1b[0], y2b[0]],
                mode="lines+markers",
                line=dict(width=3, color="red"),
                marker=dict(size=8, color="red"),
                name="Pendule B"
            ),

            # trajectoire progressive du pendule A
            go.Scatter(
                x=[x2a[0]],
                y=[y2a[0]],
                mode="lines",
                line=dict(width=2, color="blue"),
                name="Trajectoire du pendule A"
            ),

            # trajectoire progressive du pendule B
            go.Scatter(
                x=[x2b[0]],
                y=[y2b[0]],
                mode="lines",
                line=dict(width=2, color="red"),
                name="Trajectoire du pendule B"
            )
        ]
    )

    fig.update_layout(
        width=700,
        height=600,
        xaxis=dict(range=[xm, xM], autorange=False, zeroline=False),
        yaxis=dict(range=[ym, yM], autorange=False, zeroline=False, scaleanchor="x", scaleratio=1),
        title_text="Sensibilité aux conditions initiales – Double pendule",
        title_x=0.5,
        updatemenus=[
            dict(
                type="buttons",
                buttons=[
                    dict(
                        label="Play",
                        method="animate",
                        args=[
                            None,
                            {
                                "frame": {"duration": 20, "redraw": True},
                                "fromcurrent": True,
                                "transition": {"duration": 0}
                            }
                        ]
                    )
                ]
            )
        ]
    )

    fig.frames = [
        go.Frame(
            data=[
                # bras A à l'instant i
                go.Scatter(
                    x=[0, x1a[i], x2a[i]],
                    y=[0, y1a[i], y2a[i]]
                ),

                # bras B à l'instant i
                go.Scatter(
                    x=[0, x1b[i], x2b[i]],
                    y=[0, y1b[i], y2b[i]]
                ),

                # trajectoire A jusqu'à i
                go.Scatter(
                    x=x2a[:i+1],
                    y=y2a[:i+1]
                ),

                # trajectoire B jusqu'à i
                go.Scatter(
                    x=x2b[:i+1],
                    y=y2b[:i+1]
                )
            ],
            traces=[0, 1, 2, 3]
        )
        for i in range(n)
    ]

    fig.show()

# simulation
thetaA, thetaB = simulation_double_pendule()

x1a, y1a, x2a, y2a = positions(thetaA)
x1b, y1b, x2b, y2b = positions(thetaB)

# animation
animer_double_pendule(x1a, y1a, x2a, y2a, x1b, y1b, x2b, y2b, step=50)