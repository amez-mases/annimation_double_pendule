import numpy as np
import plotly.graph_objects as go
from scipy.integrate import solve_ivp

m1 = 2496
m2 = 100
mu = m2 / (m1 + m2)
C = 3.0

# position et vitesse initiales de m3
x0 = 0.5
y0 = 0.0
z0 = 0.0
vx0 = 0.0
vz0 = 0.0

t_span = (0, 50)
nb_points = 5000
t_eval = np.linspace(t_span[0], t_span[1], nb_points)

# positions des corps primaires dans le repère tournant
m1_position = np.array([-mu, 0])
m2_position = np.array([1 - mu, 0])

def distances(x, y, z):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2 + z**2)
    return r1, r2

def potentiel_effectif(x, y, z):
    r1, r2 = distances(x, y, z)
    return (1 - mu) / r1 + mu / r2 + 0.5 * (x**2 + y**2)

def vitesse_initiale_avec_C(x, y, z, vx, C):
    Omega = potentiel_effectif(x, y, z)
    v_carre = 2 * Omega - C - vx**2

    if v_carre < 0:
        raise ValueError(
            "Impossible de calculer la vitesse initiale. "
            "Choisissez une valeur de C différente."
        )

    vy = np.sqrt(v_carre)
    return vy

vy0 = vitesse_initiale_avec_C(x0, y0, z0, vx0, C)

# équations différentielles du mouvement
def equations(t, etat):
    x, y, z, vx, vy, vz = etat

    r1, r2 = distances(x, y, z)

    ax = (
        x
        + 2 * vy
        - (1 - mu) * (x + mu) / r1**3
        - mu * (x - 1 + mu) / r2**3
    )

    ay = (
        y
        - 2 * vx
        - (1 - mu) * y / r1**3
        - mu * y / r2**3
    )

    az = (
        -(1 - mu) * z / r1**3
        - mu * z / r2**3
    )

    return [vx, vy, vz, ax, ay, az]

etat_initial = [x0, y0, z0, vx0, vy0, vz0]

solution = solve_ivp(
    equations,
    t_span,
    etat_initial,
    t_eval=t_eval,
    rtol=1e-10,
    atol=1e-12
)

def animer(x, y, step=1):
    x = x[::step]
    y = y[::step]

    xx = x
    yy = y

    # positions fixes des deux corps principaux
    x1p, y1p = m1_position[0], m1_position[1]
    x2p, y2p = m2_position[0], m2_position[1]

    # bornes du graphique en incluant aussi les corps fixes
    xm = min(np.min(xx), x1p, x2p) - 1.5
    xM = max(np.max(xx), x1p, x2p) + 1.5
    ym = min(np.min(yy), y1p, y2p) - 1.5
    yM = max(np.max(yy), y1p, y2p) + 1.5

    fig = go.Figure(
        data=[
            # trajectoire progressive du troisième corps
            go.Scatter(
                x=[xx[0]],
                y=[yy[0]],
                mode="lines",
                line=dict(width=2, color="blue"),
                name="Trajectoire du troisième corps"
            ),

            # troisième corps mobile
            go.Scatter(
                x=[xx[0]],
                y=[yy[0]],
                mode="markers",
                marker=dict(color="red", size=10),
                name="Troisième corps"
            ),

            # premier corps principal fixe
            go.Scatter(
                x=[x1p],
                y=[y1p],
                mode="markers",
                marker=dict(color="black", size=12),
                name="Corps principal 1"
            ),

            # deuxième corps principal fixe
            go.Scatter(
                x=[x2p],
                y=[y2p],
                mode="markers",
                marker=dict(color="green", size=11),
                name="Corps principal 2"
            )
        ]
    )

    fig.update_layout(
        width=700,
        height=550,
        xaxis=dict(range=[xm, xM], autorange=False, zeroline=False),
        yaxis=dict(
            range=[ym, yM],
            autorange=False,
            zeroline=False,
            scaleanchor="x",
            scaleratio=1
        ),
        title_text="Trajectoire du troisième corps dans le problème restreint à trois corps",
        title_x=0.5,
        updatemenus=[
            dict(
                type="buttons",
                buttons=[
                    dict(
                        args=[
                            None,
                            {
                                "frame": {"duration": 10, "redraw": True},
                                "fromcurrent": True,
                                "transition": {"duration": 10}
                            }
                        ],
                        label="Play",
                        method="animate",
                    )
                ]
            )
        ]
    )

    fig.update(
        frames=[
            go.Frame(
                data=[
                    # trajectoire jusqu'à la position actuelle
                    go.Scatter(
                        x=xx[:k+1],
                        y=yy[:k+1]
                    ),

                    # position actuelle du troisième corps
                    go.Scatter(
                        x=[xx[k]],
                        y=[yy[k]]
                    )
                ],
                traces=[0, 1]
            )
            for k in range(len(xx))
        ]
    )

    fig.show()

x = solution.y[0]
y = solution.y[1]

animer(x, y, 3)
