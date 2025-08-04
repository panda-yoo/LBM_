import numpy as np

NX = 100
NY = 40
NO_Q = 9

buffer_nx = NX + 2
buffer_ny = NY + 2

INFLOW_VELOCITY_UX = 0.4
REYNOLDS_NUMBER = 100
CYLINDER_RADIUS = 8
CS_2 = 1.0 / 3.0

CENTER_CYLINDER = np.array([NX // 4, NY // 2])

w = np.array(
    [
        4.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
    ],
    dtype=np.float64
)

kinematic_viscosity = (INFLOW_VELOCITY_UX * CYLINDER_RADIUS) / REYNOLDS_NUMBER

OMEGA = 1.0 / (3.0 * kinematic_viscosity + 0.5)
TAU = (1 / OMEGA) - .5

cx = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
cy = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
