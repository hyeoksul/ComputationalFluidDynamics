"""Manufactured closed-boundary test of compatible divergence and gradient.

This is a projection building block, not a SIMPLEC/SIMPLER flow solver.
Prescribed boundary normal velocities remain fixed. Compatible net flux is
required; one pressure value fixes the constant-pressure nullspace.
"""
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


def project(grid, u, v, dt=1.0, rho=1.205):
    if not np.isfinite([dt, rho]).all() or min(dt, rho) <= 0:
        raise ValueError('dt and rho must be positive and finite.')
    if not np.isfinite(u).all() or not np.isfinite(v).all():
        raise ValueError('Velocities must be finite.')
    imbalance = grid.divergence(u, v)*grid.volume
    if abs(imbalance.sum()) > 1e-12*max(np.abs(imbalance).sum(), 1e-14):
        raise ValueError('Prescribed boundary velocities have incompatible net flux.')
    ny, nx = grid.shape
    matrix = lil_matrix((nx*ny, nx*ny))
    # Symmetric integrated -Laplacian; each internal face is assembled once.
    for j in range(ny):
        for i in range(nx):
            row = j*nx+i
            neighbors = []
            if i+1 < nx:
                neighbors.append((row+1, grid.dy[j]/(grid.xc[i+1]-grid.xc[i])))
            if j+1 < ny:
                neighbors.append((row+nx, grid.dx[i]/(grid.yc[j+1]-grid.yc[j])))
            for col, conductance in neighbors:
                matrix[row, row] += conductance
                matrix[col, col] += conductance
                matrix[row, col] -= conductance
                matrix[col, row] -= conductance
    rhs = -rho/dt*imbalance.ravel()
    pressure = np.zeros(nx*ny)
    # Eliminate only the gauge degree of freedom; verify its continuity row too.
    pressure[1:] = spsolve(matrix.tocsr()[1:, 1:], rhs[1:])
    pressure = pressure.reshape(ny, nx)
    corrected_u, corrected_v = u.copy(), v.copy()
    corrected_u[:, 1:-1] -= dt/rho*np.diff(pressure, axis=1)/np.diff(grid.xc)
    corrected_v[1:-1, :] -= dt/rho*np.diff(pressure, axis=0)/np.diff(grid.yc)[:, None]
    return corrected_u, corrected_v, pressure
