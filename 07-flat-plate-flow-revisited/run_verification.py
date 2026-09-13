"""Run geometric and manufactured pressure-correction checks, without a GUI."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from staggered_grid import Grid, make_grid
from pressure_projection import project


def check_grid(grid):
    ny, nx = grid.shape
    # Affine vector field: exact divergence = 2-3 = -1, including boundary cells.
    u = 2*grid.xf[None, :]+np.zeros((ny, 1))
    v = -3*grid.yf[:, None]+np.zeros((1, nx))
    divergence_error = np.max(np.abs(grid.divergence(u, v)+1))
    assert divergence_error < 1e-11
    # A face-streamfunction field is discretely divergence-free by construction.
    xx, yy = np.meshgrid(grid.xf, grid.yf)
    psi = np.sin(np.pi*xx/grid.xf[-1])*np.sin(np.pi*yy/grid.yf[-1])
    psi[[0, -1], :] = 0
    psi[:, [0, -1]] = 0
    base_u = np.diff(psi, axis=0)/grid.dy[:, None]
    base_v = -np.diff(psi, axis=1)/grid.dx
    # Add a known pressure gradient with homogeneous normal boundary gradient.
    xx, yy = np.meshgrid(grid.xc, grid.yc)
    p = np.cos(np.pi*xx/grid.xf[-1])*np.cos(np.pi*yy/grid.yf[-1])
    p -= p[0, 0]
    u, v = base_u.copy(), base_v.copy()
    dt, rho = 0.03, 1.205
    u[:, 1:-1] += dt/rho*np.diff(p, axis=1)/np.diff(grid.xc)
    v[1:-1] += dt/rho*np.diff(p, axis=0)/np.diff(grid.yc)[:, None]
    unew, vnew, recovered = project(grid, u, v, dt, rho)
    before = np.max(np.abs(grid.divergence(u, v)))
    after = np.max(np.abs(grid.divergence(unew, vnew)))
    pressure_error = np.max(np.abs(recovered-p))
    velocity_error = max(np.max(np.abs(unew-base_u)), np.max(np.abs(vnew-base_v)))
    assert after < 1e-7 and after/before < 1e-9
    assert pressure_error < 1e-9 and velocity_error < 1e-9
    np.testing.assert_array_equal(unew[:, [0, -1]], u[:, [0, -1]])
    np.testing.assert_array_equal(vnew[[0, -1]], v[[0, -1]])
    # Net-inflow errors must be rejected, not hidden by pinning a pressure cell.
    incompatible = u.copy()
    incompatible[:, 0] += 1
    try:
        project(grid, incompatible, v)
    except ValueError:
        pass
    else:
        raise AssertionError('Incompatible boundary flux was accepted.')
    return dict(nx=nx, ny=ny, first_dy_m=grid.dy[0],
                affine_divergence_error=divergence_error,
                divergence_before=before, divergence_after=after,
                pressure_recovery_error=pressure_error,
                velocity_recovery_error=velocity_error)


def main():
    out = Path(__file__).resolve().parent/'results'
    out.mkdir(exist_ok=True)
    grids = [('uniform', make_grid(24, 16)),
             ('wall_clustered', make_grid(24, 16, stretching=1.5)),
             ('wall_clustered_fine', make_grid(48, 32, stretching=1.5)),
             ('both_axes_nonuniform', Grid(np.linspace(0, 1, 25)**1.3,
                                          make_grid(24, 16, stretching=2).yf))]
    rows = [dict(case=name, **check_grid(grid)) for name, grid in grids]
    with (out/'operator_verification.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    g = make_grid(8, 6, stretching=1.5)
    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
    for x in g.xf:
        ax.axvline(x, color='0.8', linewidth=0.7)
    for y in g.yf:
        ax.axhline(y, color='0.8', linewidth=0.7)
    for xs, ys, marker, label in [(g.xc, g.yc, 'o', 'p, T: cell centers'),
                                  (g.xf, g.yc, '>', 'u: vertical faces'),
                                  (g.xc, g.yf, '^', 'v: horizontal faces')]:
        xx, yy = np.meshgrid(xs, ys)
        ax.scatter(xx, yy, s=20, marker=marker, label=label)
    ax.set(xlim=(0, 1), ylim=(0, .1), xlabel='x (m)', ylabel='y (m)',
           title='Explicit staggered locations on a wall-clustered mesh (axes not to scale)')
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.17), ncol=3)
    fig.savefig(out/'staggered_grid.png', dpi=180)
    plt.close(fig)
    for row in rows:
        print(f"{row['case']}: divergence {row['divergence_before']:.3e} -> "
              f"{row['divergence_after']:.3e}, pressure error "
              f"{row['pressure_recovery_error']:.3e}")
    print('Geometry and closed-boundary pressure-projection checks passed.')


if __name__ == '__main__':
    main()
