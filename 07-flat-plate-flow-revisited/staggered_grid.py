"""Explicit cell and face geometry; arrays are indexed [y, x]."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Grid:
    xf: np.ndarray
    yf: np.ndarray

    def __post_init__(self):
        for name in ('xf', 'yf'):
            a = np.array(getattr(self, name), dtype=float, copy=True)
            if a.ndim != 1 or len(a) < 3 or not np.all(np.isfinite(a)):
                raise ValueError('Faces must be finite one-dimensional arrays with >=2 cells.')
            if np.any(np.diff(a) <= 0):
                raise ValueError('Faces must be strictly increasing.')
            a.flags.writeable = False
            object.__setattr__(self, name, a)

    @property
    def xc(self):
        return (self.xf[:-1] + self.xf[1:])/2

    @property
    def yc(self):
        return (self.yf[:-1] + self.yf[1:])/2

    @property
    def dx(self):
        return np.diff(self.xf)

    @property
    def dy(self):
        return np.diff(self.yf)

    @property
    def shape(self):
        return len(self.yc), len(self.xc)

    @property
    def volume(self):
        return self.dy[:, None]*self.dx[None, :]

    def divergence(self, u, v):
        ny, nx = self.shape
        if u.shape != (ny, nx+1) or v.shape != (ny+1, nx):
            raise ValueError('u and v must live on x and y faces, respectively.')
        return np.diff(u, axis=1)/self.dx + np.diff(v, axis=0)/self.dy[:, None]


def make_grid(nx=40, ny=40, lx=1.0, ly=0.1, stretching=0.0):
    if not isinstance(nx, int) or not isinstance(ny, int) or min(nx, ny) < 2:
        raise ValueError('nx and ny must be integers >=2.')
    if not np.isfinite([lx, ly, stretching]).all() or min(lx, ly) <= 0:
        raise ValueError('Lengths must be positive and parameters finite.')
    if not 0 <= stretching <= 5:
        raise ValueError('stretching must lie between 0 and 5.')
    eta = np.linspace(0, 1, ny+1)
    yf = eta if stretching == 0 else (
        1-np.tanh(stretching*(1-eta))/np.tanh(stretching))
    return Grid(np.linspace(0, lx, nx+1), ly*yf)
