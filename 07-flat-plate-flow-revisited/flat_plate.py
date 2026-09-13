"""Steady staggered FVM: donor-cell momentum, SIMPLEC and energy transport.

All arrays use [y, x]. Open normal-velocity boundaries have prescribed gauge
pressure; their fluxes participate in the pressure-correction equation.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve
from staggered_grid import make_grid


@dataclass(frozen=True)
class Physics:
    rho: float = 1.205
    cp: float = 1005.0
    k: float = 0.0257
    pr: float = 0.713
    speed: float = 0.5
    cold: float = 300.0
    hot: float = 400.0
    heating_start: float = 0.25

    @property
    def mu(self):
        return self.pr*self.k/self.cp


def transport(xs, ys, widths, heights, fx, fy, diffusivity, boundaries):
    """Integrated conservative upwind operator on an orthogonal grid.

    boundaries: west/east/south/north -> (value, distance). A None value is
    zero diffusive flux with advective face value equal to the boundary cell.
    Dirichlet distance is measured from the unknown to the prescribed node.
    fx and fy are integrated advective flux coefficients oriented +x/+y.
    """
    ny, nx = len(ys), len(xs)
    rows, cols, data = [], [], []
    rhs = np.zeros((ny, nx))
    diagonal = np.zeros((ny, nx))
    for j in range(ny):
        for i in range(nx):
            r = j*nx+i
            faces = [(i-1, j, -fx[j,i], heights[j], 'west', i),
                     (i+1, j, fx[j,i+1], heights[j], 'east', nx-1-i),
                     (i, j-1, -fy[j,i], widths[i], 'south', j),
                     (i, j+1, fy[j+1,i], widths[i], 'north', ny-1-j)]
            for ni, nj, flux, area, side, interior in faces:
                if interior:
                    distance = abs(xs[ni]-xs[i]) if ni != i else abs(ys[nj]-ys[j])
                    diff = diffusivity*area/distance
                    diagonal[j,i] += diff+max(flux, 0)
                    rows.append(r); cols.append(nj*nx+ni)
                    data.append(-diff+min(flux, 0))
                else:
                    value, distance = boundaries[side]
                    if value is None:
                        diagonal[j,i] += flux
                    else:
                        if np.ndim(value):
                            value = value[j] if side in ('west','east') else value[i]
                        diff = diffusivity*area/distance
                        diagonal[j,i] += diff+max(flux, 0)
                        rhs[j,i] += (diff-min(flux, 0))*value
            rows.append(r); cols.append(r); data.append(diagonal[j,i])
    return coo_matrix((data,(rows,cols)), shape=(nx*ny,nx*ny)).tocsr(), rhs


def momentum(g, f, u, v, p, component):
    ny, nx = g.shape
    if component == 'u':
        # u unknowns include the pressure outlet; inlet u is prescribed.
        xs, ys = g.xf[1:], g.yc
        widths = np.diff(np.r_[g.xc, g.xf[-1]])
        heights = g.dy
        uf = np.c_[(u[:,:-1]+u[:,1:])/2, u[:,-1]]
        vf = np.empty((ny+1,nx))
        w = (g.xf[1:-1]-g.xc[:-1])/np.diff(g.xc)
        vf[:,:-1] = v[:,:-1]*(1-w)+v[:,1:]*w
        vf[:,-1] = v[:,-1]
        fx = f.rho*uf*heights[:,None]
        fy = f.rho*vf*widths
        bc = dict(west=(f.speed, g.dx[0]), east=(None,1),
                  south=(0.,g.yc[0]), north=(f.speed,g.yf[-1]-g.yc[-1]))
        force = (p-np.c_[p[:,1:],np.zeros(ny)])*g.dy[:,None]
    else:
        # v unknowns include the open top; wall v=0 is prescribed.
        xs, ys = g.xc, g.yf[1:]
        widths = g.dx
        heights = np.diff(np.r_[g.yc,g.yf[-1]])
        vf = np.r_[(v[:-1]+v[1:])/2, v[-1:]]
        uf = np.empty((ny,nx+1))
        w = ((g.yf[1:-1]-g.yc[:-1])/np.diff(g.yc))[:,None]
        uf[:-1] = u[:-1]*(1-w)+u[1:]*w
        uf[-1] = u[-1]
        fx = f.rho*uf*heights[:,None]
        fy = f.rho*vf*widths
        bc = dict(west=(0.,g.xc[0]), east=(None,1),
                  south=(0.,g.dy[0]), north=(None,1))
        force = (p-np.r_[p[1:],np.zeros((1,nx))])*g.dx
    A, b = transport(xs,ys,widths,heights,fx,fy,f.mu,bc)
    return A, b+force


def predict(A, rhs, old, area, relaxation):
    diagonal = A.diagonal()
    if np.any(diagonal <= 0):
        raise RuntimeError('Nonpositive momentum diagonal.')
    added = diagonal*(1/relaxation-1)
    relaxed = A+diags(added)
    # sum off-diagonal magnitudes excludes eliminated boundary values.
    denominator = np.asarray(A.sum(axis=1)).ravel()+added
    if np.any(denominator <= 0):
        raise RuntimeError('Nonpositive SIMPLEC denominator.')
    predicted = spsolve(relaxed, rhs.ravel()+added*old.ravel()).reshape(old.shape)
    return predicted, np.broadcast_to(area, old.shape)/denominator.reshape(old.shape)


def correct(g, f, u, v, du, dv):
    ny, nx = g.shape
    rows, cols, values = [], [], []
    diagonal = np.zeros(ny*nx)
    for j in range(ny):
        for i in range(nx):
            r = j*nx+i
            for neighbor, coefficient in [
                    (r+1 if i+1<nx else None, f.rho*g.dy[j]*du[j,i]),
                    (r+nx if j+1<ny else None, f.rho*g.dx[i]*dv[j,i])]:
                diagonal[r] += coefficient
                if neighbor is not None:
                    diagonal[neighbor] += coefficient
                    rows.extend([r,neighbor]); cols.extend([neighbor,r])
                    values.extend([-coefficient,-coefficient])
    rows.extend(range(nx*ny)); cols.extend(range(nx*ny)); values.extend(diagonal)
    A = coo_matrix((values,(rows,cols)),shape=(nx*ny,nx*ny)).tocsr()
    rhs = -f.rho*(np.diff(u,axis=1)*g.dy[:,None]+np.diff(v,axis=0)*g.dx)
    pc = spsolve(A,rhs.ravel()).reshape(ny,nx)
    u[:,1:] += du*(pc-np.c_[pc[:,1:],np.zeros(ny)])
    v[1:] += dv*(pc-np.r_[pc[1:],np.zeros((1,nx))])
    return u,v,pc


def solve_flow(nx=40, ny=24, stretching=1.5, height=0.1,
               tolerance=2e-6, max_iterations=1500, physics=None):
    f = physics or Physics()
    g = make_grid(nx,ny,ly=height,stretching=stretching)
    u, v, p = np.full((ny,nx+1),f.speed), np.zeros((ny+1,nx)), np.zeros((ny,nx))
    history = []
    for iteration in range(1,max_iterations+1):
        au,bu = momentum(g,f,u,v,p,'u')
        av,bv = momentum(g,f,u,v,p,'v')
        us,du = predict(au,bu,u[:,1:],g.dy[:,None],0.5)
        vs,dv = predict(av,bv,v[1:],g.dx,0.5)
        u[:,1:],v[1:] = us,vs
        u,v,pc = correct(g,f,u,v,du,dv)
        p += 0.7*pc
        if iteration % 5 == 0 or iteration == 1:
            # Reassemble the unrelaxed nonlinear equations at corrected fields.
            au,bu = momentum(g,f,u,v,p,'u')
            av,bv = momentum(g,f,u,v,p,'v')
            ru = np.max(np.abs(au@u[:,1:].ravel()-bu.ravel())/au.diagonal())/f.speed
            rv = np.max(np.abs(av@v[1:].ravel()-bv.ravel())/av.diagonal())/f.speed
            continuity = np.max(np.abs(g.divergence(u,v)))*height/f.speed
            history.append([iteration,ru,rv,continuity])
            if not np.isfinite([ru,rv,continuity]).all():
                raise RuntimeError('Nonfinite flow residual.')
            if max(ru,rv)<tolerance and continuity<1e-8:
                return g,f,u,v,p,np.array(history)
    raise RuntimeError(f'Flow failed to converge: {history[-1]}')


def solve_energy(g,f,u,v):
    wall = np.where(g.xc<f.heating_start,f.cold,f.hot)
    fx,fy = f.rho*f.cp*u*g.dy[:,None], f.rho*f.cp*v*g.dx
    bc = dict(west=(f.cold,g.xc[0]), east=(None,1),
              south=(wall,g.yc[0]),north=(None,1))
    A,b = transport(g.xc,g.yc,g.dx,g.dy,fx,fy,f.k,bc)
    T = spsolve(A,b.ravel()).reshape(g.shape)
    residual = np.max(np.abs(A@T.ravel()-b.ravel())/A.diagonal())/(f.hot-f.cold)
    qwall = f.k*(wall-T[0])/g.yc[0]
    # Same numerical boundary fluxes as the assembled energy equation.
    inlet = np.sum(fx[:,0]*f.cold+f.k*g.dy/g.xc[0]*(f.cold-T[:,0]))
    outlet = np.sum(fx[:,-1]*T[:,-1])
    top = np.sum(fy[-1]*T[-1])
    heat = np.dot(qwall,g.dx)
    imbalance = abs(inlet+heat-outlet-top)/max(abs(heat),1e-12)
    return T,qwall,residual,imbalance
