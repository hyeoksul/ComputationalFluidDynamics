"""Reproduce operator tests, SIMPLEC refinement, physical comparisons and plots."""
from pathlib import Path
import csv
import json
from time import perf_counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from flat_plate import solve_flow, solve_energy, Physics, transport
from run_verification import main as verify_operators


def write_csv(path, rows):
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer = csv.DictWriter(stream,fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)


def blasius():
    eta = np.linspace(0,12,241)
    initial = np.array([eta-1+np.exp(-eta),1-np.exp(-eta),np.exp(-eta)])
    sol = solve_bvp(lambda x,y: np.array([y[1],y[2],-.5*y[0]*y[2]]),
                    lambda a,b: np.array([a[0],a[1],b[1]-1]),
                    eta,initial,tol=1e-9,max_nodes=5000)
    assert sol.success and abs(sol.y[2,0]-.3320573362)<1e-8
    return sol


def verify_diffusion():
    # Independently prescribed affine harmonic solution on a nonuniform grid.
    from staggered_grid import Grid
    g = Grid(np.linspace(0,1,14)**1.3,np.linspace(0,1,10)**1.6)
    xs,ys = g.xc,g.yc
    bc = dict(west=(2+3*ys,xs[0]),east=(3+3*ys,1-xs[-1]),
              south=(2+xs,ys[0]),north=(5+xs,1-ys[-1]))
    A,b = transport(xs,ys,g.dx,g.dy,np.zeros((len(ys),len(xs)+1)),
                    np.zeros((len(ys)+1,len(xs))),1.,bc)
    exact = 2+xs+3*ys[:,None]
    assert np.max(np.abs(A@exact.ravel()-b.ravel()))<1e-12


def main():
    out = Path(__file__).resolve().parent/'results'
    verify_operators()
    verify_diffusion()
    similarity = blasius()
    f = Physics()
    summaries, profiles, solutions = [],[],{}
    # xi=0.25 is a face in every mesh; no partial heated cells are hidden.
    cases = [('coarse',24,16,1.5,.1),('medium',48,32,1.5,.1),
             ('fine',72,48,1.5,.1),('uniform',48,32,0.,.1),
             ('taller',48,48,1.5,.15)]
    for name,nx,ny,beta,height in cases:
        started = perf_counter()
        print(f'Solving {name}: {nx} x {ny}, H={height}, beta={beta}',flush=True)
        g,f,u,v,p,history = solve_flow(nx,ny,beta,height,physics=f)
        T,q,residual,energy_error = solve_energy(g,f,u,v)
        uc = (u[:,1:]+u[:,:-1])/2
        nu = q*g.xc/(f.k*(f.hot-f.cold))
        cf = 2*f.mu*uc[0]/(g.yc[0]*f.rho*f.speed**2)
        re = f.rho*f.speed*g.xc/f.mu
        reference_cf = 2*similarity.y[2,0]/np.sqrt(re)
        reference_nu = np.full(nx,np.nan)
        heated = g.xc>f.heating_start
        reference_nu[heated] = (.332*f.pr**(1/3)*np.sqrt(re[heated])
            /(1-(f.heating_start/g.xc[heated])**.75)**(1/3))
        comparison = (g.xc>=.35)&(g.xc<=.85)
        nu_error = np.mean(np.abs(nu[comparison]/reference_nu[comparison]-1))
        cf_error = np.mean(np.abs(cf[comparison]/reference_cf[comparison]-1))
        inlet = f.rho*np.dot(u[:,0],g.dy)
        outlet = f.rho*np.dot(u[:,-1],g.dy)
        top = f.rho*np.dot(v[-1],g.dx)
        mass_error = abs(inlet-outlet-top)/inlet
        assert mass_error<1e-9 and energy_error<1e-8 and residual<1e-10
        assert T.min()>=f.cold-1e-7 and T.max()<=f.hot+1e-7
        assert np.min(u[:,-1])>0 and np.min(v[-1])>=-1e-9, 'Backflow requires thermal inflow BC.'
        assert np.all(u[:,0]==f.speed) and np.all(v[0]==0)
        summaries.append(dict(case=name,nx=nx,ny=ny,height_m=height,stretching=beta,
            iterations=int(history[-1,0]),u_residual=history[-1,1],v_residual=history[-1,2],
            max_scaled_divergence=history[-1,3],mass_relative_error=mass_error,
            energy_relative_error=energy_error,temperature_residual=residual,
            net_wall_heat_W_per_m=float(q@g.dx),
            heated_wall_heat_W_per_m=float(q[heated]@g.dx[heated]),
            mean_Nu_relative_error=nu_error,mean_Cf_relative_error=cf_error,
            elapsed_seconds=perf_counter()-started))
        for i,x in enumerate(g.xc):
            profiles.append(dict(case=name,x_m=x,Nu=nu[i],Nu_reference=reference_nu[i],
                                 Cf=cf[i],Cf_reference=reference_cf[i],q_wall_W_m2=q[i]))
        np.savez_compressed(out/f'{name}_fields.npz',x=g.xc,y=g.yc,xf=g.xf,yf=g.yf,
                            u=u,v=v,p=p,T=T,history=history)
        write_csv(out/'flow_verification.csv',summaries)
        print(summaries[-1],flush=True)
        solutions[name] = (g,u,v,T,history,nu,cf,reference_nu,reference_cf)
    write_csv(out/'wall_profiles.csv',profiles)
    with (out/'physics.json').open('w',encoding='utf-8') as stream:
        json.dump(dict(**vars(f),mu=f.mu),stream,indent=2)
    # Accuracy gates are deliberately modest for first-order upwind and a finite domain.
    assert summaries[2]['mean_Nu_relative_error']<.15
    assert summaries[2]['mean_Cf_relative_error']<.15
    assert summaries[2]['mean_Nu_relative_error']<summaries[0]['mean_Nu_relative_error']
    g,u,v,T,history,nu,cf,refnu,refcf = solutions['fine']
    fig,axes = plt.subplots(3,1,figsize=(10,8),layout='constrained')
    for ax,field,title in zip(axes,[(u[:,1:]+u[:,:-1])/2,(v[1:]+v[:-1])/2,T],
                              ['Streamwise velocity (m/s)','Normal velocity (m/s)','Temperature (K)']):
        im = ax.pcolormesh(g.xf,g.yf,field,shading='flat',cmap='viridis')
        fig.colorbar(im,ax=ax); ax.set(xlabel='x (m)',ylabel='y (m)',title=title)
        ax.axvline(f.heating_start,color='white',ls=':',lw=1)
    fig.savefig(out/'flow_and_temperature.png',dpi=180); plt.close(fig)
    fig,axes = plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for name in ['coarse','medium','fine','uniform']:
        gg,*rest = solutions[name]
        nn,cc = solutions[name][5:7]
        axes[0].plot(gg.xc[gg.xc>.25],nn[gg.xc>.25],label=name)
        axes[1].plot(gg.xc,cc,label=name)
    axes[0].plot(g.xc[g.xc>.25],refnu[g.xc>.25],'k--',label='Unheated-start correlation')
    axes[1].plot(g.xc,refcf,'k--',label='Blasius')
    axes[0].set(xlim=(.25,1),xlabel='x (m)',ylabel='Local Nu_x')
    axes[1].set(xlim=(.1,1),ylim=(0,.018),xlabel='x (m)',ylabel='Local C_f')
    for ax in axes: ax.grid(alpha=.3); ax.legend(fontsize=8)
    fig.savefig(out/'wall_benchmarks.png',dpi=180); plt.close(fig)
    fig,ax = plt.subplots(figsize=(7,4),layout='constrained')
    for col,label in [(1,'u momentum'),(2,'v momentum'),(3,'Continuity')]:
        ax.semilogy(history[:,0],np.maximum(history[:,col],1e-17),label=label)
    ax.set(xlabel='SIMPLEC iteration',ylabel='Dimensionless residual'); ax.legend(); ax.grid(alpha=.3)
    fig.savefig(out/'simplec_convergence.png',dpi=180); plt.close(fig)
    fig,ax = plt.subplots(figsize=(6,4),layout='constrained')
    for xpos in [.25,.5,.75]:
        idx = np.argmin(abs(g.xc-xpos))
        eta = g.yc*np.sqrt(f.speed/(f.mu/f.rho*g.xc[idx]))
        ax.plot(((u[:,1:]+u[:,:-1])/2)[:,idx]/f.speed,eta,label=f'x={g.xc[idx]:.3f} m')
    eta = np.linspace(0,8,100)
    ax.plot(similarity.sol(eta)[1],eta,'k--',label='Blasius')
    ax.set(xlabel='u/U_inf',ylabel='Similarity coordinate eta',ylim=(0,8)); ax.legend(); ax.grid(alpha=.3)
    fig.savefig(out/'blasius_profiles.png',dpi=180); plt.close(fig)
    print('All operator, nonlinear convergence, conservation and benchmark checks passed.',flush=True)


if __name__ == '__main__':
    main()
