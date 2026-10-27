For the uploaded 1998 Caramana et al. paper, the key point is: **it is not proposing a new artificial-viscosity law to replace von Neumann–Richtmyer.** It is proposing a way to construct a staggered Lagrangian discretization so that momentum and internal-energy updates are *algebraically compatible*, hence total energy is conserved to roundoff. The artificial-viscosity force may be supplied separately. This is why the paper feels much more abstract than von Neumann's. 

In fact, if you reduce Caramana to 1D and choose the von Neumann artificial pressure \(q\), you get something very close to the algorithm you already coded.

Let cell \(i+\tfrac12\) lie between nodes \(i\) and \(i+1\), and define

$$
\Pi_{i+1/2}=p_{i+1/2}+q_{i+1/2}.
$$

Von Neumann's viewpoint is essentially

$$
q = \text{chosen nonlinear artificial pressure},
$$

then use \(p+q\) in both momentum and thermodynamic work. The original paper specifically constructs \(q\) as a nonlinear viscosity localized to compression, with shock thickness \(O(\Delta x)\). 

Caramana instead says: regard the forces themselves as the fundamental discrete objects.

In 1D, the force exerted by cell \(i+\tfrac12\) on its two endpoints is simply

$$
f^{\,i+1/2}_i=-\Pi_{i+1/2},\qquad
f^{\,i+1/2}_{i+1}=+\Pi_{i+1/2}
$$

(for unit cross-sectional area). In multidimensions these become vector-valued **corner forces** \( {\bf f}_p^z\). All forces acting on a node are summed to obtain its momentum equation. Caramana explicitly allows these corner forces to come from pressure, tensor stress, artificial viscosity, subzonal forces, etc.; their detailed origin is almost arbitrary, subject to momentum conservation. 

So the 1D Caramana momentum step is

$$
M_i
\left(u_i^{n+1}-u_i^n\right)
=
\Delta t
\left(
\Pi_{i-1/2}^{\sigma}
-
\Pi_{i+1/2}^{\sigma}
\right).
\tag{C1}
$$

Nothing surprising yet. This is basically the familiar staggered Lagrangian pressure-force equation.

The important Caramana step is what comes next. Define

$$
u_i^{n+1/2}
=\frac{u_i^{n+1}+u_i^n}{2}.
$$

Then **do not independently discretize**

$$
\frac{de}{dt}=-p\,\frac{dV}{dt}
$$

and hope that it matches the momentum equation. Instead define the energy change from the work of exactly the same discrete forces:

$$
m_{i+1/2}
\left(
e_{i+1/2}^{n+1}-e_{i+1/2}^{n}
\right)
=
-\Delta t\,
\Pi_{i+1/2}^{\sigma}
\left(
u_{i+1}^{n+1/2}-u_i^{n+1/2}
\right).
\tag{C2}
$$

This is the 1D version of their central equation,

$$
\Delta e_z
=
-\frac{\Delta t}{M_z}
\sum_{p\in z}
{\bf f}_p^{z,\sigma}\cdot
{\bf v}_p^{\,n+1/2}.
\tag{20}
$$



That equation is essentially **the whole paper**.

Why is the midpoint velocity essential? Because

$$
{\bf v}^{n+1/2}\cdot
({\bf v}^{n+1}-{\bf v}^{n})
=
\frac12
\left(
|{\bf v}^{n+1}|^2-|{\bf v}^{n}|^2
\right).
$$

Therefore the kinetic-energy increment produced by a corner force is

$$
\Delta K_p
=
\Delta t\,
{\bf f}_p^z\cdot{\bf v}_p^{n+1/2},
$$

while the internal-energy increment assigned to that same cell is

$$
\Delta U_z
=
-\Delta t\,
{\bf f}_p^z\cdot{\bf v}_p^{n+1/2}.
$$

They cancel **algebraically**, not merely to truncation error. That is what Caramana means by “compatible.”

Then move the nodes:

$$
x_i^{n+1}
=
x_i^n+\Delta t\,u_i^{n+1/2},
\tag{C3}
$$

recompute the cell volume from the new coordinates,

$$
V_{i+1/2}^{n+1}
=
x_{i+1}^{n+1}-x_i^{n+1},
\qquad
\rho_{i+1/2}^{n+1}
=
\frac{m_{i+1/2}}{V_{i+1/2}^{n+1}},
$$

and evaluate the EOS.

Their practical implementation is predictor–corrector. In both predictor and corrector the ordering is

$$
\boxed{\text{velocity}\rightarrow
       \text{internal energy}\rightarrow
       \text{coordinates/volume/density}}
$$

because the new velocity is required to construct \(v^{n+1/2}\) before the energy equation can be advanced. They normally do one predictor plus one corrector. 

So, compared directly with von Neumann:

|                      | von Neumann–Richtmyer                     | Caramana compatible SGH                                   |
| -------------------- | ----------------------------------------- | --------------------------------------------------------- |
| Main new idea        | Artificial pressure \(q\) to smear shocks | Compatible discrete force/work pairing                    |
| Primary object       | \(p+q\) in modified PDE                   | Corner force \({\bf f}_p^z\)                              |
| Artificial viscosity | Gives a specific nonlinear \(q\)          | Not specified by this paper; may be essentially arbitrary |
| Momentum             | Pressure/artificial-pressure difference   | Sum of corner forces                                      |
| Internal energy      | Discretize work of \(p+q\)                | **Derived from the same corner forces**                   |
| Energy conservation  | Based on conservative physics of \(p+q\)  | Algebraically enforced in the fully discrete equations    |
| Geometry             | 1D                                        | Arbitrary dimension/coordinate system                     |
| Time update          | Original cycle \(U\to V\to q\to p\)       | \(v\to e\to x,V,\rho\), predictor–corrector               |
| Why corner forces?   | Unnecessary in 1D                         | They identify exactly which zone did work on which node   |

The motivation can therefore be reduced to one sentence:

> **In multiple dimensions, “discretize \(-\nabla p\)” and “discretize \(-p\nabla\!\cdot v\)” independently is dangerous, because the discrete gradient/divergence no longer automatically satisfy the continuum integration-by-parts identity; Caramana makes them conjugate by using one common corner force for momentum and work.**

That is why the method is almost trivial in 1D. In 1D,

$$
-\nabla p \leftrightarrow
p_{i-1/2}-p_{i+1/2}
$$

and

$$
-p\nabla\cdot u
\leftrightarrow
-p_{i+1/2}(u_{i+1}-u_i)
$$

are naturally paired. Once you go to distorted quadrilaterals/hexahedra, tensor AV, subzonal pressures, etc., that pairing is no longer automatic. The corner-force abstraction restores it. The authors explicitly describe the common forces as exchanging nodal kinetic energy and zonal internal energy. 

So if you already have a working 1D VNR code, the shortest route to understanding/implementing Caramana is:

1. Keep your staggered nodes/cells and even keep your VNR \(q\).
2. Think of each cell pressure \(p+q\) as generating two endpoint forces.
3. Sum endpoint forces to update nodal velocity.
4. Compute \(u^{n+1/2}=(u^n+u^{n+1})/2\).
5. Update each cell's internal energy by **minus the work of its two endpoint forces**, rather than by an independently chosen discretization of \(p\,dV\).
6. In 2D/3D, replace those two endpoint forces by several vector corner forces; the exact same logic remains.

That is the algorithmic content you need before going into Caramana's discussions of support operators, entropy/volume consistency, area-weighted differencing, boundary conditions, etc. Those sections explain consequences and geometric subtleties; they are not another core algorithm.
