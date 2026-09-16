Denote primal cells, dual cells, and corners by $z, p, i$, respectively. Given node coordiantes $X_p$, primal cells pressure $P_z$, primal cells interorior energy $e_z$, corners mass $m_i$, dual cell velocity $v_p$. Then we know primal cells mass $M_z = \sum_{i\in z} m_i$ and dual cell mass $M_p = \sum_{i\in p}m_i$.



### 1. Dual cell momentum update

Calculate every corner's corner vertor
$$
\vec C_i = 1/2 (\vec n_i + \vec n_{i+1}),
$$
where $\vec n_i$ denotes normal vectors of the cell edges that corner coresponds, and can be calculate from node coordiantes $X$.

Primal cell pressure $P_z$ acts on corner $i$ and make corner force
$$
\vec f_i^z = P_z \vec C_i.
$$
Then four corners' resultant force push dual cell $p$ and make accelation
$$
\vec a_p = \frac 1 {M_p} \sum_{i\in p } f_i^z.
$$
So the momuntum is updated
$$
\vec v_p^{n+1} = \vec v_p^n + \tau \vec a_p.
$$

### 2. Primacl cell interior energy update

Corner's velocity is the node's velocity
$$
\vec v_i = \vec v_p, \quad i\in p.
$$
The primal cell expand and do work, losing interior energy
$$
e_z^{n+1} = e_z^n - \frac\tau {M_z} \sum_{i\in z} \vec f_i^z \cdot \vec v_i^{n+1/2},
$$
vice versa.



### 3. Grid geometry update

Update node coordinates
$$
\vec X^{n+1}_p = \vec X_p^n + \tau \vec v_p^{n+1/2},
$$
where mid-time velocity $v^{n+1/2} = 1/2(v^n + v^{n+1})$.



### 4. Primal cell pressure update

Then calculate primal cell pressure by EOS
$$
P^{n+1}_z = P(\rho_z^{n+1}, e_z^{n+1}),
$$
where density $\rho_z^{n+1} = M_z/A_z(X^{n+1})$.





Now we finish time $t=t^{n+1}$,  and have enough information to go to next time step.