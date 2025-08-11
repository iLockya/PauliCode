# PauliCode
IBM's [[144,12,12]] code
$$
f = x+x^2+y^3, \quad g = y+y^2+x^3
$$
In the thermodynamic limit.

```julia
using PauliCode
C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"]); # define the code

ϵX = PauliCode._x_excitation_map(C)
image(C.X)[1] == kernel(ϵX)[1] # topological condition
# true

Q = topological_index(C)
# 8
MobilitySubLattice(C)
#  l₀ => 12
#  m₀ => 12
```

In a finite-size lattice, if the size is not a multiple of $\ell_0$and $m_0$, the logical dimension does not saturate $2Q$. This is the phenomenon of *topological frustration*.

```julia
C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"], [12,12]); # 12 x 12 lattice
logical_qubits(C) == 16 # true

C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"], [12,6]); # 12 x 6 lattice
logical_qubits(C) == 12 # true

C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"], [12,3]); # 12 x 3 lattice
logical_qubits(C) == 8 # true

C, (x,y) = CSSCode(["x^3+y^2+y","y^3+x^2+x"], [12,7]); # 12 x 7 lattice
logical_qubits(C) == 0 # true
```

$(\ell_0, m_0)$ may not the unit cells of *mobility sublattices*. For example $\mathbb{BB}[\overline{1},\overline{1},3,-3]$ code

```julia
R, (x,y,x̄,ȳ) = polynomial_ring(GF(2), [:x,:y,:x̄,:ȳ])
I = ideal([x*x̄-1,y*ȳ-1, 1+x+x̄*ȳ^3,1+y+ȳ*x^3])
L = exponent_lattice(I,[x,y])
L.basis_matrix
# [6   360]
# [0   762]
```

The mobility sublattice $\Lambda_{\mathcal{C}} = \langle x^6y^{360}, y^{762}\rangle$.
