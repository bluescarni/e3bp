# e3bp
Python code for the solution of the three-dimensional problem of two fixed centres

The code has the following prerequisites:

* the w_elliptic library, for the computation of Weierstrass elliptic and related functions,
  https://github.com/bluescarni/w_elliptic,
* the mpmath library, http://mpmath.org/,
* the matplotlib/numpy stack.

Usage example:

```python
import e3bp
e = e3bp.e3bp(1.1,1.2,1.3,[.1,.3,.4,.1,-.1,.2])
```

This will initialise an object ``e`` corresponding to two fixed centres problem with the following parameters:

* ``a = 1.1``,
* ``mu_1 = 1.2``,
* ``mu_2 = 1.3``,
* initial cartesian position vector ``[.1,.3,.4]``,
* initial cartesian velocity vector ``[.1,-.1,.2]``.

The ``e`` object can then be used as follows:

```python
e.xi_tau(0.1) # Compute the xi coordinate at tau = 0.1
e.eta_tau(0.2) # Compute the eta coordinate at tau = 0.2
e.ell_ham_state_tau(0.3) # Compute the Hamiltonian state at tau = 0.3
```
