# Wrappers for various math functions.
def sqrt(x):
    if isinstance(x,float):
        from math import sqrt
        return sqrt(x)
    if isinstance(x,complex):
        from cmath import sqrt
        return sqrt(x)
    from mpmath import sqrt
    return sqrt(x)

def atan2(x,y):
    if isinstance(x,float) and isinstance(y,float):
        from math import atan2
        return atan2(x,y)
    from mpmath import atan2
    return atan2(x,y)

def cos(x):
    if isinstance(x,float):
        from math import cos
        return cos(x)
    if isinstance(x,complex):
        from cmath import cos
        return cos(x)
    from mpmath import cos
    return cos(x)

def sin(x):
    if isinstance(x,float):
        from math import sin
        return sin(x)
    if isinstance(x,complex):
        from cmath import sin
        return sin(x)
    from mpmath import sin
    return sin(x)

def polyval(l,x):
    from mpmath import polyval, mpf, mpc
    if any([isinstance(_,(mpf,mpc)) for _ in l + [x]]):
        return polyval(l,x)
    retval = 0.
    tmp = 1.
    for c in l[-1::-1]:
        retval += c*tmp;
        tmp *= x
    return retval

def polyroots(l):
    from mpmath import polyroots, mpc, mpf
    retval = polyroots(l,maxsteps=200,extraprec=40)
    if not any([isinstance(_,(mpf,mpc)) for _ in l]):
        retval = [float(_) if isinstance(_,mpf) else complex(_) for _ in retval]
    return retval

# Functions for coordinate transformations.
def cart2cyl(state):
    x,y,z,vx,vy,vz = state
    rho = sqrt(x**2 + y**2)
    phi = atan2(y,x)
    vrho = (x*vx+y*vy)/rho
    vphi = (vy*x-vx*y)/(rho**2)
    return rho,phi,z,vrho,vphi,vz

def cyl2cart(state):
    rho,phi,z,vrho,vphi,vz = state
    x = rho*cos(phi)
    y = rho*sin(phi)
    vx = vrho*cos(phi)-rho*vphi*sin(phi)
    vy = vrho*sin(phi)+rho*vphi*cos(phi)
    return x,y,z,vx,vy,vz

def cyl2ellcyl(state,a):
    rho,phi,z,vrho,vphi,vz = state
    xi  = (sqrt(rho**2+(z+a)**2)+sqrt(rho**2+(z-a)**2))/(2*a)
    eta = (sqrt(rho**2+(z+a)**2)-sqrt(rho**2+(z-a)**2))/(2*a)
    vxi  = 1/(2*a)*((vrho*rho+vz*(z+a))/sqrt(rho**2+(z+a)**2)+(vrho*rho+vz*(z-a))/sqrt(rho**2+(z-a)**2))
    veta = 1/(2*a)*((vrho*rho+vz*(z+a))/sqrt(rho**2+(z+a)**2)-(vrho*rho+vz*(z-a))/sqrt(rho**2+(z-a)**2))
    return xi,eta,phi,vxi,veta,vphi

def ellcyl2cyl(state,a):
    xi,eta,phi,vxi,veta,vphi = state
    rho = a*sqrt((xi**2-1)*(1-eta**2))
    z = a*xi*eta
    vrho = a*(vxi*xi*(1-eta**2)-veta*eta*(xi**2-1))/sqrt((xi**2-1)*(1-eta**2))
    vz = a*(vxi*eta+xi*veta)
    return rho,phi,z,vrho,vphi,vz

def ellcyl2ham(state,a):
    xi,eta,phi,vxi,veta,vphi = state
    pxi = a**2*vxi*(xi**2-eta**2)/(xi**2-1)
    peta = a**2*veta*(xi**2-eta**2)/(1-eta**2)
    pphi = a**2*vphi*(xi**2-1)*(1-eta**2)
    return pxi,peta,pphi,xi,eta,phi

def ham2ellcyl(state,a):
    pxi,peta,pphi,xi,eta,phi = state
    vxi = pxi*(xi**2-1)/(a**2*(xi**2-eta**2))
    veta = peta*(1-eta**2)/(a**2*(xi**2-eta**2))
    vphi = pphi/(a**2*(xi**2-1)*(1-eta**2))
    return xi,eta,phi,vxi,veta,vphi

# The integrals of motion.
def integrals(a,mu1,mu2,cart_state):
    pxi,peta,pphi,xi,eta,phi = ellcyl2ham(cyl2ellcyl(cart2cyl(cart_state),a),a)
    h = 1/(2*a**2*(xi**2-eta**2))*(pxi**2*(xi**2-1)+peta**2*(1-eta**2)) + pphi**2/(2*a**2*(xi**2-1)*(1-eta**2)) \
        -mu1/(a*(xi-eta))-mu2/(a*(xi+eta))
    hxi = -xi/a*(mu1+mu2)+pxi**2/(2*a**2)*(xi**2-1)-xi**2*h+pphi**2/(2*a**2*(xi**2-1))
    heta = -eta/a*(mu1-mu2)+peta**2/(2*a**2)*(1-eta**2)+eta**2*h+pphi**2/(2*a**2*(1-eta**2))
    return h,hxi,heta,pphi

# Get the coefficients of the f_xi polynomial.
def get_fxi_cf(a,mu1,mu2,cart_state):
    h,hxi,_,pphi = integrals(a,mu1,mu2,cart_state)
    c4 = 2*a**2*h
    c3 = 2*mu1*a+2*mu2*a
    c2 = -2*a**2*h+2*a**2*hxi
    c1 = -2*mu1*a-2*mu2*a
    c0 = -2*a**2*hxi-pphi**2
    return c4,c3,c2,c1,c0

# Get a function for the evaluation of the polynomial fxi, its derivatives, its roots
# and the roots of its derivative.
def get_fxi(a,mu1,mu2,cart_state):
    c4,c3,c2,c1,c0 = get_fxi_cf(a,mu1,mu2,cart_state)
    def retval(xi):
        return polyval([c4,c3,c2,c1,c0],xi)
    def retval_p(xi):
        return polyval([4*c4,3*c3,2*c2,c1],xi)
    def retval_pp(xi):
        return polyval([4*3*c4,3*2*c3,2*c2],xi)
    return retval,retval_p,retval_pp,polyroots([c4,c3,c2,c1,c0]),polyroots([4*c4,3*c3,2*c2,c1])

# Same above, for eta.
def get_feta_cf(a,mu1,mu2,cart_state):
    h,_,heta,pphi = integrals(a,mu1,mu2,cart_state)
    c4 = 2*a**2*h
    c3 = -2*mu1*a+2*mu2*a
    c2 = -2*a**2*h-2*a**2*heta
    c1 = 2*mu1*a-2*mu2*a
    c0 = 2*a**2*heta-pphi**2
    return c4,c3,c2,c1,c0

def get_feta(a,mu1,mu2,cart_state):
    c4,c3,c2,c1,c0 = get_feta_cf(a,mu1,mu2,cart_state)
    def retval(eta):
        return polyval([c4,c3,c2,c1,c0],eta)
    def retval_p(eta):
        return polyval([4*c4,3*c3,2*c2,c1],eta)
    def retval_pp(eta):
        return polyval([4*3*c4,3*2*c3,2*c2],eta)
    return retval,retval_p,retval_pp,polyroots([c4,c3,c2,c1,c0]),polyroots([4*c4,3*c3,2*c2,c1])

# Main class.
class e3bp(object):
    @staticmethod
    def __invariants_from_cfs(c4,c3,c2,c1,c0):
        a0,a1,a2,a3,a4 = c4,c3/4,c2/6,c1/4,c0
        g2 = a0*a4-4*a1*a3+3*a2**2
        g3 = a0*a2*a4+2*a1*a2*a3-a2**3-a0*a3**2-a1**2*a4
        return g2,g3
    def __init__(self,a,mu1,mu2,cart_state):
        from mpmath import mpf, mpc
        if not any([isinstance(_,(mpf,mpc)) for _ in [a,mu1,mu2]+cart_state]):
            from w_elliptic import we
            def make_fp(x):
                if not isinstance(x,(complex,float)):
                    return float(x)
                return x
            def pi():
                from math import pi
                return pi
        else:
            from weierstrass_elliptic import weierstrass_elliptic as we
            def make_fp(x):
                if isinstance(x,complex):
                    return mpc(x)
                if isinstance(x,mpc):
                    return x
                return mpf(x)
            def pi():
                from mpmath import pi
                return pi()
        a,mu1,mu2 = [make_fp(_) for _ in [a,mu1,mu2]]
        if isinstance(a,(mpc,complex)):
            is_complex = True
        else:
            is_complex = False
        cart_state = tuple([make_fp(_) for _ in cart_state])
        # Store the initial parameters.
        self.__params = [a,mu1,mu2,cart_state]
        # Store the initial conditions in elliptic coords.
        self.__ell_state = cyl2ellcyl(cart2cyl(self.cart_state),a)
        # Store the initial conditions in elliptic Hamiltonian coords.
        self.__ell_ham_state  = ellcyl2ham(self.ell_state,a)
        if self.__ell_ham_state[2] == 0.:
            raise ValueError("currently, the 2d case cannot be handled")
        # Store the integrals.
        self.__integrals = integrals(a,mu1,mu2,cart_state)
        # Calculate the invariants and poly quantities for xi.
        c4,c3,c2,c1,c0 = get_fxi_cf(a,mu1,mu2,cart_state)
        self.__fxi = get_fxi(a,mu1,mu2,cart_state)
        # TODO: make this a public member (and feta too?)
        self.fxi = self.__fxi
        g2_xi,g3_xi = e3bp.__invariants_from_cfs(c4,c3,c2,c1,c0)
        # And eta.
        c4,c3,c2,c1,c0 = get_feta_cf(a,mu1,mu2,cart_state)
        self.__feta = get_feta(a,mu1,mu2,cart_state)
        self.feta = self.__feta
        g2_eta,g3_eta = e3bp.__invariants_from_cfs(c4,c3,c2,c1,c0)
        # Build the wp objects.
        # NOTE: the invariants are always real by construction.
        self.__wp_xi = we(g2_xi.real,g3_xi.real)
        self.__wp_eta = we(g2_eta.real,g3_eta.real)
        # Compute the times of passage.
        tau_xi_list = [(r,a**2 * self.__wp_xi.Pinv(self.__fxi[2](r)/24 + self.__fxi[1](r)/(4*(self.__ell_state[0]-r)))[0]) for r in self.__fxi[3]]
        tau_eta_list = [(r,a**2 * self.__wp_eta.Pinv(self.__feta[2](r)/24 + self.__feta[1](r)/(4*(self.__ell_state[1]-r)))[0]) for r in self.__feta[3]]
        # Filter out the following:
        # - root values which are complex or not reachable,
        # - tau values which are complex.
        tau_xi_list = list(filter(lambda t: t[0].imag == 0 and t[0].real > 1 and t[1].imag == 0,tau_xi_list))
        tau_eta_list = list(filter(lambda t: t[0].imag == 0 and t[0].real < 1 and t[0].real > -1 and t[1].imag == 0,tau_eta_list))
        # Pick the smallest reachable root, the smallest tau_xi and adjust its sign
        # according to the sign of the initial p_xi (or dxi/dtau).
        tau_xi_list = list(min(tau_xi_list,key = lambda t: t[0]))
        if self.ell_ham_state[0] >= 0:
            tau_xi_list[1] = -tau_xi_list[1]
        tau_xi_list[1] = tau_xi_list[1].real
        # Same for eta.
        tau_eta_list = list(min(tau_eta_list,key = lambda t: t[0]))
        if self.ell_ham_state[1] >= 0:
            tau_eta_list[1] = -tau_eta_list[1]
        tau_eta_list[1] = tau_eta_list[1].real
        # Implementation of xi_tau and pxi_tau.
        xi_tau = lambda tau, xi_r, tau_xi : xi_r + self.__fxi[1](xi_r)/(4*(self.wp_xi.P((tau-tau_xi)/a**2)-self.__fxi[2](xi_r)/24))
        pxi_tau = lambda tau, xi_r, tau_xi: -a**2/(xi_tau(tau,xi_r,tau_xi)**2 - 1) * self.__fxi[1](xi_r) * self.wp_xi.Pprime((tau-tau_xi)/a**2) / (4*a**2*(self.wp_xi.P((tau-tau_xi)/a**2)-self.__fxi[2](xi_r)/24)**2)
        self.__xi_tau = lambda tau: xi_tau(tau,*tau_xi_list)
        self.__pxi_tau = lambda tau: pxi_tau(tau,*tau_xi_list)
        # Now eta/peta.
        eta_tau = lambda tau, eta_r, tau_eta : eta_r + self.__feta[1](eta_r)/(4*(self.wp_eta.P((tau-tau_eta)/a**2)-self.__feta[2](eta_r)/24))
        peta_tau = lambda tau, eta_r, tau_eta: -a**2/(1-eta_tau(tau,eta_r,tau_eta)**2) * self.__feta[1](eta_r) * self.wp_eta.Pprime((tau-tau_eta)/a**2) / (4*a**2*(self.wp_eta.P((tau-tau_eta)/a**2)-self.__feta[2](eta_r)/24)**2)
        self.__eta_tau = lambda tau: eta_tau(tau,*tau_eta_list)
        self.__peta_tau = lambda tau: peta_tau(tau,*tau_eta_list)
        # Now phi.
        # It's useful to define the general form for J1.
        def J1(u,v,wp,Ppv,zv):
            # Ppv is by definition purely real or imaginary.
            if abs(Ppv.real) > abs(Ppv.imag):
                Ppv = Ppv.real
                return 1/Ppv * (wp.ln_sigma_real(v-u) - wp.ln_sigma_real(u+v) + 2*u*zv.real)
            else:
                Ppv = Ppv.imag
                # NOTE: the -1 factor here is cancelled by the fact that we take the inverse of Ppv (so -1/(i*b) == i/b).
                return 1/Ppv * (wp.ln_sigma_imag(v-u) - wp.ln_sigma_imag(u+v) + pi() + 2*u*zv.imag)
        # Some useful shortcuts.
        pphi = self.ell_ham_state[2]
        xi_r, tau_xi = tau_xi_list
        self.xi_r = xi_r
        self.tau_xi = tau_xi
        Axi = self.__fxi[1](xi_r)
        Bxi = self.__fxi[2](xi_r)/24
        vxi1 = self.__wp_xi.Pinv(-(Axi-4*Bxi*xi_r-4*Bxi)/(4*(xi_r+1)))[0]
        Ppvxi1 = self.__wp_xi.Pprime(vxi1)
        zvxi1 = self.__wp_xi.zeta(vxi1)
        vxi2 = self.__wp_xi.Pinv(-(Axi-4*Bxi*xi_r+4*Bxi)/(4*(xi_r-1)))[0]
        Ppvxi2 = self.__wp_xi.Pprime(vxi2)
        zvxi2 = self.__wp_xi.zeta(vxi2)
        eta_r, tau_eta = tau_eta_list
        self.eta_r = eta_r
        self.tau_eta = tau_eta
        Aeta = self.__feta[1](eta_r)
        Beta = self.__feta[2](eta_r)/24
        veta1 = self.__wp_eta.Pinv(-(Aeta-4*Beta*eta_r-4*Beta)/(4*(eta_r+1)))[0]
        Ppveta1 = self.__wp_eta.Pprime(veta1)
        zveta1 = self.__wp_eta.zeta(veta1)
        veta2 = self.__wp_eta.Pinv(-(Aeta-4*Beta*eta_r+4*Beta)/(4*(eta_r-1)))[0]
        Ppveta2 = self.__wp_eta.Pprime(veta2)
        zveta2 = self.__wp_eta.zeta(veta2)
        # TODO these are for the average linear evolution of phi, see if we need to keep them.
        self.phi_xi = lambda tau: self.ell_ham_state[5] + pphi*Axi/(8*(xi_r+1)**2) * (J1((tau-tau_xi)/a**2,vxi1,self.wp_xi,Ppvxi1,zvxi1)-J1(-tau_xi/a**2,vxi1,self.wp_xi,Ppvxi1,zvxi1)) \
            - pphi*Axi/(8*(xi_r-1)**2) * (J1((tau-tau_xi)/a**2,vxi2,self.wp_xi,Ppvxi2,zvxi2)-J1(-tau_xi/a**2,vxi2,self.wp_xi,Ppvxi2,zvxi2)) + pphi*tau/(a**2*(xi_r**2-1))
        self.phi_eta = lambda tau: - pphi*Aeta/(8*(eta_r+1)**2) * (J1((tau-tau_eta)/a**2,veta1,self.wp_eta,Ppveta1,zveta1)-J1(-tau_eta/a**2,veta1,self.wp_eta,Ppveta1,zveta1)) \
            + pphi*Aeta/(8*(eta_r-1)**2) * (J1((tau-tau_eta)/a**2,veta2,self.wp_eta,Ppveta2,zveta2)-J1(-tau_eta/a**2,veta2,self.wp_eta,Ppveta2,zveta2)) \
            - pphi*tau/(a**2*(eta_r**2-1))
        self.__phi_tau = lambda tau: self.ell_ham_state[5] + pphi*Axi/(8*(xi_r+1)**2) * (J1((tau-tau_xi)/a**2,vxi1,self.wp_xi,Ppvxi1,zvxi1)-J1(-tau_xi/a**2,vxi1,self.wp_xi,Ppvxi1,zvxi1)) \
            - pphi*Axi/(8*(xi_r-1)**2) * (J1((tau-tau_xi)/a**2,vxi2,self.wp_xi,Ppvxi2,zvxi2)-J1(-tau_xi/a**2,vxi2,self.wp_xi,Ppvxi2,zvxi2)) \
            - pphi*Aeta/(8*(eta_r+1)**2) * (J1((tau-tau_eta)/a**2,veta1,self.wp_eta,Ppveta1,zveta1)-J1(-tau_eta/a**2,veta1,self.wp_eta,Ppveta1,zveta1)) \
            + pphi*Aeta/(8*(eta_r-1)**2) * (J1((tau-tau_eta)/a**2,veta2,self.wp_eta,Ppveta2,zveta2)-J1(-tau_eta/a**2,veta2,self.wp_eta,Ppveta2,zveta2)) \
            + pphi*tau/(a**2*(xi_r**2-1)) - pphi*tau/(a**2*(eta_r**2-1))
        # Real time.
        def J2(u, v, wp, Pv, Ppv, zv):
            g2 = wp.invariants[0]
            retval = -2*wp.zeta(u)-wp.Pprime(u)/(wp.P(u)-Pv)-2*u*Pv-(6 * Pv**2 - g2/2) * J1(u,v,wp,Ppv,zv)
            if abs(Ppv.real) > abs(Ppv.imag):
                return 1/Ppv.real**2 * retval
            else:
                return -1/Ppv.imag**2 * retval
        bxi = self.__wp_xi.Pinv(Bxi)[0]
        Pbxi = Bxi
        Ppbxi = self.__wp_xi.Pprime(bxi)
        zetabxi = self.__wp_xi.zeta(bxi)
        beta = self.__wp_eta.Pinv(Beta)[0]
        Pbeta = Beta
        Ppbeta = self.__wp_eta.Pprime(beta)
        zetabeta = self.__wp_eta.zeta(beta)
        self.__t_tau = lambda tau: a**2*Axi**2/16 * (J2((tau-tau_xi)/a**2,bxi,self.wp_xi,Pbxi,Ppbxi,zetabxi)-J2(-tau_xi/a**2,bxi,self.wp_xi,Pbxi,Ppbxi,zetabxi)) \
            + a**2*Axi*xi_r/2 * (J1((tau-tau_xi)/a**2,bxi,self.wp_xi,Ppbxi,zetabxi)-J1(-tau_xi/a**2,bxi,self.wp_xi,Ppbxi,zetabxi)) \
            - a**2*Aeta**2/16 * (J2((tau-tau_eta)/a**2,beta,self.wp_eta,Pbeta,Ppbeta,zetabeta)-J2(-tau_eta/a**2,beta,self.wp_eta,Pbeta,Ppbeta,zetabeta)) \
            - a**2*Aeta*eta_r/2 * (J1((tau-tau_eta)/a**2,beta,self.wp_eta,Ppbeta,zetabeta)-J1(-tau_eta/a**2,beta,self.wp_eta,Ppbeta,zetabeta)) \
            + (xi_r**2 - eta_r**2) * tau
        # Check if motion is bounded.
        self.__bounded = self.xi_tau(self.a**2*self.wp_xi.periods[0]/2+tau_xi_list[1]).real > (self.__fxi[2](tau_xi_list[0])/24).real
        if not self.bounded:
            Pi = self.wp_xi.Pinv(self.__fxi[2](tau_xi_list[0])/24)
            # We compute the tau_infs with both the positive and negative values of Pinv, as the reachable asymptote
            # will depend on the sign of the initial velocities. In any case, the reachable asymptotes will be the 2 closest
            # ones in absolute value to tau = 0.
            tau_infs = [_.real for _ in [tau_xi + self.a**2*Pi[0],tau_xi + self.a**2*Pi[1],tau_xi - self.a**2*Pi[0],tau_xi - self.a**2*Pi[1]]]
            # Sort first according to abs(tau), then order by tau.
            tau_infs = sorted(tau_infs,key=lambda _:abs(_))[0:2]
            self.__tau_bounds = sorted(tau_infs)
    def av_phi_plot(self,tfin):
        import matplotlib.pyplot as plt
        from matplotlib import rc
        from pylab import xlabel, ylabel
        import numpy as np
        rc('text', usetex=True)
        r = np.linspace(0,tfin,3000)
        phi = [self.phi_tau(_) for _ in r]
        T_xi = self.wp_xi.periods[0].real
        T_eta = self.wp_eta.periods[0].real
        Phixi = (self.phi_xi(T_xi)-self.phi_xi(0.))/T_xi
        Phieta = (self.phi_eta(T_eta)-self.phi_eta(0.))/T_eta
        phi_0 = self.ell_ham_state[5]
        phi_av = [phi_0 + (Phixi + Phieta) * _ for _ in r]
        plt.plot(r,phi,color='black')
        plt.plot(r,phi_av,'--',color='black')
        xlabel(r'$\tau$')
        ylabel(r'$\phi$')
    @property
    def a(self):
        return self.__params[0]
    @property
    def mu1(self):
        return self.__params[1]
    @property
    def mu2(self):
        return self.__params[2]
    @property
    def cart_state(self):
        return self.__params[3]
    @property
    def ell_state(self):
        return self.__ell_state
    @property
    def ell_ham_state(self):
        return self.__ell_ham_state
    @property
    def integrals(self):
        return self.__integrals
    @property
    def wp_xi(self):
        return self.__wp_xi
    @property
    def wp_eta(self):
        return self.__wp_eta
    @property
    def xi_roots(self):
        return self.__fxi[3]
    @property
    def eta_roots(self):
        return self.__feta[3]
    @property
    def xi_tau(self):
        return self.__xi_tau
    @property
    def eta_tau(self):
        return self.__eta_tau
    @property
    def pxi_tau(self):
        return self.__pxi_tau
    @property
    def peta_tau(self):
        return self.__peta_tau
    @property
    def phi_tau(self):
        return self.__phi_tau
    @property
    def t_tau(self):
        return self.__t_tau
    def ell_ham_state_tau(self,tau):
        from numpy import array
        return array([self.pxi_tau(tau),self.peta_tau(tau),self.ell_ham_state[2],self.xi_tau(tau),self.eta_tau(tau),self.phi_tau(tau)])
    @property
    def bounded(self):
        return self.__bounded
    @property
    def tau_infs(self):
        if self.bounded:
            raise ValueError("the trajectory is bounded")
        return self.__tau_bounds
    def __repr__(self):
        retval = ""
        retval += "a,mu1,mu2 = ({0},{1},{2})\n".format(self.a,self.mu1,self.mu2)
        retval += "Initial cartesian state = {0}\n".format(self.cart_state)
        retval += "Initial elliptic state = {0}\n".format(self.ell_state)
        retval += "Initial Hamiltonian state = {0}\n".format(self.ell_ham_state)
        retval += "Bounded = {0}\n".format(self.bounded)
        if not self.bounded:
            retval += "Tau infs = {0}\n".format(self.tau_infs)
        return retval
    # Numerical integration in elliptic coordinates and
    # fictitious time.
    def numerical_ell_ham(self,tau1,dtau):
        from mpmath import mpf,mpc
        from scipy.integrate import ode
        import numpy as np
        a = self.a
        mu1 = self.mu1
        mu2 = self.mu2
        is_mp = isinstance(a,(mpf,mpc))
        h = self.integrals[0]
        def f_int(_,y):
            # NOTE: the Ham variables are inverted, first momenta
            # then coordinates.
            if is_mp:
                pxi,peta,pphi,xi,eta,phi = [mpc(_) for _ in y]
            else:
                pxi,peta,pphi,xi,eta,phi = y
            pxi_tau = -(-2*xi*h-(mu1+mu2)/a-pphi**2*xi/(a**2*(xi**2-1)**2)+pxi**2*xi/a**2)
            peta_tau = -(2*eta*h-(mu1-mu2)/a+pphi**2*eta/(a**2*(1-eta**2)**2)-peta**2*eta/a**2)
            pphi_tau = 0
            xi_tau = pxi*(xi**2-1)/a**2
            eta_tau = peta*(1-eta**2)/a**2
            phi_tau = pphi/a**2*(1/(xi**2-1)+1/(1-eta**2))
            if is_mp:
                return np.array([float(_) for _ in [pxi_tau.real,peta_tau.real,pphi_tau.real,xi_tau.real,eta_tau.real,phi_tau.real]])
            else:
                return np.array([pxi_tau,peta_tau,pphi_tau,xi_tau,eta_tau,phi_tau])
        r = ode(f_int).set_integrator('dopri5',rtol=1E-14,nsteps=5000)
        if is_mp:
            r.set_initial_value(np.array([float(_.real) for _ in self.ell_ham_state]),0.)
        else:
            r.set_initial_value(np.array(self.ell_ham_state),0.)
        t_ret = []
        y_ret = []
        while r.successful() and r.t < tau1:
            r.integrate(r.t+dtau)
            t_ret.append(r.t)
            y_ret.append(r.y)
        return np.array(t_ret),np.array(y_ret)
    def c_integrals(self,ell_ham_state):
        a,mu1,mu2 = self.a, self.mu1, self.mu2
        pxi,peta,pphi,xi,eta,phi = ell_ham_state
        h = 1/(2*a**2*(xi**2-eta**2))*(pxi**2*(xi**2-1)+peta**2*(1-eta**2)) + pphi**2/(2*a**2*(xi**2-1)*(1-eta**2)) \
            -mu1/(a*(xi-eta))-mu2/(a*(xi+eta))
        hxi = -xi/a*(mu1+mu2)+pxi**2/(2*a**2)*(xi**2-1)-xi**2*h+pphi**2/(2*a**2*(xi**2-1))
        heta = -eta/a*(mu1-mu2)+peta**2/(2*a**2)*(1-eta**2)+eta**2*h+pphi**2/(2*a**2*(1-eta**2))
        return h,hxi,heta,pphi
    def plot_2d(self,tfin,npoints=3000):
        import matplotlib.pyplot as plt
        from matplotlib import rc
        from pylab import xlabel, ylabel
        import numpy as np
        if not self.bounded:
            raise ValueError("the trajectory is not bounded")
        rc('text', usetex=True)
        xi_min,xi_max = sorted([self.xi_tau(self.a**2*self.wp_xi.periods[0].real/2.+self.tau_xi),self.xi_tau(self.tau_xi)])
        eta_min,eta_max = sorted([self.eta_tau(self.a**2*self.wp_eta.periods[0].real/2.+self.tau_eta),self.eta_tau(self.tau_eta)])
        # We know that the interesting part is within the elliptical donut. Thus we establish where the outer boundary
        # of the donut interstects the two axes.
        rho_max = self.a*sqrt(xi_max**2-1)
        rho_max += rho_max/20.
        z_max = self.a*xi_max
        z_max += z_max/20.
        z,rho=np.ogrid[-z_max:z_max:npoints*1j,0:rho_max:npoints*1j]
        opts = {'colors':'black','antialiased':True,'linewidths':0.5,'linestyles':'dashed'}
        xlabel(r'$\rho$')
        ylabel(r'$z$')
        plt.contour(rho.ravel(),z.ravel(),(np.sqrt(rho**2+(z+self.a)**2)+np.sqrt(rho**2+(z-self.a)**2))/(2*self.a)-xi_max,[0],**opts)
        plt.contour(rho.ravel(),z.ravel(),(np.sqrt(rho**2+(z+self.a)**2)+np.sqrt(rho**2+(z-self.a)**2))/(2*self.a)-xi_min,[0],**opts)
        plt.contour(rho.ravel(),z.ravel(),(np.sqrt(rho**2+(z+self.a)**2)-np.sqrt(rho**2+(z-self.a)**2))/(2*self.a)-eta_max,[0],**opts)
        plt.contour(rho.ravel(),z.ravel(),(np.sqrt(rho**2+(z+self.a)**2)-np.sqrt(rho**2+(z-self.a)**2))/(2*self.a)-eta_min,[0],**opts)
        r = np.linspace(0,tfin,npoints)
        plt.plot([self.a*sqrt((self.xi_tau(_)**2-1)*(1-self.eta_tau(_)**2)) for _ in r],[self.a*self.xi_tau(_)*self.eta_tau(_) for _ in r],color='black')
    def plot(self,tfin,npoints=3000):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        import numpy as np
        fig = plt.figure()
        fig.add_subplot(121)
        self.plot_2d(tfin,npoints)
        ax = fig.add_subplot(122, projection='3d')
        # This can be useful to tweak, keep it in mind.
        # ax.set_zlim(-3,3)
        r = np.linspace(0,tfin,npoints)
        ell_state = np.array([self.ell_ham_state_tau(_) for _ in r])
        cart_state = np.array([cyl2cart(ellcyl2cyl(ham2ellcyl(s,self.a),self.a)) for s in ell_state])
        ax.plot(cart_state[:,0],cart_state[:,1],cart_state[:,2],color='black')
    def plot_panel(self,npoints=3000):
        from numpy import linspace
        import matplotlib.pyplot as plt
        from matplotlib import rc
        from pylab import xlabel, ylabel, xlim, ylim
        rc('text', usetex=True)
        fig = plt.figure()
        opts = 'k'
        if self.bounded:
            tfin = self.a**2 * (self.wp_xi.periods[0].real + self.wp_eta.periods[0].real)
            l = linspace(0,tfin,npoints)
        else:
            #T = self.tau_infs[1] - self.tau_infs[0]
            #l = linspace(self.tau_infs[0] + T/10.,self.tau_infs[1] - T/10.)
            #tfin = self.a**2 * (self.wp_xi.periods[0].real + self.wp_eta.periods[0].real)
            delta_tau = self.tau_infs[1] - self.tau_infs[0]
            l = linspace(self.tau_infs[0]+.01*delta_tau,self.tau_infs[1]-0.01*delta_tau,npoints)
        # xi.
        fig.add_subplot(221)
        plt.plot(l,[self.xi_tau(_) for _ in l],opts)
        xlabel(r'$\tau$')
        ylabel(r'$\xi$')
        xlim(l[0],l[-1])
        #ylim(1,5)
        # eta.
        fig.add_subplot(222)
        plt.plot(l,[self.eta_tau(_) for _ in l],opts)
        xlabel(r'$\tau$')
        ylabel(r'$\eta$')
        xlim(l[0],l[-1])
        #ylim(-1,1)
        # phi.
        fig.add_subplot(223)
        plt.plot(l,[self.phi_tau(_) for _ in l],opts)
        xlabel(r'$\tau$')
        ylabel(r'$\phi$')
        xlim(l[0],l[-1])
        # t.
        fig.add_subplot(224)
        plt.plot(l,[self.t_tau(_) for _ in l],opts)
        xlabel(r'$\tau$')
        ylabel(r'$t$')
        xlim(l[0],l[-1])
        #ylim(-5,5)

def test_random(ntries = 10):
    import random
    random.seed(0)
    retval = []
    for i in range(0,ntries):
        #state = [random.uniform(-10,10) for _ in range(0,6)]
        #state = [random.uniform(-1,1) for _ in range(0,6)]
        state = [random.uniform(-5,5) for _ in range(0,3)] + [random.uniform(-1,1) for _ in range(0,3)]
        a = random.uniform(0,1)
        mu1, mu2 = random.uniform(0,1), random.uniform(0,1)
        e = e3bp(a,mu1,mu2,state)
        if e.bounded:
            t_fin = (10+random.uniform(0,1))*max(e.wp_xi.periods[0].real,e.wp_eta.periods[0].real)*a**2
        else:
            t_fin = e.tau_infs[1] - e.tau_infs[1]/10.
        num = e.numerical_ell_ham(t_fin,t_fin/100.)
        t = num[0][-1]
        delta_phi = abs((e.phi_tau(t) - num[1][-1][5])/num[1][-1][5])
        retval.append((delta_phi,e))
    return sorted(retval,key=lambda _: _[0],reverse=True)

def d_circular_rho(a,mu0,mu1,z):
    from math import sqrt
    A = (mu1*(a+z))**(2./3)
    B = (mu0*(a-z))**(2./3)
    return sqrt(((a+z)**2*B-(a-z)**2*A)/(A-B))

def d_circular_v(a,mu0,mu1,z):
    from math import sqrt
    rho = d_circular_rho(a,mu0,mu1,z)
    Fg = mu0*rho/sqrt(rho**2+(z-a)**2)**3 + mu1*rho/sqrt(rho**2+(z+a)**2)**3
    mu = rho**2*Fg
    v = sqrt(mu/rho)
    return v

# Find a quasi-periodic orbit using Scipy's slsqp algorithm.
# Masses are fixed at +-1 on the z axis, one mass twice the other.
def find_qp(n=3,m=7,delta=.5,max_iter=100,epsilon=1E-12):
    import random, numpy as np
    from scipy.optimize import fmin_slsqp
    # The bounds of the search: around x = 1, y = 0, z = 1, and initial
    # velocity mostly in the y direction.
    min_bounds = [1-delta,-delta,1-delta,-delta,1-delta,-delta]
    max_bounds = [1+delta,delta,1+delta,delta,1+delta,delta]
    bounds = list(zip(min_bounds,max_bounds))
    # Initial guess.
    x0 = np.array([random.uniform(b[0],b[1]) for b in bounds])
    # Objective function.
    def obj_fun(x):
        e = e3bp(1,1.5,1,list(x))
        tmp = e.wp_xi.periods[0].real*m-e.wp_eta.periods[0].real*n
        #return np.array([tmp**2])
        return np.array([abs(tmp)])
    ret = fmin_slsqp(obj_fun,x0,bounds=bounds,acc=1E-15,iter=max_iter,epsilon=epsilon)
    return ret

def find_p(n=3,m=7,k=3,delta=.5,max_iter=100,epsilon=1E-8):
    import random, numpy as np
    from scipy.optimize import fmin_slsqp
    from math import pi, sqrt
    # The bounds of the search: around x = 1, y = 0, z = 1, and initial
    # velocity mostly in the y direction.
    min_bounds = [1-delta,-delta,1-delta,-delta,1-delta,-delta]
    max_bounds = [1+delta,delta,1+delta,delta,1+delta,delta]
    bounds = list(zip(min_bounds,max_bounds))
    # Initial guess.
    x0 = np.array([random.uniform(b[0],b[1]) for b in bounds])
    # Objective function.
    def obj_fun(x):
        e = e3bp(1,1,.05,list(x))
        if not e.bounded:
            return np.array([float('inf')])
        tmp1 = e.wp_xi.periods[0].real*m-e.wp_eta.periods[0].real*n
        P = e.wp_xi.periods[0].real*m
        tmp2 = (e.phi_tau(P) - e.ell_ham_state[5]) % ((2*pi)/k)
        #return np.array([sqrt(tmp1**2+tmp2**2)])
        return np.array([abs(tmp1) + abs(tmp2)])
    ret = fmin_slsqp(obj_fun,x0,bounds=bounds,acc=1E-15,iter=max_iter,epsilon=epsilon)
    return ret

def find_p2(n=3,m=7,k=3,delta=.5,max_iter=100,epsilon=1E-8,x0=None):
    import random, numpy as np
    from scipy.optimize import fmin_slsqp
    from math import pi, sqrt
    # The bounds of the search: around x = 1, y = 0, z = 1, and initial
    # velocity mostly in the y direction.
    min_bounds = [0,-delta,0,-delta,0,-delta]
    max_bounds = [1+delta,delta,1+delta,delta,1+delta,delta]
    bounds = list(zip(min_bounds,max_bounds))
    # Initial guess.
    if x0 is None:
        x0 = np.array([random.uniform(b[0],b[1]) for b in bounds])
    # Objective function.
    def obj_fun(x):
        e = e3bp(1,1.1,2.2,list(x))
        if not e.bounded:
            return np.array([float('inf')])
        tmp1 = e.wp_xi.periods[0].real*m-e.wp_eta.periods[0].real*n
        P = e.wp_xi.periods[0].real*m
        tmp2 = (e.phi_tau(P) - e.ell_ham_state[5]) % ((2*pi)/k)
        return np.array([sqrt(tmp1**2+tmp2**2)])
        return np.array([abs(tmp1) + abs(tmp2)])
    ret = fmin_slsqp(obj_fun,x0,bounds=bounds,acc=1E-15,iter=max_iter,epsilon=epsilon)
    return ret

# A few periodic orbits.
# 1,1.1,2.2,[1.4999990196858124,-0.49397579021723603,1.0989729488631805,-0.49740840558170107,1.4999025672627277,0.49958505894507332] n=91,m=99,k=3
# 1,1,.05,[0.56601403401148631,0.49999957576443482,1.4113903958406369,0.46573660643962256,1.0311053099651799,-0.43669906836064837] n=91,m=99,k=3
# 1,1,.05,[1.2079375966673582,-0.49332055863672453,1.1976067859456487,-0.49843514767491376,0.54822816720530665,0.49662691628363204] n=91,m=99,k=1
# 1,1,.05,[1.0461882287192672,-0.13778723108242488,1.1101217177323248,0.098704778381667607,0.70351429240402208,-0.02212071382670518] n=139,m=141,k=1
# 1,1,1.2,[0.17530743025865267,-0.0212827803021505,1.2124145998267168,-0.10172220543529761,0.096880025404197115,0.2336832259629372] n=139,m=141,k=1
# 1,.05,1,[0.81615762923705693,-0.49987930164344963,0.77532748862554002,0.12141510596757098,0.53086754182698159,0.14622712654555733] n=91,m=99,k=1
# 1,1.1,2.2,[1.499752031671306,0.49954941344236148,1.2764048974966333,-0.10574662666307875,1.4999186167025778,-0.49733228338350099] n=91,m=99,k=1

def find_const_xi(delta=.5,max_iter=100,epsilon=1E-8,x0=None):
    import random, numpy as np
    from scipy.optimize import fmin_slsqp
    from math import pi, sqrt
    # The bounds of the search: around x = 1, y = 0, z = 1, and initial
    # velocity mostly in the y direction.
    min_bounds = [0,-delta,0,-delta,0,-delta]
    max_bounds = [1+delta,delta,1+delta,delta,1+delta,delta]
    bounds = list(zip(min_bounds,max_bounds))
    # Initial guess.
    if x0 is None:
        x0 = np.array([random.uniform(b[0],b[1]) for b in bounds])
    # Objective function.
    def obj_fun(x):
        e = e3bp(1,1.1,2.2,list(x))
        if not e.bounded:
            return np.array([float('inf')])
        return abs(e.fxi[1](e.xi_r))
    ret = fmin_slsqp(obj_fun,x0,bounds=bounds,acc=1E-15,iter=max_iter,epsilon=epsilon)
    return ret

def find_const_eta(delta=.5,max_iter=100,epsilon=1E-8,x0=None):
    import random, numpy as np
    from scipy.optimize import fmin_slsqp
    from math import pi, sqrt
    # The bounds of the search: around x = 1, y = 0, z = 1, and initial
    # velocity mostly in the y direction.
    min_bounds = [0,-delta,0,-delta,0,-delta]
    max_bounds = [1+delta,delta,1+delta,delta,1+delta,delta]
    bounds = list(zip(min_bounds,max_bounds))
    # Initial guess.
    if x0 is None:
        x0 = np.array([random.uniform(b[0],b[1]) for b in bounds])
    # Objective function.
    def obj_fun(x):
        e = e3bp(1,1.1,2.2,list(x))
        if not e.bounded:
            return np.array([float('inf')])
        return abs(e.feta[1](e.eta_r) - .001)
    ret = fmin_slsqp(obj_fun,x0,bounds=bounds,acc=1E-15,iter=max_iter,epsilon=epsilon)
    return ret

def find_const_xi_eta(delta=.5,max_iter=100,epsilon=1E-8,x0=None):
    import random, numpy as np
    from scipy.optimize import minimize
    from math import pi, sqrt
    # The bounds of the search: around x = 1, y = 0, z = 1, and initial
    # velocity mostly in the y direction.
    min_bounds = [0,-delta,0,-delta,0,-delta]
    max_bounds = [1+delta,delta,1+delta,delta,1+delta,delta]
    bounds = list(zip(min_bounds,max_bounds))
    # Initial guess.
    if x0 is None:
        x0 = np.array([random.uniform(b[0],b[1]) for b in bounds])
    # Objective function.
    def obj_fun(x):
        e = e3bp(1,.5,.499,list(x))
        if not e.bounded:
            return np.array([float('inf')])
        return sqrt(e.fxi[1](e.xi_r)**2+e.feta[1](e.eta_r)**2)
    ret = minimize(obj_fun,x0,method="SLSQP",bounds=bounds,tol=1E-15,options={'maxiter':max_iter,'disp':True})
    return ret

# Spring orbit.
# 1,.5,.5,[1.4612799385642152,-0.13596728798792429,0.002542664866964655,0.057450558290888409,0.61745533212591741,9.7341720515643732e-05]
# Kind of sping-like, but different masses.
# 1,.5,.4,[0.028224118013957633,0.45628485533987939,-2.4226329961135072e-08,0.40307849158372189,-1.5050973245126132e-08,-0.075904839286226641]
# These stay in a very small region in the rho/z diagram.
# 1,.5,.499,[0.56824380011393494,-0.0057610562889010945,0.94949708415126532,0.0069907409680348585,0.93261770059829563,-0.00019701094269089766]
# 1,.5,.499,[0.66800950005677395,0.045505046330103549,0.9381421278118599,-0.040558462392587896,0.87954322201140478,-0.01549380118798635]
# This one is basically a circular orbit displaced on the z axis (as the ones above as well, probably).
# 1,.5,.499,[0.90005170183260064,0.028752907817175204,0.84047543280198744,-0.023874116354788769,0.76044103757279002,0.0014508691116051184]

def find_const_xi_eta_r(delta=.5,max_iter=100,epsilon=1E-8,x0=None,n=100):
    retval = []
    for _ in range(n):
        try:
            ret = find_const_xi_eta(delta,max_iter,epsilon,x0)
            retval.append(ret)
        except:
            raise
    return retval

#from PyGMO.problem import base as _pbase
#class qp_prob(_pbase):
    #def __init__(self,n=3,m=7):
        #super(qp_prob,self).__init__(6)
        #self.set_bounds(-2,2)
        #self.n = n
        #self.m = m
    #def _objfun_impl(self,x):
        #e = e3bp(1,2,1,list(x))
        #tmp = e.wp_xi.periods[0].real*self.m-e.wp_eta.periods[0].real*self.n
        #return (tmp**2,)

#class p_prob(_pbase):
    #def __init__(self,n=3,m=7,k=3):
        #super(p_prob,self).__init__(6)
        #delta = .5
        #self.set_bounds([1-delta,-delta,1-delta,-delta,1-delta,-delta],[1+delta,delta,1+delta,delta,1+delta,delta])
        #self.n = n
        #self.m = m
        #self.k = k
    #def _objfun_impl(self,x):
        #from math import pi, sqrt
        #e = e3bp(1,1,.05,list(x))
        #tmp1 = e.wp_xi.periods[0].real*self.m-e.wp_eta.periods[0].real*self.n
        #P = e.wp_xi.periods[0].real*self.m
        #tmp2 = (e.phi_tau(P) - e.ell_ham_state[5]) % ((2*pi)/self.k)
        #return (sqrt(tmp1**2+tmp2**2),)

#def find_quasi_periodic(n,m):
    #from PyGMO.algorithm import scipy_slsqp
    #from PyGMO import island
    #prob = qp_prob(n,m)
    #algo = scipy_slsqp()
    #isl = island(algo,prob,20)
    #pop = algo.evolve(isl.population)
    #pop = algo.evolve(pop)
    #return pop.champion.f,pop.champion

#def find_periodic(n,m,k):
    #from PyGMO.algorithm import scipy_slsqp
    #from PyGMO import island
    #prob = p_prob(n,m,k)
    #algo = scipy_slsqp()
    #isl = island(algo,prob,20)
    #pop = algo.evolve(isl.population)
    #pop = algo.evolve(pop)
    #pop = algo.evolve(pop)
    #pop = algo.evolve(pop)
    #pop = algo.evolve(pop)
    #return pop.champion.f,pop.champion

#(0.5975967480103391,-1.3338365760059723,0.8720813209186803,1.5476676696694698,0.28450034681617586,2.8000268923310587) 1 1 1
#(-0.9471035756992736,-0.5646627652559865,-0.6977162142941262,-0.7835596640482966,-1.086049384461206,0.03826612625616655) 1 1 1
#(0.5128698969888743,-0.4996325712869821,0.5022560614711243,-0.4999955600124257,1.3877141085192108,0.27792642325022415) 1 1 0.1
#(0.6025587736272164,-0.3794756224426397,1.251129375638713,-0.46980563205199743,1.4738485324711814,-0.4697756558851869) 1 1 0.1
#(0.7350903691315228,0.25397998840230623,0.5850432509821148,-0.43140559115358235,1.3195748725913112,-0.47529674389957627) 1 1 0.1 91 99
#(1.1142692720185303,-0.4664080662038235,1.4409265564824143,-0.468268650165366,0.5793512521216965,0.4975713849386918) 1 1 0.05 91 99
#(0.5964967500398608,-0.46658969392422284,1.2199091480279232,-0.21345640283724415,1.0745552544242645,-0.08189741866819046) 1 1 0.05, 139 141
#(0.5080708298333626,0.4974656543529942,0.6910639145170152,0.4002982142305068,0.9449717205577095,0.4997852434425768) 1 1 0.05, 139 141
#(1.258686479754558,-0.4201802370015161,1.2543488848121276,-0.49913249267287796,0.5002564030871568,0.49994544721370154) 1 1 0.05, 91 99
