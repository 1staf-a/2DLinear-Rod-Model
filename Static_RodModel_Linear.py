import numpy as np
import matplotlib.pyplot as plt

# limiting this to only work for cantilever beam with force at one end

class Rod: #2D cylindrical rod
    """
    Basic 2D rod for simulating the linear elastic model. Can also be used for non-linear model.
    """

    m = 1                                       # mass per unit length (not used here)

    L = 1.0                                     # default rod length to 1[m]
    dx = 0.1                                    # minimum distance between cross-sections [m]
    n = L / dx                                  # number of sections
    x = np.arange(0, L + dx, dx)                # parametric length [m]

    I = 1                                       # moment of inertia of cross-sections [kg mm^2]
    E = 1                                       # elastic modulus [N/mm2]
    dy = np.zeros_like(x)                       # deflection of cross-sections in [m]

    kappa = np.zeros_like(x)                    # curvature of each cross-section (ignore for now)
    r = [x,dy]                                  # position vector of cross-section centreline [m,m]

    f = np.ones_like(x)                         # constant distributed load [N/m] (Points +j)

    F=0                                         # force at free end [N] (Points -j)
    Q=0                                         # moment at free end [Nm] (Points +k)

    V = np.zeros_like(x)                        # internal shear loading [N]
    M = np.zeros_like(x)                        # internal moment [Nm]

def boundary_condition(rod_i,load_mag=1,load_var=0,Conc_load=0,Conc_moment=0):
    """
    Apply the boundary conditions to the rod
    :param rod_i: provide a rod class object
    :param load_mag: scalar multiple(N/mm) for the distributed load being applied.
    :param load_var: enter the order of the variation of the loads wrt length. zero meaning constant, 1 meaning linear variation and so on...
    :param Conc_load: concentrated load(N) being applied at the free end
    :param Conc_moment: concentrated moment(Nm) being applied at the free end
    :return:
    """
    rod_i.F=Conc_load
    rod_i.Q=Conc_moment
    rod_i.f= load_mag * rod_i.x ** load_var

    pass

def rod_properties(rod_i,elasticity=1,MoI=1,L=1,dx=0.1):
    """
    Change the properties of the rod_i
    :param rod_i: provide a rod class object to be modified
    :param elasticity: new elasticity
    :param MoI: new Moment of Inertia of cross-sections
    :param L: new Length of Rod
    :param dx: distance of rod segments
    :return:
    """
    # functionality that can be added later
    # d = 2.82842712474619                      # setting diameter(mm) to moment of inertia to 1mm^4
    # I = 1 / 2 * m * d * d / 4                   # moments of inertia of each cross-section
    # m = np.ones(len(x))                         # mass per unit length


    rod_i.E=elasticity
    rod_i.I=MoI
    rod_i.L=L
    rod_i.dx=dx
    rod_i.x=np.arange(0, rod_i.L + rod_i.dx, rod_i.dx)

def linear_model(rod_i):
    """
    Solves the linear elastic rod model for a perpendicular distributed(force) or concentrated load(at free end: force or moment)
    :param rod_i: inputs the rod object which has all the defined loads applied to it.
    :return: the deformed rod object.
    """
    rod_i.dy=1/(rod_i.E*rod_i.I)*(rod_i.f * rod_i.x ** 4 / 24 + (rod_i.F - rod_i.f * rod_i.L) * rod_i.x ** 3 / 6 +
                                  (rod_i.Q - rod_i.F * rod_i.L + rod_i.f * rod_i.L ** 2 / 2) * rod_i.x ** 2 / 2)
    rod_i.r = [rod_i.x, rod_i.dy]
    rod_i.V = rod_i.f*rod_i.x
    rod_i.V[-1] -= rod_i.F
    rod_i.M = rod_i.f*rod_i.x**2/2
    rod_i.M[-1] -= rod_i.Q
    pass

def plot_rod(rod_i):
    fig=plt.figure()
    fig.set_size_inches(10,10)
    ax1 = plt.subplot(3,1,1)

    ax1.grid()

    plt.title("Rod Under Loading")
    max_def = np.amax(np.abs(rod_i.dy))
    plt.scatter(rod_i.x, rod_i.dy)
    plt.plot(rod_i.x, rod_i.dy)
    plt.tick_params('x', labelbottom=False)
    plt.ylabel('Y [m]')
    plt.ylim(max_def*-2, max_def*2)

    ax2 = plt.subplot(312, sharex=ax1)
    ax2.grid()
    plt.plot(rod_i.x,rod_i.V)
    plt.ylim(rod_i.V[0], rod_i.V[-1])
    plt.tick_params('x', labelbottom=False)
    plt.ylabel('Shear [N]')

    ax3 = plt.subplot(313, sharex=ax1)
    ax3.grid()
    plt.plot(rod_i.x,rod_i.M)
    plt.ylim(rod_i.M[0], rod_i.M[-1])
    plt.ylabel('Moment [Nm]')
    plt.xlabel('X [m]')
    plt.show()

if __name__ == '__main__':
    rod = Rod()
    #Vary only the parameters in this main function if you wish to simulate a different environment
    boundary_condition(rod,load_mag=1,load_var=0,Conc_load=0,Conc_moment=0)
    rod_properties(rod,elasticity=1,MoI=1,L=1,dx=0.1)
    linear_model(rod)
    plot_rod(rod)






