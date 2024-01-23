from __future__ import print_function
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.animation as animation
from mshr import *
from PIL import Image, ImageSequence

def Helmholtz(k, f, typeOfBoundary, g, u_D):
    n = mesh.num_vertices()
    d = mesh.geometry().dim()
    mesh_coordinates = mesh.coordinates().reshape((n, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0], mesh_coordinates[:, 1], triangles)
    if typeOfBoundary == 2:
        def boundary(x, on_boundary):
            return on_boundary and (near(x[0], 0, tol) or near(x[0],1, tol))
    else:
        def boundary(x, on_boundary):
            return on_boundary
#    if k <= 0:
#        P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
#        R = FiniteElement("Real", mesh.ufl_cell(), 0)
#        W = FunctionSpace(mesh, P1 * R)
#        # Define variational problem
#        (u, c) = TrialFunction(W)
#        (v, d) = TestFunctions(W)
#        #f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
#        #g = Expression("-sin(5*x[0])", degree=2)
#        a = (inner(grad(u), grad(v)) + c*v + u*d)*dx
#        L = f*v*dx + g*v*ds

#        # Compute solution
#        w = Function(W)
#        solve(a == L, w)
#        (u, c) = w.split()

#        # Plot solution
#        #plot(u, interactive=True)
#        z = np.asarray([u(point) for point in mesh_coordinates])
#        plt.tripcolor(triangulation, z, edgecolors='k')
#        plt.savefig('u_vertex.png')
#        plt.colorbar()
#        error_L2 = errornorm(u_D, u, 'L2')
#        print('L2-error = ', error_L2)
#        return

    #function space
    V = FunctionSpace(mesh, 'P', 1)

    bc = DirichletBC(V, u_D, boundary)

    # defining variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
         
    a = dot(grad(u), grad(v))*dx + k*u*v*dx
    if typeOfBoundary == 2:
        L = f*v*dx + g*v*ds
    else:
        L = f*v*dx

    # computing solution
    u = Function(V)
    solve(a == L, u, bc)
    
#    plot(u, interactive=True)
#    
#    vtkfile = File('1_solution.pvd')
#    vtkfile << u
#    
    plt.figure()
    z = np.asarray([u(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, edgecolors='k')
    plt.savefig('u_vertex.png')
    plt.colorbar()

    # computing L2-error
    error_L2 = errornorm(u_D, u, 'L2')

    # printing errors
    print('L2-error = ', error_L2)
    
    return u

    
    
def Heat(aa, f, typeOfBoundary, g, u_D):
    n = mesh.num_vertices()
    d = mesh.geometry().dim()
    mesh_coordinates = mesh.coordinates().reshape((n, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0], mesh_coordinates[:, 1], triangles)
    V = FunctionSpace(mesh, 'P', 1)
    u_n = interpolate(u_D, V)
    t = 0
    vtkfile = File('2_solution.pvd')
    imlist = []
    fig = plt.figure()
    for n in range(num_steps):
        # updating current time
        t += dt
        u_D.t = t
        u = Helmholtz(1/(aa*dt), (u_n/(aa*dt) + (1/aa)*f), typeOfBoundary, g/aa, u_D)
        #plot(u)
        vtkfile << u
        u_e = interpolate(u_D, V)
        z = np.asarray([u(point) for point in mesh_coordinates])
        im = plt.tripcolor(triangulation, z, edgecolors='k', vmin=0, vmax = 5)
        imlist.append([im])
        # computing L2-error
        error_L2 = errornorm(u_e, u, 'L2')
        print(' t = ', t, ', error = ', error_L2)
        u_n.assign(u)
    an = animation.ArtistAnimation(fig, imlist, interval=50, blit=True, repeat_delay=300)
    
    an.save('line1.gif', dpi=80, writer='imagemagick')

    #interactive()

    
tol = 1E-14
# Creating mesh
domain = Circle(Point(0.0,0.0),1.0)
mesh = generate_mesh(domain, 64)
n = mesh.num_vertices()
d = mesh.geometry().dim()

T = 2.0
num_steps = 100
dt = T / num_steps
alpha = 3
beta = 1.2

#Helmholtz(0, Expression('x[0]', degree=2), 2, Expression('cos(atan2(x[1],x[0]))', degree=2), Expression('x[0]', degree=2))
#Helmholtz(0, Expression('0', degree=2), 1, 0, Expression('x[0] + x[1]', degree=2))
#Helmholtz(1, Expression('x[0] + x[1]', degree=2), 1, 0, Expression('x[0] + x[1]', degree=2))
#Helmholtz(1, Expression('x[0]', degree=2), 1, 0, Expression('x[0]', degree=2))
#Helmholtz(1, Expression('x[0] + x[1]', degree=2), 2, Expression('cos(atan2(x[1],x[0])) + sin(atan2(x[1],x[0]))', degree=2), Expression('x[0] + x[1]', degree=2))
#Helmholtz(1, Expression('x[0]', degree=2), 2, Expression('cos(atan2(x[1],x[0]))', degree=2), Expression('x[0]', degree=2))
#Heat(1, Constant(beta - 2 - 2*alpha), 1, 0, Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0))
#Heat(1, Expression('beta - (6*x[0] + x[0]*x[0]*x[0] + 6*alpha*x[1] + alpha*x[1]*x[1]*x[1])', degree=2, alpha=alpha, beta=beta), 1, 0, Expression('x[0]*x[0]*x[0] + alpha*x[1]*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0))
#Heat(1, Expression('beta - x[0] - x[1]', degree=2, beta = beta), 2, Expression('cos(atan2(x[1],x[0])) + sin(atan2(x[1],x[0]))', degree=2, alpha=alpha), Expression('x[0] + x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0))
Heat(1, Expression('beta - x[0]', degree=2, alpha=alpha, beta=beta), 2, Expression('cos(atan2(x[1],x[0]))', degree=2), Expression('x[0] + beta*t', degree=2, alpha=alpha, beta=beta, t=0))

