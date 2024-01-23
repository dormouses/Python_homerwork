from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from dolfin import *
import numpy as np
from mshr import *
from PIL import Image, ImageSequence
import matplotlib.animation as animation


def SolvH(alpha, f, type, g, u_D):
    n = mesh.num_vertices()
    d = mesh.geometry().dim()
    mesh_coordinates = mesh.coordinates().reshape((n, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0], mesh_coordinates[:, 1], triangles)
    
    if type == 1:
        def boundary(x, on_boundary):
            return on_boundary
    else:
        def boundary(x, on_boundary):
            return on_boundary and (near(x[0], 0, 1E-10) or near(x[0],1, 1E-10))

    #    if alpha == 0:
    #        P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    #        R = FiniteElement("Real", mesh.ufl_cell(), 0)
    #        W = FunctionSpace(mesh, P1 * R)
    #        # Define variational problem
    #        (u, c) = TrialFunction(W)
    #        (v, d) = TestFunctions(W)
    #        a = (inner(grad(u), grad(v)) + c*v + u*d)*dx
    #        L = f*v*dx + g*v*ds
    
    #        # Compute solution
    #        w = Function(W)
    #        solve(a == L, w)
    #        (u, c) = w.split()
    
    #        # Plot solution
    #        #plot(u, interactive=True)
    #        z = np.asarray([u(point) for point in mesh_coordinates])
    #        plt.tripcolor(triangulation, z, edgecolors='alpha')
    #        plt.savefig('u_vertex.png')
    #        plt.colorbar()
    #        error_L2 = errornorm(u_D, u, 'L2')
    #        return error_L2
    
    #function space
    V = FunctionSpace(mesh, 'P', 1)
    bc = DirichletBC(V, u_D, boundary)
    
    # defining variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    
    a = dot(grad(u), grad(v))*dx + alpha*u*v*dx
    if type == 2:
        L = f*v*dx + g*v*ds
    else:
        L = f*v*dx

    # computing solution
    u = Function(V)
    solve(a == L, u, bc)

    #######

    plt.figure()
    z = np.asarray([u(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, edgecolors='k')
    plt.savefig('u_vertex.png')
    plt.colorbar()
    
    # computing L2-error
    error_L2 = errornorm(u_D, u, 'L2')
    #print (error_L2)
    return (error_L2, u)
    #return

def Heat(alpha, f, type, g, u_D):
    dt = 1.0 / 10
    n = mesh.num_vertices()
    d = mesh.geometry().dim()
    mesh_coordinates = mesh.coordinates().reshape((n, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0], mesh_coordinates[:, 1], triangles)
    V = FunctionSpace(mesh, 'P', 1)
    u_n = interpolate(u_D, V)
    t = 0
    ######
    imlist = []
    fig = plt.figure()
    er=0
    for n in range(10):
        # updating current time
        t += dt
        u_D.t = t
        u = SolvH(1/(alpha*dt), (u_n/(alpha*dt) + (1/alpha)*f), type, g/alpha, u_D)[1]
        #####
        u_e = interpolate(u_D, V)
        z = np.asarray([u(point) for point in mesh_coordinates])
        im = plt.tripcolor(triangulation, z, edgecolors='k', vmin=0, vmax = 5)
        imlist.append([im])
        # computing L2-error
        error_L2 = errornorm(u_e, u, 'L2')
        print ('t=', t, 'er=',error_L2 )
        er+=error_L2
        u_n.assign(u)
    #ans = animation.ArtistAnimation(fig, imlist, interval=50, blit=True, repeat_delay=300)
    #ans.save('ans.gif', dpi=80, writer='imagemagick')
    #ans.save('ans.mp4', dpi=80, writer='ffmpeg')
    return (er, u)

# Creating mesh
domain = Circle(Point(0.0,0.0),1.0)
mesh = generate_mesh(domain, 64)


#print (SolvH(1, Expression('x[0] + x[1]', degree=2), 1, Expression('cos(atan2(x[1],x[0]))+sin(atan2(x[1],x[0]))', degree=2), Expression('x[0] + x[1]', degree=2)) [0] )
#print (SolvH(1, Expression('2 * sin (x[0])', degree=2), 2, Expression('cos(atan2(x[1],x[0]))*cos(x[0])', degree=2), Expression('sin(x[0])', degree=2))[0] )
#print (SolvH(1, Expression('x[0]', degree=2), 1, 0, Expression('x[0]', degree=2))[0] )
#print (SolvH(1, Expression('2 * cos (x[0])', degree=2), 2, Expression('-cos(atan2(x[1],x[0]))*sin(x[0])', degree=2), Expression('cos(x[0])', degree=2))[0] )


#print (Heat(1, Constant(1), 1, Expression('2*x[0]*cos(atan2(x[1],x[0]))+2*x[1]*sin(atan2(x[1],x[0]))', degree=2), Expression('x[0]*x[0] + x[1]*x[1] + t', degree=2, t=0)) [0])
print (Heat(1, Constant(-8), 1, Expression('6*x[0]+6*x[1]', degree=2), Expression('3*x[0]*x[0] + 3*x[1]*x[1] + 4*t', degree=2, t=0))[0] )
#print (Heat(2, Expression('2*sin(x[0])', degree=2), 2, Expression('cos(atan2(x[1],x[0]))*cos(x[0])', degree=2, alpha=alpha), Expression('sin(x[0])', degree=2, t=0))[0] )
#print  (Heat(3, Expression('3*cos(x[0])', degree=2), 2, Expression('sin(atan2(x[1],x[0]))*sin(x[1])', degree=2, alpha=alpha), Expression('cos(x[1])', degree=2, t=0)) [0])


