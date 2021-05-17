from sympy import *
from sympy.matrices import Matrix

# x = Symbol('x')
# y = Symbol('y')

x0 = Symbol('x0')
y0 = Symbol('y0')
x1 = Symbol('x1')
y1 = Symbol('y1')
x2 = Symbol('x2')
y2 = Symbol('y2')
x3 = Symbol('x3')
y3 = Symbol('y3')

xi = Symbol('xi') #chi
eta = Symbol('eta')

#shape functions
N0 = (1-xi)*(1-eta)/4	
N1 = (1+xi)*(1-eta)/4	
N2 = (1+xi)*(1+eta)/4	
N3 = (1-xi)*(1+eta)/4	

#transformacion de coordenadas, va de xi, eta a coordenadas x,y
x = x0 * N0 + x1 * N1 + x2 * N2 + x3 * N3
y = y0 * N0 + y1 * N1 + y2 * N2 + y3 * N3

#Matriz jacobiana
dx_dxi = x.diff(xi)
dx_deta = x.diff(eta)
dy_dxi = y.diff(xi)
dy_deta = y.diff(eta)

J = Matrix([
	[dx_dxi, dx_deta],
	[dy_dxi, dy_deta]
	])

#print(pretty(J))

#Jinv = J.inv()
#print(pretty(Jinv))

detJ = simplify(J.det())
print(pretty(detJ))