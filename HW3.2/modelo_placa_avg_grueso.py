import numpy as np
from quad4_avg import quad4, quad4_post, quad_line_load

fid = open("MODELO 2D - grueso.msh", "r")

LINE_ELEMENT = 1
TRI_ELEMENTO = 2
QUAD_ELEMENT = 3

Empotrado = 16
BordeNatural = 17 
Placa = 14
Extremos = 15

#geometric parameters
a = 2;
b = 4;


while True: 
	line = fid.readline()

	if line.find("$Nodes") >= 0:
		break 

Nnodes = int(fid.readline())

xy = np.zeros([Nnodes, 2])

for i in range(Nnodes):
	line = fid.readline()
	sl = line.split()
	xy[i,0] = float(sl[1])
	xy[i,1] = float(sl[2])

print(f" xy = {xy}")

print(Nnodes)

while True: 
	line = fid.readline()

	if line.find("$Elements") >= 0:
		break 

Nelements = int(fid.readline())

print(Nelements)

conec = np.zeros((Nelements, 4), dtype=np.int32)

fixed_nodes = []
nodos_borde_nat = []
Nquads = 0
Quadrangles = []

for i in range(Nelements):
	line = fid.readline()
	sl = line.split()
	element_number = np.int32(sl[0]) - 1
	element_type = np.int32(sl[1])
	physical_grp = np.int32(sl[3])
	entity_number = np.int32(sl[4])

	if element_type == LINE_ELEMENT and \
		physical_grp == Empotrado: #linea
		n1 = np.int32(sl[5]) - 1
		n2 = np.int32(sl[6]) - 1
		fixed_nodes += [n1, n2]
    
	if element_type == LINE_ELEMENT and \
		physical_grp == BordeNatural: #linea
		n1 = np.int32(sl[5]) - 1
		n2 = np.int32(sl[6]) - 1
		nodos_borde_nat.append([n1, n2])

	if element_type == QUAD_ELEMENT and \
		(physical_grp == Placa or physical_grp == Extremos):
		n0 = np.int32(sl[5]) - 1
		n1 = np.int32(sl[6]) - 1
		n2 = np.int32(sl[7]) - 1
		n3 = np.int32(sl[8]) - 1


		conec[element_number, :] = [n0, n1, n2, n3]

		Quadrangles.append(element_number)
		Nquads += 1

print(conec[15])

print(conec) #NOTAR QUE EL CONEC SALE CON PUROS 0 PERO DENTRO DEL FOR SI TIENE LOS DISTINTOS VALORES
#print(Quadrangles)


print("Fin del archivo")

NDOFs = 2*Nnodes

properties = {}

rho = 2500.
g = 9.81

properties = {}
properties["E"] = 20e9 #del profe, por mientras
properties["nu"] = 0.25
properties["bx"] = 0.
properties["by"] = 0.
properties["t"] = 4e-3 #4 mm

K = np.zeros((NDOFs ,NDOFs))
f = np.zeros((NDOFs, 1))

for e in Quadrangles: 

	ni = conec[e, 0]
	nj = conec[e, 1]
	nk = conec[e, 2]
	nl = conec[e, 3]

	print(f"e = {e}, ni = {ni}, nj = {nj}, nk = {nk}, nl = {nl}")

	xy_e = xy[[ni, nj, nk, nl], :]

	ke, fe = quad4(xy_e, properties)

	#print (f"ke = {ke}")

	# Node k --> [3*k, 3*k+1, 3*k +2]

	d = [2*ni, 2*ni+1, 2*nj, 2*nj+1, 2*nk, 2*nk+1, 2*nl, 2*nl+1] #global dofs from local axes

	#DSM
	for i in range(8):
		p = d[i]
		for j in range(8):
			q = d[j]
			K[p, q] += ke[i,j]
		#f[p] += fe[i]

properties_load = {}
properties_load["t"] = properties["t"]
properties_load["tx"] = 1000/(properties_load["t"]*b)
properties_load["ty"] = 0


for nn in nodos_borde_nat:
    ni = nn[0]
    nj = nn[1]
    
    xy_e = xy[[ni, nj], :]
    print(f"ni = {ni} nj = {nj} xy_e = {xy_e}")
    
    ke, fe = quad_line_load(xy_e, properties_load)
    
    d = [2*ni, 2*ni+1, 2*nj, 2*nj+1] #globalDOFs from local DOFs
    
    #DSM
    for i in range(4):
        p = d[i]
        f[p] += fe[i]
        
print(f"= {f}")


fixed_nodes = np.unique(fixed_nodes)

constrained_DOFs = []

for n in fixed_nodes:
	constrained_DOFs += [2*n, 2*n+1]

free_DOFs = np.arange(NDOFs)
free_DOFs = np.setdiff1d(free_DOFs, constrained_DOFs)

print(f"fixed_nodes = {fixed_nodes}")
print(f"constrained_DOFs = {constrained_DOFs}")
print(f"free_DOFs = {free_DOFs}")


#print(K)
#plt.matshow(K)
#plt.show()

for i in range(NDOFs):
    if abs(K[i,i])<1e-8:
        print(f"i={i} K[i,i] = {K[i,i]}")
    
    
Kff = K[np.ix_(free_DOFs, free_DOFs)]
Kfc = K[np.ix_(free_DOFs, constrained_DOFs)]
Kcf = K[np.ix_(constrained_DOFs, free_DOFs)]
Kcc = K[np.ix_(constrained_DOFs, constrained_DOFs)]

ff = f[free_DOFs]
fc = f[constrained_DOFs]
print(f"ff = {ff}")


# Solve
from scipy.linalg import solve
u = np.zeros((NDOFs,1))

u[free_DOFs] = solve(Kff, ff)

#Get reaction forces 
R = Kcf @ u[free_DOFs] + Kcc @ u[constrained_DOFs]	- fc 

print (f"u = {u}")
print (f"R = {R}")


import matplotlib.pyplot as plt

factor = 1e5

uv= u.reshape([-1, 2])  #*factor
print(f"uv = {uv}")
plt.plot(xy[:,0] + factor*uv[:,0], xy[:,1]+ factor*uv[:,1], '.')

for e in Quadrangles:  #ingresa cada elemento e con su respecitov ke, fe y lo ensambla con el DSM
	ni = conec[e,0]
	nj = conec[e,1]
	nk = conec[e,2]
	nl = conec[e,3]
    
	xy_e = xy[[ni, nj, nk, nl, ni],:] + factor*uv[[ni, nj, nk, nl, ni],:]  #esto para muchos elementos se va a trabar en matplotlib, visualizar en gmsh

	plt.plot(xy_e[:,0], xy_e[:,1], 'k')

plt.axis('equal')
plt.show()


from gmsh_post import write_node_data, write_node_data_2, write_element_data

nodes = np.arange(1, Nnodes+1)

write_node_data("ux_avg_grueso.msh", nodes, uv[:,0], "Despl. X" )
write_node_data("uy_avg_grueso.msh", nodes, uv[:,1], "Despl. Y" )

write_node_data_2("desplazamientos_avg_grueso.msh", nodes, uv[:,0], uv[:,1], "Despl" )

#calculo de tensiones

sigmaxx = np.zeros(Nquads+1)
sigmayy = np.zeros(Nquads+1)
sigmaxy = np.zeros(Nquads+1)

i = 0
for e in Quadrangles:  #ingresa cada elemento e con su respecitov ke, fe y lo ensambla con el DSM
	ni = conec[e,0]
	nj = conec[e,1]
	nk = conec[e,2]
	nl = conec[e,3]
	xy_e = xy[[ni, nj, nk, nl, ni],:] 

	uv_e = uv[[ni, nj, nk, nl],:]  #esto para muchos elementos se va a trabar en matplotlib, visualizar en gmsh

	u_e = uv_e.reshape((-1)) 

	Ee,  sigmae = quad4_post(xy_e, u_e, properties)

	sigmaxx[i] = sigmae[0] 
	sigmayy[i] = sigmae[1]
	sigmaxy[i] = sigmae[2]

	i += 1

	#plt.plot(xy_e[:,0], xy_e[:,1], 'k')


elements = np.array(Quadrangles)+1
write_element_data("sigma_x_avg_grueso.msh", elements, sigmaxx, "Sigma_x" )
write_element_data("sigma_y_avg_grueso.msh", elements, sigmayy, "Sigma_y" )
write_element_data("sigma_xy_avg_grueso.msh", elements, sigmaxy, "Sigma_xy" )