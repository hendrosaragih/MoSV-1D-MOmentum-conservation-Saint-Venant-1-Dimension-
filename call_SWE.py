import numpy as np
from initial_Condition import *
def call_SWE(xleft, xright, Nx, dx, tfin, A0, xC, xwide, hC, h0, hshall, hwide) :
	# 0. INPUT PARAMETERS
	# wave initial condition
	eta_list = list()
	time_list = list()
	x  = [0]*(Nx+1)
	x  = np.linspace(xleft, xright, Nx + 1)
	
	# bathymetry
	# hC      = 750.
	# h0      = 20.
	# hshall  = 1.
	# hwide   = 100.
	# time
	dt      = 0.0001
	# parameters
	g       = 9.81

	# 1. PREPARATION 
	h  		= [0]*(Nx+1)
	H  		= [0]*(Nx+1)	
	H_star  = [0]*(Nx+2)  # total depth in 'half-grid'
	H_bar  	= [0]*(Nx+2)  # total depth at half grid
	u 		= [0]*(Nx+2)  # velocity in 'half-grid'
	u_star  = [0]*(Nx+1)
	u_du 	= [0]*(Nx+2)  # advection term in half-grid
	u_new  	= [0]*(Nx+2)  # updated velocity u
	q  		= [0]*(Nx+2)  # q flux = H*u at half grid
	q_bar  	= [0]*(Nx+1)  # q at full grid
	eta  	= [0]*(Nx+1)  # old eta
	eta_new = [0]*(Nx+1)  # updated eta

	eta = initial_Condition(x, A0, xC, xwide)
	h = h0 - hshall*np.exp(-((x-hC)/hwide)*((x-hC)/hwide)) # 2.2. Bathymetry/bottom
	H = h + eta                                            # 2.3. total depth in full grid

	t = 0
	step = 0

	while (t <= tfin) :
		t = t + dt
		step += 1
		for j in range(1, Nx+1) :
			if (u[j] >= 0) :
				H_star[j]     = H[j-1]
				# alpha_star[j] = alpha[j-1]
				# beta_star[j]  = beta[j-1]
			else :
				H_star[j]     = H[j]
				# alpha_star[j] = alpha[j]
				# beta_star[j]  = beta[j]
			q[j] = H_star[j]*u[j]
		H_star[0]       = H_star[1]
		H_star[Nx+1]    = H_star[Nx]
		q[0]    = 0
		q[Nx+1] = 0


		for j in range(0, Nx+1) :
			if(j==0 or j==Nx) :
				eta_new[j] = eta[j] - (dt/dx)*(q[j+1]-q[j])
			else :
				eta_new[j] = eta[j] - (dt/dx)*(q[j+1]-q[j])
			H[j] = eta_new[j] + h[j]
			q_bar[j] = 0.5*(q[j+1] + q[j])
			if(q_bar[j] >= 0) :
				u_star[j]=u[j]
			else :
				u_star[j]=u[j+1]
		# 3.2. Calculate 2nd eq. (momentum eq.)
		for j in range(1, Nx+1) :
			H_bar[j] = 0.5*(H[j] + H[j-1])
			if(H_bar[j] == 0) :
				u_new[j] = 0
			else :
				u_du[j] = 1./(H_bar[j]*dx)*((q_bar[j]*u_star[j] - q_bar[j-1]*u_star[j-1]) - (u[j]*(q_bar[j] - q_bar[j-1])))
				u_new[j] = u[j] - (dt/dx*g*( eta_new[j] - eta_new[j-1] )) - (dt*u_du[j])
		u_new[0] = 0
		u_new[Nx+1] = 0

		# updating values of H, eta, u
		eta   = eta_new
		H     = h + eta
		u     = u_new
		if (step % 2 == 0) :
			eta_list.append(eta_new[:])
			time_list.append(t)
	return eta_list
#--------------------------------------------------------------------------------------------------------------------------------------