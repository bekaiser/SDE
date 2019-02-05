# Function bank
# Bryan Kaiser
# 10/16/16


# =============================================================================
# function declaration 

function advection(adv,q,Kc,Lc,u,v,Xshift,Yshift)
	Xshift[:,2:Nx] = r2r(u.*q,FFTW.RODFT10,2)[:,1:Nx-1]./(2.0*float(Nx)) 
	# fft on Ny number of rows
	Yshift[2:Ny,:] = r2r(v.*q,FFTW.RODFT10,1)[1:Ny-1,:]./(2.0*float(Ny)) 
	# fft on Nx number of columns
	adv[:,:] = r2r(Xshift.*Kc,FFTW.REDFT01,2)+r2r(Yshift.*Lc,FFTW.REDFT01,1)
	# flux conserving advection = Jacobian(psi,q)
	return adv
end

function inversion(psi,PSI,Q,Ks_2,Ls_2,LR_m2)
	PSI[:,:] = -Q.*((Ks_2+Ls_2+LR_m2).^(-1))
	psi[:,:] = r2r(PSI,FFTW.RODFT01,1:2) # DST-III, inverse
	return psi,PSI
end

function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
    	m, n = length(vy), length(vx)
    	vx = reshape(vx, 1, n)
    	vy = reshape(vy, m, 1)
    	(repmat(vx, m, 1), repmat(vy, 1, n))
end

function RKcoeff(q,Q,psi,PSI,k,Ks,Ls,Kc,Lc,LR_m2,beta,kappa,r,Xshift,Yshift,u,v,plnt,adv,diff,btmfr,hypdiff,phi)
	Q[:,:] = r2r(q,FFTW.RODFT10,1:2)./(4.0*float(Nx)*float(Ny)) # DST-II
	psi[:,:],PSI[:,:] = inversion(psi,PSI,Q,Ks_2,Ls_2,LR_m2)
	u[:,:],v[:,:],plnt[:,:] = velocities_beta(u,v,plnt,psi,Kc,Lc,beta,Xshift,Yshift)
	adv[:,:] = advection(adv,q,Kc,Lc,u,v,Xshift,Yshift) 
	# advection, planetary vorticity (1st derivatives) 
	diff[:,:] = r2r(-Q.*(Ks_2+Ls_2),FFTW.RODFT01,1:2).*kappa 
	# diffusion of qgpv, kappa must be in units of m^2/s
	hypdiff[:,:] = r2r(-Q.*(Ks_6+Ls_6),FFTW.RODFT01,1:2).*phi 
	# diffusion of potential vorticity
	if r > 0.0
	btmfr[:,:] = r2r(-PSI.*(Ks_2+Ls_2),FFTW.RODFT01,1:2).*r # bottom friction
	k[:,:] = -adv-plnt+diff-btmfr+wind+hypdiff # Runge-Kutta coefficient
	else
	k[:,:] = -adv-plnt+diff+wind+hypdiff # Runge-Kutta coefficient
	end
	return k
end

function sde_rk4(T,z_dt,A)
	# RK4 k coefficient generator. Uses the algorithm
	# described by Hansen & Penland (2006).

	# Here, the spatial correlation matrix G is 
	# integrated for uncorrelated, autonomous
	# spatial correlations. The dynamical operator
	# A is also autonomous, so A*T is autonomous.

	# input variables:
	# size(T) = [1 N], N = number of state variables
	# t = a scalar variable; the time
	# dt = the time step.

	k = A*T+z_dt;
	return k
end

function velocities_beta(u,v,plnt,psi,Kc,Lc,beta,Xshift,Yshift)
	Xshift[:,2:Nx] = r2r(psi,FFTW.RODFT10,2)[:,1:Nx-1]./(2.0*float(Nx)) 
	# fft on Ny number of rows
	Yshift[2:Ny,:] = r2r(psi,FFTW.RODFT10,1)[1:Ny-1,:]./(2.0*float(Ny)) 
	# fft on Nx number of columns
	v[:,:] = r2r(Xshift.*Kc,FFTW.REDFT01,2) # DCT-III (IDCT, dpsi/dx)
	u[:,:] = -r2r(Yshift.*Lc,FFTW.REDFT01,1) # DCT-III (IDCT, -dpsi/dy)
	plnt[:,:] = r2r(Xshift.*Kc,FFTW.REDFT01,2).*beta # DCT-III,  planetary vorticity 
	return u,v,plnt
end

function xcov(X,Y)
	C_XX = X*Y'; # covariance matrix
	N = length(X); Nd = 2*N-1; 
	rho = zeros(Nd); # unbiased autocorrelation vector
	lags = zeros(Nd); 
	stderr = zeros(Nd); # standard error formula for the auto-covariance
	for lag = 1:Nd
    		rho[lag] = mean(diag(C_XX,lag-N));
    		stderr[lag] = sqrt(var(diag(C_XX,lag-N))/length(diag(C_XX,lag-N)));
    		lags[lag] = lag-N;
	end
	return rho,lags,stderr
end


