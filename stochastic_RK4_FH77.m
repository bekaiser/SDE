% stochastic integration (4th-order Runge-Kutta) example
% Bryan Kaiser
% 1/21/16

% This script generates a sea surface temperature anomaly signal by stochastically integrating 
% the Frankignoul & Hasselmann (1977) stochastic forcing model of ocean-atmosphere heat flux
% interaction. The constructed signal is stochastically integrated by the Runge Kutta 4th-order 
% method for Stratonovich stochastic differential equations as shown by Hansen & Penland (2006).

% This script must be run in a folder containing the files "mcolormaps.m" and "sde_RK4.m" 

close all; clear all; clc;
run('mcolormaps.m');
movieflag = 1;

%==========================================================
% signal

dt = 1; % days, time step
Nt = 2^6; % number of time steps (i.e. observations) 90 years if each time step = 1 day
Nx = 10; % square matrix (i.e. square horizontal domain) side dimension
N = Nx^2; % number of state variables (i.e. grid points)
decay_time = 15; % days 

% A, the autonomous linear dynamical operator 
lambda = -ones(1,N)./decay_time; % 1/days
A = diag(lambda); % 1/days, linear, autonomous dynamical operator

% nonnormal if non-zero:
nonnormality = sum(sum(A*(A')-(A')*A))
% numerical abscissa
numerical_abscissa = max(max(eig((A+A')/2)))
% spectral abscissa
spectral_abscissa = max(max(real(eig(A))))

% spatial correlation matrix:
eta = 1; 
G = eye(N).*eta; % no spatial correlation


%==========================================================
% The Stratonovich stochastic differential equation is advanced using the RK4
% method shown in Hansen & Penland (2006).

T = nan(N,Nt); % (modeled) SSTa for all grid points and times
Q = nan(N,Nt);
T0 = randn(N,1); % initial SST 
T(:,1) = T0; % initial SST 
t = 0; % starting time
for n = 2:Nt
    
    z_dt = randn(N,1).*sqrt(dt);
    % z_dt = random number from a normal distribution 
    % with variance dt and mean 0. (Brownian)
    
    k1 = sde_RK4(T(:,n-1),z_dt,A,G);
    k2 = sde_RK4(T(:,n-1)+k1*dt/2,z_dt,A,G);
    k3 = sde_RK4(T(:,n-1)+k2*dt/2,z_dt,A,G);
    k4 = sde_RK4(T(:,n-1)+k3*dt,z_dt,A,G);
    
    % advance in time:
    T(:,n) = T(:,n-1)+(k1+k2.*2+k3.*2+k4).*(dt/6); % T'
    t = t+dt;
end

absolute_maximum_T = max(max(abs(T)))
total_sampled_time = dt*Nt; % days

%==============================================================================
% results movie

if movieflag == 1 % show the movie ....

Ny = Nx;
Lx = 100; % km, domain size
Ly = Lx; % km
Lxcenter = 0.0; % x value @ the center of the grid
Lycenter = 0.0; % y value @ the center of the grid
dx = Lx/Nx; % km, uniform longitudinal grid spacing
dy = Ly/Ny; % km, uniform longitudinal grid spacing
x = [0.5*dx:dx:dx*Nx]-(Lx/2.0-Lxcenter); % km, centered uniform grid 
y = [0.5*dy:dy:dy*Ny]-(Ly/2.0-Lxcenter); % km, centered uniform grid
[X,Y] = meshgrid(x,y);

% reshape T into the 2D spatial grid:
Tmovie = reshape(T,Nx,Nx,Nt);
Nmovie = Nt; % unitless number of time steps  

count = 1; % initial movie counter for watching
figure
for n = 2:Nmovie
	contourf(X,Y,squeeze(Tmovie(:,:,n)),100); shading flat
	xlabel('X'); ylabel('Y'); %1title('nabla2(psi) error'); 
    title(['t=' num2str(floor(n*dt)) ' days']);
	colorbar; colormap(dawn/255);  colormap(flipud(colormap)); 
    %caxis([-0.05,0.05]);

    % watch:
    M(count)=getframe;
 	count = count+1;

end

end % movieflag == 1
