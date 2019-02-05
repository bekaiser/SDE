# 1D stochastic integration
# Bryan Kaiser 
# 11/4/16

# using DataArrays
# Pkg.build("PyCall")
# using PyPlot
# using PyCall
# using Base.FFTW
include("./functions.jl")

# use the Richardson extrapolation to compute for infinite sampling!
# 3 sims with the same lambda for a larger and larger time period
# is there a link between R_TT(tau) and the Richardson extrapolation?

# =============================================================================
# set up parameters

dt = 1.0; # days, time step
Nt = 2^13; # number of time steps 
decay_time0 = [0.5,1.0,1.25,1.5]; # 1.75 2.0 2.25 2.5 2.75 3.0 3.5 4.0 4.5 5.0 6.0 10.0]; # days 
plot_flag = 0;
NE = 1000; # number of simulations to ensemble average
decay_time = 0.0; # initialization
lambda_computed = 0.0; # initialization
lambda_error = 0.0; # initialization
ratio_computed = 0.0; # initialization
A = 0.0; # initialization

for i = 1:length(decay_time0)

decay_time = decay_time0[i]; 
lambda_computed = 0.0; lambda_error = 0.0;
A = -(decay_time^(-1.0)); # 1/days, A = lambda, the linear autonomous 
# dynamical operator (feedback factor defined by Frankignoul & Hasselmann 1977)
println("Prescribed atmospheric time scale (1/q) = $dt days")
println("Prescribed ocean time scale (1/lambda) = $decay_time days")
println("Prescribed lambda = $A 1/days")

tic()
for j = 1:NE
#println("trial $j");
#println("$(toq()) s");
#tic()


# =============================================================================
# The Stratonovich stochastic differential equation is advanced using the RK4
# method shown in Hansen & Penland (2006).

#Nlag = 30;
Nlag = Int32(round(decay_time+0.01))+1; # size of the lag window to average over
# the minimum decay_time is 0.5 days
N = Nt;
T = zeros(N); # C, (modeled) sea surface temperature anomaly
q_dt = zeros(N); # forcing
t = zeros(N); # time
t[1] = dt/2.0; # days, starting time
for n = 2:N
     	# Gaussian white noise forcing:
    	q_dt[n-1] = randn()*sqrt(dt);
	# magnitude?
    	# a random number from a normal distribution with variance dt and 
	# mean 0.
    
 	# advance in time:
    	k1 = sde_rk4(T[n-1],q_dt[n-1],A);
    	k2 = sde_rk4(T[n-1]+k1*dt/2.0,q_dt[n-1],A);
    	k3 = sde_rk4(T[n-1]+k2*dt/2.0,q_dt[n-1],A);
    	k4 = sde_rk4(T[n-1]+k3*dt,q_dt[n-1],A);
    	T[n] = T[n-1]+(k1+k2.*2+k3.*2+k4).*(dt/6.0);
    	t[n] = t[n-1]+dt;

end

# center the data:
q_dt = q_dt[2:N-1]; T = T[2:N-1]; # collocate in time
L = t[N-1]-t[2]+dt; # total time elapsed
t = t[2:N-1]; 


# =============================================================================
# compute Q

# wavenumbers  
N = length(T); k = Array{Float64}(N); k[1] = 0.0;
k[2:Int32(N/2)+1] = collect(1:Int32(N/2)).*(2.0*pi/L) # rad/km
k[Int32(N/2)+2:N] = -collect(Int32(N/2)-1:-1:1).*(2.0*pi/L) # rad/km

TFFT = fft(T);
kabs = abs(k);
for m=1:length(k)
    if kabs[m] >= pi/2.0
        TFFT[m] = 0+0*im;
    end
end

# derivative
Q = real(ifft(TFFT.*k.*im)); # K/day

# throw away end points
Q = Q[100:N-100]; T = T[100:N-100]; q_dt = q_dt[100:N-100]; t = t[100:N-100]; # (redefined)
N = length(T); # redefine
L = t[N]-t[1]+dt; # total time elapsed (redefined)
#println("total time = $L days")

# signal plot
if plot_flag == 1; fig = figure()
subplot(311); plot(t,T,"r"); xlabel("t"); ylabel("T"); 
subplot(312); plot(t,Q,"b"); xlabel("t"); ylabel("Q"); 
subplot(313); plot(t,q_dt,"g"); xlabel("t"); ylabel("q"); 
show(); end;


# =============================================================================
# covariances

RTT,lagsTT,stderr = xcov(T,T);
topbarTT = RTT+stderr;
bottombarTT = RTT-stderr;
#lags = lags.*dt;
RTQ,lagsTQ,stderrTQ = xcov(T,Q);
RQQ,lagsQQ,stderrQQ = xcov(Q,Q);

# covariance plot
if plot_flag == 1; fig = figure()
plot(lagsTT,RTT,"b",label="RTT"); 
plot(lagsTT,topbarTT,"g--",label="error"); 
plot(lagsTT,bottombarTT,"g--");
xlabel("lag (days)"); legend(); show(); end;

# find that lag of zero crossing (right side) and integrate autocovariance:
# integration on the left side (lag<0)
locT0 = find(RTT -> RTT <= 0.0,RTT[1:N-1])[end];
NT0 = Float32(length(lagsTT[locT0:N]));
hT = (lagsTT[N+1]-lagsTT[locT0])/NT0;
tauo = (hT/2.0*(RTT[locT0]+(sum(RTT[locT0+1:N-1]))*2.0+RTT[N]))/RTT[N];

# find that lag of zero crossing (right side) and integrate autocovariance:
# integration on the left side (lag<0)
locQ0 = find(RQQ -> RQQ <= 0.0,RQQ[1:N-1])[end];
NQ0 = Float32(length(lagsQQ[locQ0:N]));
hQ = (lagsQQ[N+1]-lagsQQ[locQ0])/NQ0;
taua = (hQ/2.0*(RQQ[locQ0]+(sum(RQQ[locQ0+1:N-1]))*2.0+RQQ[N]))/RQQ[N];

ratio_computed = (tauo/taua+(Float32(j)-1.0)*ratio_computed)/Float32(j);


# =============================================================================
# feedback factor estimate

lambda = -RTQ./RTT;

if plot_flag == 1; fig = figure()
plot(lags,lambda,"k"); 
xlabel("lag (days)"); axis([-2.0*decay_time,2.0*decay_time,-3.0*A,3.0*A]); 
show(); end

# lambda error (running average)
lambda_computed = (mean(lambda[N-Nlag:N-1])+(Float32(j)-1.0)*lambda_computed)/Float32(j);
lambda_error = abs(A-lambda_computed);

end # ensemble simulations

println("Time to compute $NE ensemble averaged realizations: $(round(toq()))s");
println("Time scale ratio estimate = $(round(ratio_computed,3))");
println("Lambda estimate = $(round(lambda_computed,6)) 1/days");
println("Lambda estimate error = $(round(lambda_error,6)) 1/days");
println("Lambda estimate percent error $(round(lambda_error/A*100.0,3))%\n");

end # for the given decay time
