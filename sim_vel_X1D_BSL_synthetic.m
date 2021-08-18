function sim = sim_vel_X1D_BSL_synthetic(X,plant,Tmax)

delt = plant.delt;

Tjump = ceil(X(1)/delt);
uamp = X(2);
alpha = X(3);
beta = X(4);

% adjust high level parameters
start = 0; % starting position of the hand start=0;
nstep = ceil(Tmax/delt); % number of timesteps to simulate

% initialize the state vector, x
x = 0:nstep;
x(1,1) = start; % set starting hand position

% simulate velocity with parameters X
% input is acc or vel
y1 = [];
y2 = [];
for i = 1:nstep+1
    if x(i) < Tjump
        y1(1,i) = alpha*x(i)+beta;  %%%%%%%%%%%%%%%%% 05/14/2020
    else
        y2(1,i) = alpha*x(i)+beta + uamp*power(max(x(i)-Tjump,0),2);
    end
end
if isempty(y2) == 1
    Y = y1;
else
    Y = [y1 y2(Tjump+1:end)];
end


% Y = uamp*power(max(x-Tjump,0),3)*delt; % vel = qubic, acc = quad
% Y = uamp*power(max(x-Tjump,0),2)*delt;

% save data to sim
sim.x = Y;
sim.acc = diff(Y)/delt;
sim.T = Tjump;
sim.delt = delt;
sim.plant = plant;
sim.X = X;




