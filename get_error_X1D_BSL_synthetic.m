function [e2,A] = get_error_X1D_BSL_synthetic(X,y,plant,Tmax)

Tmax_sim = (size(y,2))*plant.delt + .5; % always simulate 500 ms further than the data
% Tmax_sim = (size(y,2))*plant.delt; 

% perform simulations
sim = sim_vel_X1D_BSL_synthetic(X,plant,Tmax_sim); 

imax = ceil(Tmax/plant.delt); % max timestep to compare between model and data
if(imax>size(y,2))
    error('Tmax too long')
end

time = plant.delt;

% calculate error between simulation and acual data 
% e2 = nanmean((y(1:imax)-sim.acc(1:imax)).^2);
e2 = nanmean((y(1:imax)-sim.x(1:imax)).^2);

% A = y(1:imax)/sim.acc(1:imax); 
% e2 = nanmean((y(1:imax)-A*sim.acc(1:imax)).^2);


      
%     figure(3); clf; hold on  
%     plot(time*(1:length(y(1:imax))),y(1:imax),'m')
%     plot(time*(1:length(sim.acc(1:imax))),sim.acc(1:imax),'b')
%     plot(time*(1:length((y(1:imax)-sim.acc(1:imax)))),(y(1:imax)-sim.acc(1:imax)),'r')
% 
%     plot(y,'m')
%     plot(sim.acc,'b')
%     plot((y(1:imax)-sim.acc(1:imax)),'r')

sim.e2 = e2;

