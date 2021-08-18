function response = getResponses(x,pmax,dt)
% returns mean responses as:
%   response.pos
%   response.vel
%   response.acc
%   response.time

%x(isnan(x)) = pmax; % eliminate nans

pos = x;
pos(isnan(pos)) = pmax; % replace nans in position with target location
response.pos = mean(pos)';

vel = diff(x')'/dt;
vel(isnan(vel)) = 0; % replace nans in velocity with 0s

vel = savgolayFilt(vel,3,7);

response.vel = mean(vel)'; % compute smoothed, averaged velocity

response.acc = savgolayFilt(diff(response.vel)'/dt,3,7)';

%dAll.posResponse_large(subj,:,c) = nanmean([-d{subj}.Bi{c}{1}.CrX_post ; d{subj}.Bi{c}{5}.CrX_post]);
response.time = 1000*dt*([1:size(x,2)]-1)'-100; % ms not seconds, and subtract measurement 100 ms delay
%response.vel = diff(savgolayFilt(response.pos',3,7))'/dt;

%response.vel(isnan(response.vel)) = 0;

%response.acc = diff(savgolayFilt(response.vel',3,7))'/dt;

response.vel_var = std(vel)';