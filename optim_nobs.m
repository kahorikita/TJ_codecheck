% bootstrap for each subject
% Uni/Bi * Bigjump(3cm)/Small(1.5cm)
% First five are the responses for -3, -1.5, 0, 1.5, 3 cm jumps. six is no jumps
% Need 'getResponses.m'

clear
clc

tic

%%%%%%%%%%%%%%% change here
cond = 0; % 0:speed, 1:non speed constraint 
tflag = 1; % 0:uni, 1:bi
jflag = 1; % 0:small(1.5cm), 1:big(3cm)
searchmode = 0; % 0:fmincon, 1:fminsearch
initmumode = 1; % 0:org init mu, 1:precise

numofbs = 100; % number of bootstrap
%%%%%%%%%%%%%%%

if cond == 0 && tflag == 0 % speed & uni
    numofblk = 1;
    numofsubject = 8;
    load 'BimanualSkillData_speed.mat'
elseif cond == 0 && tflag == 1 % speed & bi
    numofblk = 5;
    numofsubject = 8;
    load 'BimanualSkillData_speed.mat'
elseif cond == 1 && tflag == 0 % not speed & uni
    numofblk = 2;
    numofsubject = 13;
    load 'BimanualSkillData.mat'
elseif cond == 1 && tflag == 1 % not speed & bi
    numofblk = 7;
    numofsubject = 13;
    load 'BimanualSkillData.mat'
end 

Hz = 130;
dt = 1/130;
delt = .001;
plant.delt = delt;

if initmumode == 0
    Mustock = [0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16]; 
else
    Mustock = [0.05 0.0525 0.0550 0.0575 0.06 0.0625 0.0650 0.0675 0.07 0.0725 0.0750 0.0775 0.08 0.0825 0.0850 0.0875 0.09 0.0925 0.0950 0.0975...
               0.1 0.1025 0.1050 0.1075 0.11 0.1125 0.1150 0.1175 0.12 0.1225 0.1250 0.1275 0.13 0.1325 0.1350 0.1375 0.14 0.1425 0.1450 0.1475...
               0.15 0.1525 0.1550 0.1575 0.16 0.1625 0.1650 0.1675];
end
% set serch area for mu
Astock = 0.0025; 
lb = [0 0 -10 -10]; % [mu A alpha beta]
ub = [0.5 50 10 10];
% lb = [0 0 -Inf -Inf];
% ub = [Inf Inf Inf Inf];

% fmincon
Aeq = [];
beq = [];

cols(:,:,1) = [ 0 210 255; 255 210 0; 0 0 0; 210 0 255]/256;
cols(:,:,2) = [ 0 155 255; 255 100 0; 0 0 0; 155 0 255]/256;
cols(:,:,3) = [ 0 100 255; 255 0 0; 0 0 0; 100 0 255]/256;

for subj = 1:numofsubject
    
    for blk = 1:numofblk
        
        %%%%%%%%%%%%%%% optimization with experiment data %%%%%%%%%%%%%%%%%
        % combine velocity for each Task & Jump
        if tflag == 0 && jflag == 1 % Uni & Big
            if cond == 1
                r = getResponses([-d{subj}.Uni{blk}{1}.CrX_post ; d{subj}.Uni{blk}{5}.CrX_post],.03,dt);
            else
                r = getResponses([-d{subj}.Uni{1}.CrX_post ; d{subj}.Uni{5}.CrX_post],.03,dt); % speed
            end
        elseif tflag == 0 && jflag == 0 % Uni & Small
            if cond == 1
                r = getResponses([-d{subj}.Uni{blk}{2}.CrX_post ; d{subj}.Uni{blk}{4}.CrX_post],.015,dt);
            else
                r = getResponses([-d{subj}.Uni{2}.CrX_post ; d{subj}.Uni{4}.CrX_post],.015,dt); % speed
            end
        elseif tflag == 1 && jflag == 1 % Bi & Big
            r = getResponses([-d{subj}.Bi{blk}{1}.CrX_post ; d{subj}.Bi{blk}{5}.CrX_post],.03,dt);
        elseif tflag == 1 && jflag == 0 % Bi & small
            r = getResponses([-d{subj}.Bi{blk}{2}.CrX_post ; d{subj}.Bi{blk}{4}.CrX_post],.015,dt);
        end

        
        y_130 = r.vel';
        y = resample(y_130,1/delt,Hz); 
        edata.y{subj,blk} = y;
        
        % calculate aplha(slope) and beta(intercept) from initial 100ms
        alpha = (y(100)-y(1))/(100-1);
        beta = y(1);
        
        % get velocity threshold for simulation duration
        [Vy,Vind] = max(y(1:800)); % % find peak on velocity in initial 800ms
        [Ty,Tind] = find(y(1:800) > Vy/2); % simulation duration: jumponset~1/2 peak vel
        
        if isempty(Tind) == 1 || Tind(1) < 10
            
            % if there is no peak or peak comes too early, just plot velocity and no optimization
%             figure(subj+50); hold on;
%             if tflag == 0 % Uni
%                 subplot(1,2,blk); hold on;
%             elseif tflag == 1 % Bi
%                 subplot(2,3,blk); hold on; % Bi
%             end
%             plot(delt*(1:length(y)),y,'color',cols(4,:,1),'linewidth',1.5) % time back to sec
%             xlabel('Time (s)','FontSize',10)
%             ylabel('Velocity (m/s)','FontSize',10)
        else
            
            % simulation duration
            Tmax = (Tind(1)-1)*delt;
            imax = ceil(Tmax/plant.delt); % max timestep to compare between model and data
            
            % try different initial value for mu
            for i = 1:length(Mustock)
                
                Muinit = Mustock(i);
                Ainit = Astock; % initial A is fixed but optimized
                Xinit = [Muinit Ainit alpha beta];
                
                % optimization
                f_targ = @(X) get_error_X1D_BSL_synthetic(X,y,plant,Tmax);
                
                if searchmode == 0
                    temp = fmincon(f_targ,Xinit,[],[],Aeq,beq,lb,ub);
                else
                    temp = fminsearchbnd(f_targ,Xinit,lb,ub);
                end

                sim = sim_vel_X1D_BSL_synthetic(temp,plant,Tmax);
                error = nanmean((y(1:imax)-sim.x(1:imax)).^2);
                Xopt(i,:) = [Xinit temp error];
            end
            
            % find the best mu
            [M,I] = min(Xopt(:,9)); % 9:error
            sim = sim_vel_X1D_BSL_synthetic(Xopt(I,5:8),plant,Tmax);
            edata.opty{subj,blk} = sim.x(:,1:imax);
            edata.optparams{subj,blk} = Xopt(I,:);
            
        end
        
%         figure(subj); hold on;
%         if tflag == 0
%             subplot(1,2,blk); hold on; % Uni
%         elseif tflag == 1
%             subplot(2,4,blk); hold on; % Bi
%         end
%         plot(delt*(1:length(y)),y,'color',cols(4,:,1),'linewidth',1.5)
%         plot(delt*(1:imax),sim.x(:,1:imax),'color',cols(1,:,1),'linewidth',1.5)
%         plot([Xopt(I,5),Xopt(I,5)],[-2,3],'k','linewidth',1.5);
%         xlabel('Time (s)','FontSize',10)
%         ylabel('Velocity (m/s)','FontSize',10)
%         xlim([0 0.5]);
%         ylim([-0.2 0.2])   

        figure(20); hold on;
        if tflag == 0
            subplot(1,2,blk); hold on; % Uni
        elseif tflag == 1
            subplot(2,4,blk); hold on; % Bi
        end
        plot(1,Xopt(I,5),'ko')
%         xlabel('Time (s)','FontSize',10)
%         ylabel('Velocity (m/s)','FontSize',10)
%         xlim([0 0.5]);
%         ylim([-0.2 0.2]) 
          
    end
        
end

toc

