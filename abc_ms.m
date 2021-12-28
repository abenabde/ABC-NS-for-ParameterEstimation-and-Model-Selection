function theta = abc_ms(um_1,sig,tol,accuracy)
tt = cputime;
n_sim=100; % Number of simulations
nparam_1 = 2; % number of unknown parameters in M_1
nparam_2 = 3; % number of unknown parameters in M_2
threshold_vec(1,1) = tol; % Initial value of the tolerance threshold 
pop_size = 1000; % Population size
f_0 = 1.1; % Enlargement factor
beta_0 = 0.6; % proportion of remaining 'alive' particles
alpha_0 = 0.3; % proportion of dropped particles
Av = 100/(length(um_1)*std(um_1)^2); % normalisation constant
population_1 = zeros(pop_size,n_sim*nparam_1); % matrix to store the populations for M_1
nmse_1 = ones(pop_size,n_sim); % matrix to store the NMSE values for model 1
nmse_2 = ones(pop_size,n_sim);% matric to store the NMSE values for model 1
population_2 = zeros(pop_size,n_sim*nparam_2); % matrix to store the populations for M_2
ix = zeros(n_sim,2);
ia = zeros(n_sim,2);
ic = zeros(n_sim,2);
i_new = zeros(n_sim,1);
acceptance = zeros(n_sim,1); % store the acceptance rate
adap_prior = zeros(n_sim,2); % 
adap_prior(1,1:2) = [0.5 .5]; % Equal prior probabilities
weight_1 = ones(pop_size,n_sim);
weight_2 = ones(pop_size,n_sim);
count1 = 0;count2 = 0;iter = 0;
sim_0 = 0;rej=0;
sim=1;
fprintf('~~~~~~~~~~~~~~~~~~~~ABC-NS algorithm: Simulation has started~~~~~~~~~~~~~~~~~~~~\n ')
%============================================================
%
%%       STEP 1 : First loop / GENERATE THE FIRST POPULATION
% 
%============================================================
while iter < pop_size
       sim_0 = sim_0 + 1;
       ap=find(1==mnrnd(1,[adap_prior(1,:)]'));
       if ap == 1
           % Generate particles from uniform priors
          theta_2s = unifrnd([0.002 20],[0.2 80]);
          % Simulate from the linear model
          uv = c_simulate(1,theta_2s(1,1),theta_2s(1,2),0);
          NMSE = Av*sum((um_1 - uv).^2); % evaluate the NMSE
           if NMSE <= threshold_vec(1,1)
              pop_1(count1+1,:) = theta_2s;
              zk1(count1+1,:) =  NMSE;
              count1 = count1 +1;
           end
       else
            % Uniform prios on the unknown model parameters
            theta_2s = unifrnd([0.002 20 500],[0.2 80 1500]);
           % Simulate from the cubic model     
           uv = c_simulate(1, theta_2s(1,1), theta_2s(1,2), theta_2s(1,3)); % simulate from the cubic model
           NMSE = Av*sum((um_1 - uv).^2); % evaluate the NMSE
           if NMSE <= threshold_vec(1,1)
              pop_2(count2+1,:) = theta_2s;
              zk2(count2+1,:) =  NMSE;
              count2 = count2 +1;
           end     
       end       
        iter = count2 + count1;
end
%==========================================================================
% store population 1
population_1(1:count1,1:nparam_1) = pop_1;
population_2(1:count2,1:nparam_2) = pop_2;
nmse_1(1:count1,1) = zk1;
nmse_2(1:count2,1) = zk2;
acceptance(1,1) = pop_size/sim_0;
ix(1,1:2) = [count1 count2];
[svc,ind1] = sort([zk1;zk2],'descend');
threshold_vec(1,2) = svc(alpha_0*pop_size); % define the next tolerance threshold
% Assign a weight for each particle in model 1
for i=1:ix(1,1)
       if zk1(i,1)>=threshold_vec(1,2)
           w1(i,:)=0;
       else 
           w1(i,:)=(1/threshold_vec(1,1))*(1 - (zk1(i,1)./threshold_vec(1,1)).^2);
       end
end
% Assign a weight for each particle in model 2
for i=1:ix(1,2)
       if zk2(i,1)>=threshold_vec(1,2)
           w2(i,:)=0;
       else 
           w2(i,:)=(1/threshold_vec(1,1))*(1 - (zk2(i,1)./threshold_vec(1,1)).^2);
       end
end
%--------------------------------------------------------------------------
%% Select the active particles that they will remain for both models
%--------------------------------------------------------------------------
n_ac = round(ix(1,1)*beta_0);
[seq_1,idx1] = datasample([pop_1],round(ix(1,1)*beta_0),'Weights',[w1'],'Replace',false);
active_pop_1 = seq_1;
obj_pop_1 = zk1(idx1);
[seq_2,idx2] = datasample([pop_2],round(ix(1,2)*beta_0),'Weights',[w2'],'Replace',false);
active_pop_2 = seq_2;
obj_pop_2 = zk2(idx2);
active_1(1:round(ix(1,1)*beta_0),1:nparam_1) = active_pop_1;
active_obj1(1:round(ix(1,1)*beta_0),1) = obj_pop_1;  % value of the distance in the first population
active_2(1:round(ix(1,2)*beta_0),1:nparam_2) = active_pop_2;
active_obj2(1:round(ix(1,2)*beta_0),1) = obj_pop_2;  % value of the distance in the first population 
ia(1,1:2) = [round(ix(1,1)*beta_0) round(ix(1,2)*beta_0)];
ic(1,1:2)=[round(ix(1,1)*0.4) round(ix(1,2)*0.4)];
i_new(1,1) = sum(ic(1,:));
%--------------------------------------------------------------------------
%% Compute the properties of the best ellipsoid
mu_1 = mean(active_pop_1); % mass center
B_1 = cov(active_pop_1); % covariance matrix
D_1=nparam_1; % Number of the parameters in model 1
const_1 = pi^(D_1/2)/gamma(D_1/2 + 1);
VS_1 = const_1*sqrt(det(B_1));
[Bs, mus_1, VEs, ns] = calc_ellipsoid(active_pop_1, VS_1);
Bsn_1 = Bs*f_0; % Enlarge the ellipse
mu_2 = mean(active_pop_2); % mass center
B_2 = cov(active_pop_2);   % covariance matrix
D_2=nparam_2; % Number of parameters in model 2
const_2 = pi^(D_2/2)/gamma(D_2/2 + 1);
VS_2 = const_2*sqrt(det(B_2));
[Bs, mus_2, VEs, ns] = calc_ellipsoid(active_pop_2, VS_2);
Bsn_2 = Bs*f_0; % Enlarge the ellipse
fprintf('\nABC-NS population #%d -- Tolerance threshold %4.2f',sim,threshold_vec(1,sim))
%============================================================
%
%%         STEP 2 : Seconf loop / GENERATE THE NEXT POPULATIONS
% 
%============================================================
for sim=2:100    
    %------------------------------------------
    count1 = 0; count2 = 0; iter = 0;  
         z1 = []; z2 = []; 
         w1 = []; wn1 = []; 
         w2 = []; wn2 = []; 
         comp_pop_1 = [];comp_pop_2=[];
         pop_1=[];pop_2=[];
         obj1=[];obj2=[];
         active_pop_1 = []; active_pop_2 = [];
         obj_pop_1 = []; obj_pop_2 = [];
         svc=[];svc1=[];
         ww1 = [];ww2=[];seq1=[];idx1=[];
         seq2=[];idx2=[];
    %----------------------------------------
       if ia(sim-1,2) == 0
             adap_prior = [1 0];
           elseif  ia(sim-1,1) == 0
               adap_prior = [0 1];
         else
             adap_prior = [0.5 0.5];
       end
   %********************************************************************
   adap_prior ;
   rej = 0;
     while iter < 400  
           m=find(1==mnrnd(1,[adap_prior(1,:)]'));
           rej = rej +1 ;
           if m == 1
             
%                disp('Model 1 is selected')
                sm=0;
             while sm<1
                 % Sample inside the ellipse
             theta_new = draw_from_ellipsoid(Bsn_1, mus_1, 1);
               %----------------Uniform Kernel ----------------------------
               if (theta_new(1,1)>=0.002 && theta_new(1,1)<=0.2 && theta_new(1,2)>=20 &&  theta_new(1,2)<=80) 
                   sm = sm + 1;
                   break
               end       
             end   
               % Simulate from the linear model
               uv = c_simulate(1,theta_new(1,1),theta_new(1,2),0);
               % Evaluate the NMSE 
              NMSE = Av*sum((um_1 - uv).^2);
             
               if  NMSE <= threshold_vec(1,sim)
                   comp_pop_1(count1+1,:) = theta_new;
                   z1(count1+1,:) =  NMSE;
                  count1 = count1 +1;
               end      
           else    
%                disp('Model 2 is selected');
                sm=0;
             while sm<1
                 % Sample inside the ellipse
             theta_new = draw_from_ellipsoid(Bsn_2, mus_2, 1);            
               %----------------Checking step----------------------------
               if (theta_new(1,1)>=0.002 && theta_new(1,1)<=0.2 && theta_new(1,2)>=20 &&  theta_new(1,2)<=80 && ...
                   theta_new(1,3)>=500 &&  theta_new(1,3)<=1500) 
                   sm = sm + 1;
                   break
               end       
             end            
%             Simulate from the cubic model Model 2
              uv = c_simulate(1,theta_new(1,1),theta_new(1,2),theta_new(1,3));
              % Evaluate the NMSE 
              NMSE = Av*sum((um_1 - uv).^2);
             
               if  NMSE <= threshold_vec(1,sim)
                   comp_pop_2(count2+1,:) = theta_new;
                   z2(count2+1,:) =  NMSE;
                   count2 = count2 +1 ;
               end                     
           end
 
           iter = count1 + count2;
       
     end 
     % complete the population for M_1
     pop_1 = [active_1(1:ia(sim-1,1),2*sim-3:2*sim-2);comp_pop_1];
     obj1 =[active_obj1(1:ia(sim-1,1),sim-1);z1];
     % complete the population for M_2
     pop_2 = [active_2(1:ia(sim-1,2),3*sim-5:3*sim-3);comp_pop_2];
     obj2 =[active_obj2(1:ia(sim-1,2),sim-1);z2];         
      [svc,ind1] = sort([obj1;obj2],'descend');
      % Define the next tolerance threshold
      threshold_vec(1,sim+1) = svc(alpha_0*pop_size)  ; 
     % Update the weights of the particles for Model 1
            for vv=1:length(obj1)
                if obj1(vv,1)>=threshold_vec(1,sim)
                   ww1(vv,:)=0;
                else 
                   ww1(vv,:)=(1/threshold_vec(1,sim))*(1 - (obj1(vv,1)./threshold_vec(1,sim)).^2);
                end
            end
             % Update the weights of the particles for Model 2
             for vv=1:length(obj2)
                if obj2(vv,1)>=threshold_vec(1,sim)
                   ww2(vv,:)=0;
                else 
                   ww2(vv,:)=(1/threshold_vec(1,sim))*(1 - (obj2(vv,1)./threshold_vec(1,sim)).^2);
                end
             end
                   ix(sim,:) = [size(pop_1,1) size(pop_2,1)];
                   ia(sim,:) = [round(size(pop_1,1)*beta_0) round(size(pop_2,1)*beta_0)];
                   ic(sim,:)=[length(z1) length(z2)];
                   i_new(sim,1) = sum(ic(sim,:));
     % Select the active particles     
               sd = length(find(ww1>0));
               sert = round(ix(sim,1)*beta_0);
               var1 = sert-sd;
               % Checking step
                 if sert<=sd 
                     KM = sert;
                 else 
                     KM = sd ;
                 end
               
                 if (KM == 0)
                     active_pop_1 = [];
               obj_pop_1 = [];
                 else
               [seq1,idx1] = datasample([pop_1],KM,'Weights',[ww1],'Replace',false);
               active_pop_1 = seq1;
               obj_pop_1 = obj1(idx1);              
                 end
               % checking step
               sert2 = round(ix(sim,2)*beta_0);
               if var1>0
                   sert2 = sert2 + var1 ;
               end
               
                [seq2,idx2] = datasample([pop_2],sert2,'Weights',[ww2],'Replace',false);
               active_pop_2 = seq2;
               obj_pop_2 = obj2(idx2);
                 % checking step
                 if length(active_pop_1)<ia(sim,1)
                     ia(sim,1)=ia(sim,1)- var1;
                     ia(sim,2) = ia(sim,2) + var1;
                 end
                                          
              active_1(1:ia(sim,1),2*sim-1:2*sim) = active_pop_1;
              active_obj1(1:ia(sim,1),sim) = obj_pop_1;
     
              active_2(1:ia(sim,2),3*sim-2:3*sim) = active_pop_2;
              active_obj2(1:ia(sim,2),sim) = obj_pop_2;
     
     
      %----Update Bs and mus    
                  mu_1 = mean(active_pop_1);
                  B_1 = cov(active_pop_1);
                  D_1=nparam_1;
                  
                  const_1 = pi^(D_1/2)/gamma(D_1/2 + 1);
                  VS_1 = const_1*sqrt(det(B_1));
                  [Bs, mus_1, VEs, ns] = calc_ellipsoid(active_pop_1, VS_1);
                  Bsn_1 = Bs*f_0;
                  %
                  mu_2 = mean(active_pop_2);
                  B_2 = cov(active_pop_2);
                  D_2=nparam_2;

                  const_2 = pi^(D_2/2)/gamma(D_2/2 + 1);
                  VS_2 = const_2*sqrt(det(B_2));
                  [Bs, mus_2, VEs, ns] = calc_ellipsoid(active_pop_2, VS_2);
                  Bsn_2 = Bs*f_0;
           
       population_1(1:ix(sim,1),2*sim-1:2*sim) = pop_1;
       nmse_1(1:ix(sim,1),sim) = obj1;
       
       population_2(1:ix(sim,2),3*sim-2:3*sim) = pop_2;
       nmse_2(1:ix(sim,2),sim) = obj2;               
       acceptance(sim,1) = 400/rej;
       op = ix(1:sim,:);     
      [svc1,ind2] = sort([z1;z2;obj_pop_1;obj_pop_2]);    
     % Stopping condition
     if threshold_vec(1,sim) - threshold_vec(1,sim+1)<accuracy
          break
      end
   
   fprintf('\nABC-NS population #%d -- Tolerance threshold %4.2f',sim,threshold_vec(1,sim))   
   
end
fprintf('\n~~~~~~~~~~~~~~~~~~~~ABC-NS algorithm has converged ~~~~~~~~~~~~~~~~~~~~\n ')
cpu = (cputime - tt)/60;
fprintf('\nCPU time : %2.4f minutes', cpu) 
theta=mean(pop_2);
% Desactivate for postprocessing
% fname = sprintf('Output_ABCNS_%d.mat',sim);
% save(fname)