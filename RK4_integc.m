function [tout, yout] = RK4_integcq(m,c,k_2,k_3)
global sig 
% The Runge Kutta method is used for integration
%--------------------------------------------------------------------------
h = 0.01;   % step size
t0=0;tf=5;  % Initial and final times
y0=[0 0]'; % Initial conditions
n= ceil((tf-t0)/h); % number of steps
%--------------------------------------------------------------------------
t=t0; y=y0; % Initialize t and y
output = [t0 y0']; % first row of matrix of printed values
w = [t0,y0'];
%
for i=1:n
    k1=h*c_model(t,y,m,c,k_2,k_3,sig(i));
    k2=h*c_model(t+h/2,y+k1/2,m,c,k_2,k_3,sig(i));
    k3=h*c_model(t+h/2,y+k2/2,m,c,k_2,k_3,sig(i));
    k4=h*c_model(t+h,y+k3,m,c,k_2,k_3,sig(i));
    z = y + (1/6)*(k1+2*k2+2*k3+k4);
    t = t + h;    
    y = z;    
    w = [w;t z'];    
end
tout = w(:,1); % time
yout = w(:,2); % displacement response