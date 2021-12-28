function y = c_simulate(m,c,k_2,k_3)
global sig 
% Function SIMULATE. Simulates response of a  SDOF system
% y[]            : Matrix of displacement and velocity response values.
%
% sig            : Array of force values.
% m              : Mass of SDOF system (equal to 1, supposed to be known).
% c              : Damping 
% k_2, k_3       : Stiffness coefficients.
% Integrate forward using RK4. 
[t,y] = RK4_integc(m,c,k_2,k_3);