function dydt = c_model(t,y,m,c,k_2,k_3,sr)
global sig 
f = sr; % Random excitation
dydt = [y(2);-(k_2/m)*y(1)-(c/m)*y(2)-(k_3/m)*y(1).^3+(f/m)];
end

