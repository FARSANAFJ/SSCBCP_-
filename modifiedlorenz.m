clear all; close all;
beta=[10;28;8/3;-3.6;5.2];
x0 = [0; 1; 20; 20]; tspan = 0.001:0.001:50;
 %options=odeset()
[t,x] = ode45(@(t,x)mod_lorenz_ode(t,x,beta),tspan,x0);
%%%%%%%key generatiom
k1=x(:,1); k2=x(:,2);k3=x(:,3);k4=x(:,4);
key_sub1=bitxor(bitxor(int64(k1),int64(k2)),bitxor(int64(k3),int64(k4)));
save('key_sub1','key_')