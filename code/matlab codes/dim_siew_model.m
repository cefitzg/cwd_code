clear
clc
close all

%Use Keener's plotting options 
 formatSpecF = '%6.2f\n';
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 1.7);

%dimensional system

% specify parameters
d = 0.1;
b = 0.6;
r=b-d;
mui = 0.6;
mue=0.2;
game = 1;  %wlog
gami = 0.1;
eps = 0.1 ;
K = 30;
rhoiW = 0.01; 
omega=30; %bifurcation parameter
rho = rhoiW; 
alpha=0.05; 
rw=1; 
Kw=3.0;

%specify the output points
t_start = -15000;
tstep = 0.1;
t_end = 1000;
tend=t_end;
tspan = [t_start:tstep:t_end]';

f=@(t,P)[r*P(1).*(1- (P(1)+P(2))/K) - game*P(1).*P(3) - gami*P(1).*P(2) - rho*P(1).*P(4); ...
         game*P(1).*P(3) + gami*P(1).*P(2) - mui*P(2) - omega*rho*P(2).*P(4); ...
         eps*P(2) - mue*P(3); ...
         rw*P(4).*(1-P(4)/Kw) + alpha*(rho*P(1).*P(4) + omega*rho*P(2).*P(4))];


s0 = [K;0.01;0;3];  %  format s0 = [U;V;...] as a column vector
% integrate the ode's
[T,S] = ode23s(f,tspan, s0);
X = S(:,1);
Y = S(:,2);
E = S(:,3);
W = S(:,4); 


% plot the time course of the solution
figure(1)
plot(T,X/K, 'linewidth',2)
hold on 
plot(T,Y/K, 'linewidth',2)
plot(T,E, 'linewidth',2)
%plot(T,W/K2, 'linewidth',2,'Color','#702963')
%legend('boxoff')
legend('S/K','I/K','E','fontsize',20,'Interpreter', 'latex')
xlabel('t (years)','fontsize',20)
axis([0 200 0 1.5])
title( strcat('\omega = ',sprintf(formatSpecF,omega)),'fontsize',18)
hold on
% 

% plot the time course of the solution
figure(2)
%plot(T,X/K, 'linewidth',2)
%hold on 
%plot(T,Y/K, 'linewidth',2)
%plot(T,E, 'linewidth',2)
plot(T,W/Kw, 'linewidth',2,'Color','#702963')
%legend('boxoff')
legend('W/K_{w}','fontsize',20)
xlabel('t (years)','fontsize',20)
axis([0 200 1 1.05])
title( strcat('\omega = ',sprintf(formatSpecF,omega)),'fontsize',18)
hold on
% 