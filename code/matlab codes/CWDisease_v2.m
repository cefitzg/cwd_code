% chronic wasting disease  ode model 
function  desolver
clear all 
global game gami K eps  mue mui  rhosW rhoiW  r rhofac gam xiW xifac tend Wmx

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 1.7);

% using new notation and new model with SI term
 formatSpecF = '%6.2f\n';
 % specify parameters
 d = 0.1;
 b = 0.6;
 r=b-d;
 
 mui = 0.6;
  
 mue=0.2;
 game = 1;  %wlog
 gami = 0.1;
  
 eps = 0.1 ;
  rhofac = 0.5;
 xifac=0; 
 
  K = 30;
   gam = eps*game+gami*mue;
 % the list of rho_iW values
plist = [0,0.1,0.15, 0.5 ];
 
 % figure(10)
 plot(K*ones(length(plist),1),plist,'*')
 hold on

for j = 1:length(plist)
    rhoiW = plist(j);   
    rhosW = rhoiW*rhofac;
    xiW = xifac*rhoiW;


 %specify the output points
t_start = -1500;
tstep = 0.1;
t_end = 100;
tend=t_end;
tspan = [t_start:tstep:t_end];

s0 = [K;0.01;0];  %  format s0 = [U;V;...] as a column vector
% integrate the ode's
[T,S] = ode23s(@deRHS,tspan, s0, odeset('maxstep',1));  
X = S(:,1);
Y = S(:,2);
E = S(:,3);
 % plot the time course of the solution  
figure(j)
plot(T,X/K,T,Y/K  , 'linewidth',2)
legend('boxoff')
legend('S/K','I/K','fontsize',20)
 xlabel('t (years)','fontsize',20)
  axis([0 tend 0 1])
  title( strcat('\rho_iW = ',sprintf(formatSpecF,rhoiW)),'fontsize',18)
  hold on

  % plot the I-S phase portrait of the solution
  figure(4+j)
plot( X/K, Y/K , 'linewidth',2)
%legend('S/K','I/K','E','fontsize',20)
 xlabel('S/K','fontsize',20)
 ylabel('I/K')
  title( strcat('\rho_iW = ',sprintf(formatSpecF,rhoiW)),'fontsize',18)
  
end
 
 %the bifurcation diagram
 
% the endemic curve
 
   rhoiW = 0:0.01:0.99;
   rhosW =  rhoiW*rhofac;
   xiW = xifac*rhoiW;
  
  Kt = r*mue*( mui  +  rhoiW)./( eps.*(1+xiW)*game.*(r - rhosW) + gami*mue*r -  gami*mue*rhosW);
  
  figure(10)
plot( Kt, rhoiW , Kt,  r*ones(length(K),1)/rhofac)
hold on
 
 
xlabel('K')
ylabel('\rho_iW')
text(60,0.1,'III','fontsize',18)
text(60,0.5,'II','fontsize',18)

text(10, 0.9 ,'I','fontsize',18)

text(10,1.2,'IV','fontsize',18)
   %the Hopf bifurcation curve
  
  [kp,rhoiW] = meshgrid(0:.25:150,0:.01:1);
   
  Hc=Hpf(kp,rhoiW);
 
 figure(10)
 contour(kp,rhoiW,Hc,[0,0] ,'linewidth',2)
 
  axis([0 100 0 1]) 
 
hold off
 
end

%the right hand side for ode simulation:
function s_prime=deRHS(t,s)
global game gami K eps  mue mui  rhosW rhoiW  r xiW
S = s(1);
II = s(2);
E = s(3);
 
fS= r*S*(1- (S+II)/K)-game*S*E -gami*S*II -rhosW*S;
                
fII=game*S*E+ gami*S*II -mui*II-rhoiW*II;
                  
fE=eps*(1+xiW)*II-mue*E;

 s_prime = [fS,fII, fE]'  ;
end

% the Hopf bifurcation curve
function out = Hpf(K,rhoiW)
global  game gami  eps  mue mui r rhofac xifac 
 rhosW=rhoiW*rhofac;
 xiW = xifac*rhoiW;
    gam = eps*(1+xiW)*game+gami*mue;
 M = rhoiW+mui;
   beta = (M.*(3*mue+r - rhosW) + M.^2 +  mue^2);
  % this formula assumes xiW=0;
   c3=eps*game.*gam.^2.*(r-rhosW).*(gam-M*gami);
   c2=-gam.*r.*(gami*mue.*(2*M.^2*gami*mue-3*gam.*M.^2-3*gam.*M*mue)+beta.*gam.^2);
   c1=-mue*r^2.*(gami*mue*(M.^2*gami*mue-3*gam.*M.^2-2*gam.*M*mue)+beta.*gam.^2);
   c0=-M.*r^3*mue^3.*(gam-M*gami);
   out = K.^3.*c3 +K.^2.*c2+K.*c1+c0;
 
end
 

