%make figures for CWD from xpp data files
clear all 
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
% specify parameters:
global r mui gam eps K mue rhorat
  
 r = 0.5;
 mui = 0.6;
 mue = 0.25;
 game = 0.8;
 gami = 0.1;
  
 eps = 0.1 ;
  rhorat = 0.5;
  gam = eps*game+gami*mue;
 xifac=0; 
  K=75;

formatSpecF = '%6.2f\n';

load('diagram5.dat')  % this comes from xpp
A=diagram5;
load('diagram6.dat')
B = diagram6;

indxc=find(A(:,4)==3);
%
% first things first:  Plot the periodic solutions
[a,j] = min(A(indxc,1));
  %   find   the hopf  bifurcation point:
  [hpt,hj] =max(A(indxc,1))
   [S,I]=ssend(hpt)
figure(1)
 plot(A(indxc(1:j-1),1),A(indxc(1:j-1),2),'g*' , ...
     A(indxc(1:j-1),1),A(indxc(1:j-1),3),'g*',  hpt,S,'k*') 
hold on
 figure(2)
  plot(B(indxc(1:j-1),1),B(indxc(1:j-1),2),'g*', ...
      B(indxc(1:j-1),1),B(indxc(1:j-1),3),'g*',  hpt,I,'k*')
hold on

rhoWpt= (gam*K*r - mui*mue*r)/(mue*r + gam*K*rhorat); 
rhoiWmx =r/rhorat;
figure(1)
[S,II] = ssend([0,hpt])
plot([0,hpt],S,'k--')
figure(2)
plot([0,hpt],II,'k--')
figure(1)
 
[S,II] = ssend([hpt,rhoWpt])
plot([hpt,rhoWpt],S,'r','linewidth',2)
figure(2)
plot([hpt,rhoWpt],II,'r','linewidth',2)

 
% now the disease free solution
 
 
[S,I] = ss0([0,rhoWpt]);
figure(1)
plot([0,rhoWpt],S,'k--','linewidth',2)
 [S,I] = ss0([ rhoWpt,rhoiWmx]);
 plot([rhoWpt,rhoiWmx],S,'r','linewidth',2)
 figure(1)
[S,I] = ss0([rhoWpt,rhoiWmx])
plot([rhoWpt,rhoiWmx],S,'r','linewidth',2)
plot([rhoiWmx,1.2*rhoiWmx],[0,0],'r','linewidth',2)
figure(2)
%plot([rhoWpt,rhoiWmx],I,'r','linewidth',2)
plot([rhoWpt,1.2*rhoiWmx],[0,0],'r','linewidth',2)
 plot([0,rhoWpt],[0,0],'k--')
% now finish off the plots
figure(1)
axis([0 1.2*rhoiWmx 0 15])
xlabel('\rho_iW')
ylabel('S')
plot([hpt,hpt],[0,K],'k')
plot([rhoWpt,rhoWpt],[0,K],'k')
plot([rhoiWmx,rhoiWmx],[0,K],'k')
 
text(0.03,12,'III','fontsize',18)
text(0.55,12,'II','fontsize',18)
text(0.97,12,'I','fontsize',18)
text(1.1,12,'IV','fontsize',18)
 %title( strcat('K = ',sprintf(formatSpecF,K)),'fontsize',18)
hold off

figure(2)
axis([0 1.2*rhoiWmx 0 15])
xlabel('\rho_iW')
ylabel('I')
plot([hpt,hpt],[0,K],'k')
plot([rhoWpt,rhoWpt],[0,K],'k')
text(0.03,12,'III','fontsize',18)
text(0.55,12,'II','fontsize',18)
text(0.97,12,'I','fontsize',18)
text(1.1,12,'IV','fontsize',18)
plot([rhoiWmx,rhoiWmx],[0,K],'k')
text(1.1,75,'IV','fontsize',18)
% title( strcat('K = ',sprintf(formatSpecF,K)),'fontsize',18)
hold off
  
function [S,I]=ss0(rhoiW)
% the disease free state
global r K rhorat

rhosW = rhorat*rhoiW;

S = K*(1-rhosW/r);
I = zeros(length(rhoiW),1);
end

function [S,I]=ssend(rhoiW)
% the disease endemic state
global r K gam rhorat mue mui

 rhosW = rhorat*rhoiW;

S =  mue/gam*(mui+rhoiW);
 
I =  mue*(r*(K*gam-mui*mue)-(K*gam*rhosW+r*mue*rhoiW))./((K*gam+mue*r)*gam);

end


   