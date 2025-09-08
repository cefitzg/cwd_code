% Gillespie simulation for the SIR process for CWD (Chronic Wasting
% disease)

% in this file we add a delayed control

% the reactions are:
% reaction 1 birth  s-> s+1 rate r
% reaction 2 death s -> s-1 rate d
% reaction 3 infection S + E -> I  rate gamma_e
% reaction 4 infection S + I ->2I  rate gamma_i
% reaction 5   death i -> i-1  rate mui
% reaction 6 deer predation S+W-> W rate rhos 
% reaction 7  infected predation I+W-> W rate rhoi
 
%meanwhile E changes gradually dE/dt = eps I -mue E  but continuously
 clear all 
 global r K Kp gami game eps mue mui rhosW rhoiW  
 formatSpecF = '%6.2f\n';
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);
 
 % specify parameters
  
b=0.6;% per year basic deer birth rate
d = 0.1;
r=b-d;
 
K = 30;  % effective carrying capacity % per 100 km^2 deer natural carrying capacity
Kp=b*K/(b-d);  % zero birth rate level
game = 1; %per year per prion mass (wlog)
gami = 0.1;
rhorat = 0.25; 
mui = 0.6; % per year infected deathrate % 
mue = 0.2;
cntrl=0;  % control variable for the deterministic model 

cntrlthresh = 30;  % threshhold for control
eps = 0.1 ; %prion mass/infected density %per year
 
 % now some simulations:
 A0 = 100; %km^2  This is the basic size for which K is measured
A = 2000;  % (km)^2 the area for which this simulation is being done
 
WRlist =  [0, 0.8, 1.5 ];  %these are the values used for rho_iW
  for nr = 1:  length(WRlist)
     rhoiW= WRlist(nr);
 rhosW = rhorat*rhoiW;
 %specify the output points
t_start = 0;
tstep =  0.1;
t_end = 100;
tspan = [t_start:tstep:t_end];
% first a deterministic simulation
s0 = [K;0.01;0];  %  format s0 = [U;V;...] as a column vector

[Td,S] = ode23s(@deRHS,tspan, s0, odeset('maxstep',1));  
X = S(:,1); 
Y = S(:,2);
E = S(:,3);
  
 figure(4*(nr-1)+1)
 nt = find(Td>=0);
 plot(Td(nt),S(nt)*A/A0,'b--', Td(nt),Y(nt)*A/A0,'r--')
 %axis([0 t_end 0 200])
 hold on

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now the stochastic simulation

Kt = 50;  % number of trials

%set the rate constants for reactions:
% reaction 1 birth of S  -> S+1 rate r
% reaction 2 death S   -> S-1 rate R
% reaction 3 infection S + E -> I  s->s-1, i-> i+1 rate game
% reaction 4 infection S + I --> 2I  s->s-1, i->i+1  rate gami
% reaction 5   death I -> 0  i->i-1  rate mui 
% reaction 6 wolf predation of susceptibles  s->s-1
% reaction 7 wolf predation of infecteds i->i-1
  
c(1) = b;  %s->s+1
R =  A0/(K*A);
c(2) = d; %s->s-1
c(3) = game; %S + E -> I  s->s-1, i->i+1
c(4) = gami *A0/A; %s->s-1, i->i+1
c(5) = mui;  %  i->i-1 
c(6) = rhosW ; %s->s-1 due to predation 
c(7) = rhoiW ;  %i->i-1 due to predation
 
 
%Specify the change matrix

Ch = [1,0,0;-1,0,0;-1,1,0;-1,1,0;0,-1,0; ...
    -1,0,1;0,-1,1];  %What happens to s, and i, when a reaction occurs:

%initialize the state space
clear Es Ss Is  T 

Kp = K*A/A0;

s = Kp *ones(Kt,1);
%Kp=K*A/A0 is the carrying capacity
  % start S at the carrying capacity
e = zeros(Kt,1); %start with no e
i = ones(Kt,1);  % there is exactly one introduced infected;
cd = zeros(Kt,1); %track total consumed deer
 cntrlon = zeros(Kt,1); % when this is one, there is predation

ct= zeros(Kt,1); % track total number infected by cwd
Ss = s ; %This   keeps track of the   trajectories
Is = i ;
Cntrlon = cntrlon;  %This saves the control variable
 T = zeros(Kt,1);  %This keeps track  of the transition times
j = 1; % count the number of reaction steps
rk = [];
 
 
 while (min(T(:,j))<t_end)  % make sure each simulation runs to at least t_end
    j = j+1;
    s = max(s,1); % this is to prevent s from going extinct
    % check to see if the disease has spread; if it has introduce wolves
   for jf = 1:Kt
       if(cntrlon(jf)==0)
           if (i(jf)>= cntrlthresh)
               cntrlon(jf)=1;
               cntrlsv(jf) = j; %this is the time step at which the change was made
           end
       end
   end


    % first calculate the maximal value of e for the future
    Estr = max(e,eps*i*A0/(A*mue));
    
    h(:,1) = max(c(1)*s.*(1-(s+i)*R),0) ; %reaction rate 1
    h(:,2) = c(2)*s; 
    h(:,3) = c(3)*s.*Estr; % Using Poisson thinning
    h(:,4) = c(4)*i.*s;
    h(:,5) = c(5)*i; % 
    h(:,6) = c(6)*cntrlon.*s;
    h(:,7) = c(7)*cntrlon.*i;
     
    hc = cumsum(h')'; % the cumulative sum of h
    H = sum(h')';
   
    rn = rand(Kt,2); %find 2 random numbers for each trajectory
    
    T(:,j) = T(:,j-1);  % add the current time to T
    
     delt=- log(rn(:,1))./H;
    T(:,j) =delt+T(:,j-1); % time of next reaction
    
    % use current value of i to update the e concentration 
    e  = eps*i*A0/(A*mue)+(e -eps*i*A0/(A*mue)).*exp(-mue*delt);
    for k = 1:Kt
    rk  = min(find(rn(k,2) <=hc(k,:)/H(k))); % this determines which reaction occurs
    if (rk  == 3)  %for Poisson thinning
        pstr = e(k)/Estr(k); % check the probability of reacting
        rp = rand(1,1);
        if(rp<pstr)
            s(k) = s(k) + Ch(rk,1); % update s, and i
            i(k) = i(k) + Ch(rk,2);
            ct(k) = ct(k) + 1; % another infection by prions
        end
    else 
       
    s(k) = s(k) + Ch(rk,1); % update s and i
    i(k) = i(k) + Ch(rk,2);
    
     cd(k) = cd(k) +Ch(rk,3);  %predation by wolves
    end
    end
     
    % save the values of the  trajectories
    Ss(:,j) = s ;
    Is(:,j) = i ;
    Es(:,j) = e;
    Cntrlon(:,j) = cntrlon;
    Cds(:,j)=cd;
    
    end
  
  kj = Kt; % use the first several trials to plot   sample trajectories
  
  % now make some plots
  %time course of one trajectory
  figure(4*(nr-1)+1)
  tp=find(T(1,:)<t_end);
  j=1;
 plot(T(j,tp),Ss(j,tp),'b' ,T(j,tp),Is(j,tp) ,'r',T(j,tp),Kp*Cntrlon(j,tp)/10,'--','linewidth',2)
title( strcat('\rho_iW = ',sprintf(formatSpecF,rhoiW)),'fontsize',18)
 
 xlabel('t(years)')  
 hold off
 
 %time course of E(T)
  figure(4*(nr-1)+2)
  % E(t)
 for j = 1:kj 
 semilogy(T(j,:) ,Es(j,:) ,'linewidth',2)
 hold on
 end
 axis([0 t_end 0 5])
 hold off

 xlabel('t (years)', 'fontsize', 20)
 ylabel('E', 'fontsize', 20)
title('E Trajectories','fontsize',20)
title( strcat('\rho_iW = ',sprintf(formatSpecF,rhoiW)),'fontsize',18)
 hold off
 
 endx=find(Es(:,end)>.001);
 nprion(nr)=length(endx);
 Mprion(nr) = mean(Es(endx,end));
 Ww(nr)=rhoiW;
    
% 
 
% phase portrait of trajectories
 figure(4*(nr-1)+3)
 for j = 1:kj
 plot(Ss(j,tp),Is(j,tp),'linewidth',2)
 hold on
 end
 %plot(Ss(1,cntrlsv(1)) ,Is(1,cntrlsv(1)) ,'*','linewidth',3)
 xlabel('s','fontsize',20)
 ylabel('i','fontsize',20)
 title( strcat('\rho_iW = ',sprintf(formatSpecF,rhoiW)),'fontsize',18)
 axis([0 Kp 0 350])
 hold off

 
  % find the disease survival probability
kt0=max(find(Td<=0));
ktmx = length(Td);
ktd=ktmx-kt0;

  for kk = 1:ktd
      ne0=0;

      for j=1:Kt
         nn= max(find(T(j,:)<=Td(kk+kt0)));
         if(Es(j,nn)>1.e-4)
             ne0=ne0+1;
         end

      end
      pe0(kk) = ne0;
  end
  % plot the disease survival probability
  figure(25 )
  plot(Td(kt0+1:ktmx),pe0/Kt,'*')
hold on
  end
  legend('boxoff')
  legend('\rho_iW=0','\rho_iW=0.8','\rho_iW=1.5','location','southwest')
  xlabel('t')
  ylabel('Probability of disease survival')
 
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %the right hand side for ode simulation:
function s_prime=deRHS(t,s)
global r K gami game eps mue mui rhoiW   rhosW    
S = s(1);
II = s(2);
E = s(3);
 
 fS= r*S*(1- (S+II)/K)-game*S*E -gami*S*II-  rhosW*S ;
                
fII=game*S*E +gami*S*II-mui*II- rhoiW*II ;
                  
fE=eps*II-mue*E;
 

 s_prime = [fS,fII, fE  ]';  

end   