%modified Gillespie algorithm for S-I-E-W model 
close all; clear; clc

%updated 6-6 to fix carrying capacity of wolves and get r10 and r11 uniform
%with how the deer is implemented. Also changed eps to eps1. eps is not protected but rather unsafe to use in MATLAB.

rng(1)

global r K Kp gami game eps1 mue mui rho alp omega kmw K2 rw rwb rwd
formatSpecF = '%6.2f\n';
set(0,                           ...
    'defaultaxesfontsize', 20,   ...
    'defaultaxeslinewidth', 2.0, ...
    'defaultlinelinewidth', 2.0, ...
    'defaultpatchlinewidth', 0.7);

% specify parameters
d = 0.1;
b = 0.6;
r=b-d;
mui = 0.6;
mue=0.2;
game = 1;  %wlog
gami = 0.1;
eps1 = 0.1 ;
K = 30;  % effective carrying capacity
Kp=b*K/(b-d);  % zero birth rate level
alp=0.05;
K2=3;
rwb=0.5; %updated 6-6
rwd=0.2; %updated 6-6
rw=rwb-rwd; %updated 6-6
K2p=rwb*K2/(rwb-rwd);  % zero birth rate level for wolves 
rho=0.01;
%omega is our bifurcation parameter

% now some simulations:

A0 = 100; %km^2  This is the basic size for which K is measured
A = 2000;  % (km)^2 the area for which this simulation is being done

tic

WRlist =  [1:10:500]; %vary omega
for nr = 1:length(WRlist)
    omega=WRlist(nr)

    % --- Simulation Setup ---
    t_start = -15000;
    tstep = 0.1;
    t_end = 200;
    tspan = t_start:tstep:t_end;

    % Initial conditions [S; II; E; W]
    s0 = [K; 0.01; 0; 3];

    % --- Solver ---
    % Passing parameters using an anonymous function is cleaner than 'global'
    %Gemini was used to sort plotting bug
    [Td, S] = ode23s(@(t, s) deRHS(t, s, r, K, game, gami, eps1, mue, mui, rho, alp, rw, K2, omega), tspan, s0, odeset('maxstep', 1));

    % --- Extract Results ---
    X = S(:, 1); % Susceptible
    Y = S(:, 2); % Infected
    E = S(:, 3); % Environment/Prions
    W = S(:, 4); % Wolves/Predator

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now the stochastic version

    Kt = 200;  % number of trials

    %set the rate constants for reactions:
    % reaction 1 birth of sus. deer, s-> s+1 {carrying capacity folded in
    % here}
    % reaction 2 death of sus. deer, s->s-1 
    % reaction 3 indirect infection of sus. deer, s->s-1, i-> i+1
    % reaction 4 direct infection of sus. deer, s->s-1, i->i+1
    % reaction 5 death of inf. deer,  i->i-1
    % reaction 6 wolf predation of sus. deer (no impact on w), s->s-1
    % reaction 7 wolf predation of inf. deer (no impact on w), i->i-1
    % new reactions due to dynamic W
    % reaction 8 birth of wolf from predation of sus. deer w->w+1, s->s-1 {new}
    % reaction 9 birth of wolf from predation of inf. w->w+1, i->i-1 {new}
    % reaction 10 birth of wolf, w->w+1 {new} {carrying capacity folded in
    % here} 
    % reaction 11 death of wolf {new}

    %specify reaction rates
    c(1) = b;  %birth rate of S
    R =  A0/(Kp*A); %carrying capacity of S corrected for volume (second order reaction); Note: need to use Kp here so everything matches up.
    c(2) = d; %death rate of S
    c(3) = game; %indirect transmission rate
    c(4) = gami*A0/A; %direct transmission rate corrected for volume (second order reaction)
    c(5) = mui;  %mortality rate of I
    c(6) = (1-alp)*rho*A0/A ; %S predation rate corrected for volume (second order reaction), no change to W in this case.
    c(7) = (1-alp)*omega*rho*A0/A ; %I predation rate corrected for volume (second order reaction), no change to W in this case.
    c(8) = alp*rho*A0/A; %S predation rate corrected for volume (second order reaction), change to W in this case.
    c(9) = alp*omega*rho*A0/A; %I predation rate corrected for volume (second order reaction), change to W in this case.
    R2 =  A0/(K2p*A); %carry capacity for wolves corrected for volume (second order reaction) note that this changes to R2P exactly as the deer case. 
    c(10) = rwb; %birth rate of W
    c(11) = rwd; %death rate of W

    %Specify the change matrix {index = [s,i,total consumed deer,w]}
    Ch = [1,0,0,0;-1,0,0,0;-1,1,0,0;-1,1,0,0;0,-1,0,0; ...
        -1,0,1,0;0,-1,1,0;...
        -1,0,1,1; 0,-1,1,1; 0,0,0,1; 0,0,0,-1];

    %initialize the state space
    clear Es Ss Is  T Ws
    s = K*A/A0 *ones(Kt,1); %K*A/A0 is the carrying capacity, start s at the carrying capacity
    e = zeros(Kt,1); %start with no e
    i = ones(Kt,1);  % there is exactly one introduced infected;
    cd = zeros(Kt,1); %track total consumed deer, starts at 0.
    w = K2*A/A0*ones(Kt,1);  %set wolves to the carrying capacity

    ct= zeros(Kt,1); % track total number infected by cwd
    Ss = s ; %This   keeps track of the trajectories
    Is = i ;
    T = zeros(Kt,1);  %This keeps track  of the transition times
    j = 1; % count the number of reaction steps
    rk = [];


    while (min(T(:,j))<t_end)  % make sure each simulation runs to at least t_end
        j = j+1;
        s = max(s,1); % this is to prevent s from going extinct

        % first calculate the maximal value of e for the future
        Estr = max(e,eps1*i*A0/(A*mue));

        h(:,1) = max(c(1)*s.*(1-((s)+i)*R),0) ; %carrying capacity is folded in here. 
        h(:,2) = c(2)*s;
        h(:,3) = c(3)*s.*Estr; % Using Poisson thinning
        h(:,4) = c(4)*i.*s;
        h(:,5) = c(5)*i; %
        h(:,6) = c(6)*s.*w;
        h(:,7) = c(7)*i.*w;
        h(:,8) = c(8)*s.*w;
        h(:,9) = c(9)*i.*w;
        h(:,10) = max(c(10)*w.*(1-w*R2),0); %carrying capacity is folded in here. 
        h(:,11) = c(11)*w;

        hc = cumsum(h')'; % the cumulative sum of h
        H = sum(h')';

        rn = rand(Kt,2); %find 2 random numbers for each trajectory

        T(:,j) = T(:,j-1);  % add the current time to T

        delt=- log(rn(:,1))./H;
        T(:,j) =delt+T(:,j-1); % time of next reaction

        % use current value of i to update the e concentration
        e  = eps1*i*A0/(A*mue)+(e -eps1*i*A0/(A*mue)).*exp(-mue*delt);
        for k = 1:Kt
            if H(k)==0
                disp('warn H is zero')
            end
            rk  = min(find(rn(k,2) <=hc(k,:)/H(k))); % this determines which reaction occurs
            rk_list(j,k) = rk;
            if (rk  == 3)  %for Poisson thinning
                pstr = e(k)/Estr(k); % check the probability of reacting
                rp = rand(1,1);
                if(rp<pstr)
                    s(k) = s(k) + Ch(rk,1); % update s, and i
                    i(k) = i(k) + Ch(rk,2);
                    ct(k) = ct(k) + 1; % another infection by prions
                end
            else

                s(k) = s(k) + Ch(rk,1); % update s
                i(k) = i(k) + Ch(rk,2); % update i
                w(k) = w(k) + Ch(rk,4); % update w
                cd(k) = cd(k) +Ch(rk,3);% update total deer consumed
            end
        end

        % save the values of the  trajectories
        Ss(:,j) = s ;
        Is(:,j) = i ;
        Es(:,j) = e;
        Ws(:,j) = w;

        Cds(:,j)=cd;

    end

    kj = Kt; % use the first several trials to plot   sample trajectories

    %disp('This is the number of trials with E < 0.001')

    erad_num(nr)=length(find(Es(:,end)<0.001 & Ws(:,end)>0 & Is(:,end)==0))/Kt;  %keep track of the number of trials that hit the disease-free state (i=E=~0), s and w nonzero.


    % %plot a single trajectory of S and I
    % 
    % figure(6*(nr-1)+1)
    % tp=find(T(1,:)<t_end);
    % j=1;
    % plot(T(j,tp),Ss(j,tp),'b' ,T(j,tp),Is(j,tp) ,'r','linewidth',2)
    % hold on
    % plot(Td, X*A/A0, 'b--', 'DisplayName', 's'); hold on;
    % plot(Td, Y*A/A0, 'r--', 'DisplayName', 'i');
    % title( strcat('\omega = ',sprintf(formatSpecF,omega)),'fontsize',18)
    % axis([0 200 0 500])
    % 
    % xlabel('t(years)')
    % 
    % hold off
    % 
    % alt_ss_s = -(K*rw*(K2*rho-r))/(K*K2*alp*rho^2+r*rw)*A/A0; 
    % alt_ss_i=0; 
    % 
    % %all the s, i, trajectories
    % figure(6*(nr-1)+2)
    % for j = 1:kj
    %     tp=find(T(j,:)<t_end);
    %     %j=1:1;
    %     plot(T(j,tp),Ss(j,tp),'Color',[0,0,1,0.1], 'lineWidth',2)
    %     hold on
    %     plot(T(j,tp),Is(j,tp),'Color',[1,0,0,0.1],'linewidth',2)
    %     plot(Td, X*A/A0, 'b--', 'DisplayName', 's');
    %     plot(Td, Y*A/A0, 'r--', 'DisplayName', 'i');
    %     plot(Td, alt_ss_s*ones(length(Td),1), 'Color','y', 'LineStyle','-.', 'DisplayName', 's');
    %     plot(Td, alt_ss_i*ones(length(Td),1), 'r', 'LineStyle','-.','DisplayName', 'i');
    %     title( strcat('\omega = ',sprintf(formatSpecF,omega)),'fontsize',18)
    % end
    % axis([0 200 0 600])
    % hold off
    % 
    % 
    % %plot many trajectories of E
    % figure(6*(nr-1)+3)
    % for j = 1:kj
    %     semilogy(T(j,:) ,Es(j,:) ,'linewidth',2)
    %     hold on
    %     plot(Td, E, 'cyan','LineStyle','--','LineWidth',4,'DisplayName', 'w');
    % end
    % axis([0 t_end 1e-20 2])
    % hold off
    % 
    % xlabel('t (years)', 'fontsize', 20)
    % ylabel('E', 'fontsize', 20)
    % %title('E Trajectories','fontsize',20)
    % title( strcat('\omega = ',sprintf(formatSpecF,omega)),'fontsize',18)
    % hold off
    % 
    % % phase portrait (s vs i).
    % figure(6*(nr-1)+4)
    % plot(Ss(1,tp)*A0/A,Is(1,tp)*A0/A,X,Y,'--','linewidth',2)
    % xlabel('S','fontsize',20)
    % ylabel('I','fontsize',20)
    % title( strcat('\omega = ',sprintf(formatSpecF,omega)),'fontsize',18)
    % 
    % 
    % %single trajectory of w
    % 
    % figure(6*(nr-1)+5)
    % tp=find(T(1,:)<t_end);
    % j=1;
    % plot(T(j,tp),Ws(j,tp),'b','linewidth',2)
    % hold on
    % plot(Td, W*A/A0, 'k', 'DisplayName', 'w');
    % %plot(Td, alt_ss_w*ones(length(Td),1), 'k','LineStyle','-.', 'DisplayName', 'w');
    % title( strcat('\omega = ',sprintf(formatSpecF,omega)),'fontsize',18)
    % axis([0 200 0 800])
    % %axis([0 t_end 0 500])
    % 
    % alt_ss_w = K2*r*(K*alp*rho + rw)/(K*K2*alp*rho^2 + r*rw)*A/A0; 
    % 
    % %all trajectories of w
    % 
    % figure(6*(nr-1)+6)
    % for j = 1:kj
    %     tp=find(T(j,:)<t_end);
    %     %j=1:1;
    %     plot(T(j,tp),Ws(j,tp),'Color',[0,0,1,0.1])
    %     hold on
    %     plot(Td, W*A/A0, 'k', 'DisplayName', 'w');
    %     plot(Td, alt_ss_w*ones(length(Td),1), 'Color','y','LineStyle','-.', 'DisplayName', 'w');
    %     title( strcat('\omega = ',sprintf(formatSpecF,omega)),'fontsize',18)
    % end
    % axis([0 200 0 800])
    % hold off

end

figure(25)
plot(WRlist,1-erad_num,'*')
xlabel('\omega')
ylabel('Probability of disease survival')

toc



%ODE simulation
function s_prime = deRHS(t, s, r, K, game, gami, eps1, mue, mui, rho, alp, rw, K2, omega)
S = s(1);
II = s(2);
E = s(3);
W = s(4);

fS = r*S*(1 - (S+II)/K) - game*S*E - gami*S*II - rho*S*W;
fII = game*S*E + gami*S*II - mui*II - rho*omega*II*W;
fE = eps1*II - mue*E;
fW = rw*W.*(1- W/K2)+alp*(rho*S*W + omega*rho*II*W);

s_prime = [fS; fII; fE; fW];
end