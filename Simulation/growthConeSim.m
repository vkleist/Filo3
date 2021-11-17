%growthConeSim
close all
clear all

%Gillespie ensemble setup
N_gil = 50;
T = 3600;%final time (60 hours = 3600 minutes = 216000 seconds)


temperature = 'T29C';%'T18C'; 'T25C';'T29C'

% developmental to real time factor
Scale.T18C = 2.046454768;
Scale.T25C = 1;
Scale.T29C = 0.897310513;

T = round(T.*Scale.(temperature));

file = strcat('EnsembleData/DataSimple_',temperature);

% Options for saving ensemble data
time_disc = 1;%take a snapshot every minute
time_points = ceil(T/time_disc);
ens_Data  = zeros(N_gil,time_points,3); %ensemble data filopodia

[null1,null2,r1, B50, c2, c3, thalf, h,c4,c5] = getParameters(temperature);% 

% Stoichiometric matrix  <--Changed: eliminated all reaction regarding the filopodia
%  null null r1     r2      r5    r3    r4    
S = [0   0   1      -1      0     -1     0;%sB
     0   0   0      0       -1    1      0;%LB
     0   0   0      0       0     0      1];%S

%start of Gillespie run
for i=1:N_gil
tic
t = 0; %start time

%Initialize filopodia/stab filo/bilbous/synapse
sF = 0; %short lived filopodia
LF = 0;%long lived folopodia
%actual state
X = zeros(3,1);
%sBul_LifeTimes = nan(100,1);

%setup storage

% Update state
counter2= 1;

%start of growth cone simulation
while t < T
    
    %compute propensity and their sum

    a = getPropensities(t,X,null1,null2,r1, B50, c2, c3, thalf, h,Scale.(temperature),c4,c5);
    a0 = sum(a);
    % 1. time update
    tau = (1/a0)*log(1/rand);%time step
    t = t+tau; %time update
    % 2. sample reaction
    j = find(rand <= cumsum(a)/a0,1);
    
    X_before = X;
    if X(3) < 0
        disp('before update')
        return
    end
  
    % 3. state update
    X = X + S(:,j);
    
    if X(3) < 0
        return
    end

    %time data for ensemble matrices
    t_now = ceil(t);
    t_before = ceil(t-tau);
    t_passed = t_before:time_disc:t_now-1;
    
    if t_now > T
        break;
    end
    
    %Matrices updates
    counter2 = counter2+1;
    
    %update ensemble matrices
    ens_Data(i,t_now,1) = X(1);
    ens_Data(i,t_now,2) = X(2);
    ens_Data(i,t_now,3) = X(3);

    
    if t_before <= 1
        continue;
    end
    
    aux = ens_Data(i,t_before,1);
    ens_Data(i,t_passed,1) = aux;
    %
    aux = ens_Data(i,t_before,2);
    ens_Data(i,t_passed,2) = aux;
    %
    aux = ens_Data(i,t_before,3);
    ens_Data(i,t_passed,3) = aux;
    

    
end

toc

%fill in last bits of ensemble matrices
t_before = ceil(t-tau);
t_passed = t_before:time_disc:T;

%
aux = ens_Data(i,t_before,1);
ens_Data(i,t_passed,1) = aux;
%
aux = ens_Data(i,t_before,2);
ens_Data(i,t_passed,2) = aux;
%
aux = ens_Data(i,t_before,3);
ens_Data(i,t_passed,3) = aux;
end

time_out = (1:T)./(60*Scale.(temperature))+40;
save(file,'ens_Data','time_out');

plotSimulation(time_points,ens_Data,temperature,Scale.(temperature));


function a = getPropensities(t,X,null1,null2,r1, B50, c2, c3, thalf, h,Scale,c4,c5)

    a = zeros(1,7);
    a(1)  = null1; %birth of short lived filo
    a(2)  = null2; %birth of long lived filo
    
    a(3)  = r1*(time_dampening(t/Scale)./time_dampening(20*60))*time_enhancing(t/Scale,thalf,h)./time_enhancing(20*60,thalf,h)*feedback(B50,X(2)); % Filopodium-to-bulbous transition
    a(4)  = c2*X(1); % short-bulbous death
    a(5)  = c5.*X(2);% % long-bulbous death
    a(6)  = X(1)*c3;%short-to-long bulb
    a(7)  = c4*X(2); % bulbous to synapse  
end

%Function that influences the enhancement of bulbous birth
function enhance = time_enhancing(t,thalf,h)
scale = 2^(-h);
enhance = scale.*(1+tanh(3./thalf.*(t-thalf))).^(h);
end

%Function that influences the damepning of filopodia birth
function damp = time_dampening(t)
    p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
    damp = polyval(p,t);
    damp = max(1e-4,damp);
end

function f = feedback(B50,B)
    f = B50./(B50 + B);
end

function [null1,null2,r1, B50, c2, c3, thalf, h,c4,c5] = getParameters(param)

load('AllParameters');

null1 = Parameters.(char(param)).null1;%birth           
null2 = Parameters.(char(param)).null2;%birth
%bulbous dynamics
r1 = Parameters.(char(param)).r1;%-> sB
B50 = Parameters.(char(param)).B50;%-> sB
c2 = Parameters.(char(param)).c2;%sB -> 0
c3 = Parameters.(char(param)).c3;%sB -> synB
c4 = Parameters.(char(param)).c4;%synB -> synB + synapse
c5 = Parameters.(char(param)).c5;%synB -> 0
% b) increase in propensity to form bulbous
thalf = Parameters.(char(param)).thalf;
h = Parameters.(char(param)).h;

end


  