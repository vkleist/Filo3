%growthConeSim
close all
clear all

%Gillespie ensemble setup
N_gil = 50;
T = 3600;%final time (60 hours = 3600 minutes = 216000 seconds)

%Bulbous tips need to be at least 40 min old to be synaptogenic
t_thres = 0;%40;%min

%FeedbackOption: the following options are available: 'FB','wtFB','noFB','wtr5'; 
%first is 'product inhibition of the mutant', second is 'product inhibition
%with wt parameter B50', third is 'no feedback', last is feedback with wt
%parameters B50 and k
FeedBackOption = 'wtFB';%
mutant = 'T25C';%'T18C'; 'T25C';'T29C'

% developmental to real time 	factor
Scale.T18C = 2.046454768;
Scale.T25C = 1;
Scale.T29C = 0.897310513;

T = round(T.*Scale.(mutant));

Type = mutant;
file = strcat('EnsembleData/DataSimple_',mutant);
disp('Attention: rate Bulb -> Syn changed to simple first order reaction, line 220')

% Options for saving ensemble data
time_disc = 1;%take a snapshot every minute
time_points = ceil(T/time_disc);
ens_Data  = zeros(N_gil,time_points,3); %ensemble data filopodia

[c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, p, thalf, h,c6,c7] = getParameters(mutant,FeedBackOption);


% Stoichiometric matrix         <-r_null is for numerical correctness
%    r_null r1 r2   r3    r4    r5    
S = [0      1   -1  0     -1     0;%sB
     0      0   0   -1    1      0;%LB
     0      0   0   0     0      1];%S

%start of Gillespie run
for i=1:N_gil
tic
t = 0; %start time

%actual state
X = zeros(3,1);
sBul_LifeTimes = nan(100,1);

%setup storage
counter2 = 0;
SizeStoreVector = 5000;
Time = nan(1,SizeStoreVector);
XStore = nan(length(X),SizeStoreVector);

% Update state
counter2= counter2+1;
Time(counter2) = 0;
XStore(:,counter2) = X;

counter = 0;
%start of growth cone simulation
while t < T
    
    %compute propensity and their sum
    Vec = sBul_LifeTimes(~isnan(sBul_LifeTimes));
    if isempty(Vec)
        Vec = -inf;
    end
    a = getPropensities(t,X,c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, p, thalf, h,t_thres,Vec,mutant,Scale.(mutant),c6,c7);
    a0 = sum(a);
    % 1. time update
    r1 = rand;
    tau = (1/a0)*log(1/r1);%time step
    t = t+tau; %time update
    % 2. sample reaction
    r2 = rand;
    j = find(r2 <= cumsum(a)/a0,1);
    
    X_before = X;
  
    % 3. state update
    X = X + S(:,j);
    
       
    
    % ---- new block------ 
    % reaction R4a, sB died
    if j == 6 
        idx = find(~isnan(sBul_LifeTimes),1,'last');
        sBul_LifeTimes(idx) = nan;%update lifetimes
    end
    
    %Update Bulbous Lifetimes
    sBul_LifeTimes = sBul_LifeTimes+tau;
    
    %emergence of a new bulbous; reaction R3, initialize its lifetime
    if j == 1
        idx = find(isnan(sBul_LifeTimes),1);
        sBul_LifeTimes(idx) = 0;
    end
    
    %check whether sB -> LB
    if j == 4
        Vec = sBul_LifeTimes(~isnan(sBul_LifeTimes));
        r3 = rand;
        idx2 = find(~isnan(sBul_LifeTimes));
        z = find(r3 <= cumsum(Vec)/sum(Vec),1);
        sBul_LifeTimes(idx2(z)) = nan;
    end

    %----------

    
    %time data for ensemble matrices
    t_now = ceil(t);
    t_before = ceil(t-tau);
    t_passed = t_before:time_disc:t_now-1;
    
    if t_now > T
        break;
    end
    
    %Matrices updates
    counter2 = counter2+1;
    if counter2 > SizeStoreVector
        Time = [Time,nan(1,5000)];
        XStore= [ XStore,nan(5,5000)];
        SizeStoreVector = SizeStoreVector+5000;
    end
    Time(counter2) = t;
    XStore(:,counter2) = X;
    
    %update ensemble matrices
    ens_Data(i,t_now,3) = X(1);
    ens_Data(i,t_now,4) = X(2);
    ens_Data(i,t_now,5) = X(3);

    
    if t_before <= 1
        continue;
    end
    
    %if t_passed isn't 0, fill in data for the time between
    aux = ens_Data(i,t_before,1);
    ens_Data(i,t_passed,1) = aux;
    %
    aux = ens_Data(i,t_before,2);
    ens_Data(i,t_passed,2) = aux;
    %
    aux = ens_Data(i,t_before,3);
    ens_Data(i,t_passed,3) = aux;
    %
%     aux = ens_Data(i,t_before,4);
%     ens_Data(i,t_passed,4) = aux;
%     %
%     aux = ens_Data(i,t_before,5);
%     ens_Data(i,t_passed,5) = aux;
    

    
end

toc

%fill in last bits of ensemble matrices
t_before = ceil(Time(counter2));
t_passed = t_before:time_disc:T;

aux = ens_Data(i,t_before,1);
ens_Data(i,t_passed,1) = aux;
%
aux = ens_Data(i,t_before,2);
ens_Data(i,t_passed,2) = aux;
%
aux = ens_Data(i,t_before,3);
ens_Data(i,t_passed,3) = aux;
%
% aux = ens_Data(i,t_before,4);
% ens_Data(i,t_passed,4) = aux;
% %
% aux = ens_Data(i,t_before,5);
% ens_Data(i,t_passed,5) = aux;
end


%TODO: include T, t_disc (tau),...
time_out = (1:T)./(60*Scale.(mutant))+40;
try save(file,'ens_Data','time_out')
catch
    mkdir('./EnsembleData/');
    save(file,'ens_Data','time_out')
end
    
%Trimming
Time = Time(1:counter2);
XStore = XStore(:,1:counter2);

plotGrowthConeSimSimple_BulbLifeTimeThreshold(Time,XStore,time_points,ens_Data,mutant,Type,Scale.(mutant));


function a = getPropensities(t,X,c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, p, thalf, h,t_thres,Vec,mutant,Scale,c6,c7)

    a = zeros(1,6);
    a(1)  = c1_sF*time_dampening(t/Scale); %birth of short lived filo
    a(2)  = c3*(time_dampening(t/Scale)./time_dampening(20*60))*time_enhancing(t/Scale,thalf,h)./time_enhancing(20*60,thalf,h)*feedback(B50,X(2)); % Filopodium-to-bulbous transition
    a(3)  = c4*X(1); % short-bulbous death
    a(4)  = c7.*X(2);%1/120; % long-bulbous death
    a(5)  = X(1)*c5;%short-to-long bulb
    a(6)  = c6*X(2);%1/(100*Scale)*X(4); % bulbous to synapse
   
    
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

function [c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, p, thalf, h,c6,c7] = getParameters(param,FeedBackOption)

load('AllParametersSimple3');

c1_sF = Parameters.(char(param)).c1_sF;%birth           
c1_LF = Parameters.(char(param)).c1_LF;%birth
c2_sF = Parameters.(char(param)).c2_sF;%death
c2_LF = Parameters.(char(param)).c2_LF;%death
%bulbous dynamics
c3 = Parameters.(char(param)).c3;%rF2B
B50 = Parameters.(char(param)).B50;%rB2S
c4 = Parameters.(char(param)).c4;%rB20
c5 = Parameters.(char(param)).c5;%rB2S
c6 = Parameters.(char(param)).c6;%rB2S
c7 = Parameters.(char(param)).c7;%rB2S
% time dependent functions
% a) dampening of filo dynamics
p = Parameters.(char(param)).p;
% b) increase in propensity to form bulbous
thalf = Parameters.(char(param)).thalf;
h = Parameters.(char(param)).h;

end


  