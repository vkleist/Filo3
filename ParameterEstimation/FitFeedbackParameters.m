close all
clear all

Mutants =  {'T18C';'T25C';'T29C'};
mutant = 1;
MutantName = strrep(char(Mutants(mutant)),'_','-');
Colors = [.5 .5 .5;0 0 1;1 0 0;0 1 0];
Ref_yellow = [0.9 0.8 0.2];%[0.9100 0.4100 0.1700];

FSP = 20; % finite state projection
n = 6;

c = Colors(mutant,:);

Mutants(mutant)

% % developmental to real time 	factor
Scales.T18C = 2.046454768;
Scales.T25C = 1;
Scales.T29C = 0.897310513;

Scale = Scales.(char(Mutants(mutant)));

t = 20*60;%P60
%
halfmax = 1000;
exponent = 1;
enhance = enhancing(t,halfmax,exponent);
%
p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
damp = dampening(p,t);

StartParamFB = log([0.022*7 0.1133 0.0282]);%[1/rand rand rand];%

%-------created by EstimateBulbLifeTimesNew
                   %18°C            25°C            29°C    
DeathRate_sB    = [0.0706           0.0706          0.0706];%death rate of short-lived bulbous; c2
DeathRate_synB  = [0.0060           0.0060          0.0060]; %death rate of long-lived bulbous; c5

%

d_sB = DeathRate_sB(mutant);
d_synB = DeathRate_synB(mutant);

%-------created by "PlotBulbousTipsPerLife....m"
                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_synBulb = [0.012  0.28    0.59   0.094    0.025  0       0; %18°C
                      0.15   0.56    0.24   0.04     0      0       0; %25°C
                      0.44   0.56    0.0065 0        0      0       0]; %29°C  

          
                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_sBulb = [0.76   0.24     0.0029  0       0       0       0; %18°C
                    0.69   0.27     0.033   0       0       0       0; %25°C
                    0.62   0.35     0.035   0       0       0       0];%29°C
       
                    
                % if LifeTimeThreshold = 0 in counting routine
                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_AllBulb = [0.012 0.22    0.51    0.21     0.038  0.013    0; %Control                                               
                      0.076 0.51    0.3     0.082    0.022  0.013    0;  %Atg6_mutant                                          
                      0.23  0.56    0.2     0.01     0      0        0];  %Atg7_mutant
                      
%
edges = 0:7;
x_pos = 0:n;
Data_synBulb = ExpDensity_synBulb(mutant,:);
Data_sBulb = ExpDensity_sBulb(mutant,:);
Data = [Data_sBulb , Data_synBulb];

figure(1)
hold on
bar(edges(1:end-1),ExpDensity_synBulb(mutant,:)','EdgeColor','none','FaceColor',c);
hold on

figure(2)
hold on
bar(edges(1:end-1),ExpDensity_sBulb(mutant,:)','EdgeColor','none','FaceColor',c);
hold on

figure(3)
hold on
bar(edges(1:end-1),ExpDensity_AllBulb(mutant,:)','EdgeColor','none','FaceColor',c);
hold on


%% Data fitting feedback (product inhibition on f1 or f2)
    
StartParam = StartParamFB;%[r,rand,rand];

options = optimset('Display','iter');
x = fminsearch(@fitfun,StartParam,options,Data,FSP,d_sB,d_synB,n,mutant);
   
c3_P60 = exp(x(1))
c5 = exp(x(2))
B50 = exp(x(3))

 
L = makeGenerator(d_sB,d_synB,c3_P60,B50,FSP,c5);

maxIDx = length(Data);

M = L';
[V,D] = eig(M);
idx = find(round(diag(D).*1e5)./1e5 == 0);
statdistr = V(:,idx)./sum(V(:,idx));
TupelM = Idx2Tupel(FSP);
%add together sBulb + synBulb
P1 = zeros(n+1,1);
P2 = zeros(n+1,1);
for j = 0:n
    idx1 = TupelM(:,1) == j;
    P1(j+1) = sum(statdistr(idx1));
    idx2 = TupelM(:,2) == j;
    P2(j+1) = sum(statdistr(idx2));
end

%prob sB
idx = TupelM(:,1)>0;
Prob_sB = sum(statdistr(idx));

%prob synB
idx = TupelM(:,2)>0;
Prob_synB = sum(statdistr(idx));

figure(1)
hold on
%stairs(x_pos-0.5,P2,'--','Linewidth',4,'Color','r')
stairs(x_pos-0.5,P2,'-','Linewidth',4,'Color',[0.2 0.2 0.2])
%if mutant == 1
    legend('data','no feedback', 'feedback')
%end
Path = './Figures/';
Name = strcat('Fig_BulbousFit_Synaptogenic',MutantName,'.tiff');
print(1,'-dtiff',strcat(Path,Name,'.tiff'))

figure(2)
hold on
%stairs(x_pos-0.5,P1,'--','Linewidth',4,'Color','r')
stairs(x_pos-0.5,P1,'-','Linewidth',4,'Color',[0.2 0.2 0.2])
%if mutant == 1
    legend('data','no feedback', 'feedback')
%end
Path = './Figures/';
Name = strcat('Fig_BulbousFit_Transient',MutantName,'.tiff');
print(2,'-dtiff',strcat(Path,Name,'.tiff'))

Prob_allB = nan(n+1,1);
%prob all Bulbs
for j = 0:n
    idx = sum(TupelM,2)==j;
    Prob_allB(j+1) = sum(statdistr(idx));
end
figure(3)
stairs(x_pos-0.5,Prob_allB,'-','Linewidth',5,'Color',[0.2 0.2 0.2])
%bar(edges(1:end-1),[ExpDensity(mutant,:)',P]);
set(gca,'FontSize',18)
ylabel('probability','FontSize',20)
xlabel('number bulbous tips/(GC and time)','FontSize',20)
title(strcat(MutantName,': all bulbs'),'FontSize',20)
%if mutant == 1
    legend('data','fit','FontSize',18);%, 'feedback')
% % %end
Path = './Figures/';
Name = strcat('Fig_BulbousFit_All',MutantName);
print(3,'-dtiff',strcat(Path,Name,'.tiff'))
print(3,'-depsc2',strcat(Path,Name,'.eps'))


function res = fitfun(x,Data,FSP,d_sB,d_synB,n,mutant)

    c3 = exp(x(1));%
    c5 = exp(x(2));%
    B50 = exp(x(3));%

    L = makeGenerator(d_sB,d_synB,c3,B50,FSP,c5);

    M = L';
    [V,D] = eig(M);
    idx = find(round(diag(D).*1e5)./1e5 == 0);
    statdistr = V(:,idx)./sum(V(:,idx));
    
    % add together sBulb + synBulb
    TupelM = Idx2Tupel(FSP);
    % add together sBulb + synBulb
    P1 = zeros(1,n+1);
    P2 = zeros(1,n+1);
    for j = 0:n
        idx1 = TupelM(:,1) == j;
        P1(j+1) = sum(statdistr(idx1));
        idx2 = TupelM(:,2) == j;
        P2(j+1) = sum(statdistr(idx2));
    end
    
    res = KL(Data,[P1, P2]);
end

function L = makeGenerator(d_sB,d_synB,c3,B50,FSP,c5)
    L = zeros((FSP+1)^2,(FSP+1)^2);
    
    for i = 0:FSP % number short lived bulb
        for j = 0:FSP %number long lived buld
            idx = i*(FSP+1)+(j+1);%current state
            
            if i > 0 %sB -> *
                idx_death_i = (i-1)*(FSP+1)+(j+1);%(i-1,j)
                L(idx,idx_death_i) = d_sB*i;
            end
            
            if j > 0 %synB -> *
                idx_death_j = (i)*(FSP+1)+(j);%(i,j-1)
                L(idx,idx_death_j) = (d_synB)*j;
            end
            
            if i < FSP %F -> sB +1
                idx_birth_i = (i+1)*(FSP+1)+(j+1);%(i+1,j)
                L(idx,idx_birth_i) = c3.*feedback(B50,j);%c3.*feedback(B50,(i+j)); %c3.*feedback(B50,j);%
            end
            
            if j < FSP && i > 0% sB -> synB
                idx_birth_j = (i-1)*(FSP+1)+(j+2);%(i-1,j+1)
                L(idx,idx_birth_j) = (i)*c5;%
            end

        end
    end
    L = L - diag(sum(L,2));
end


function f = feedback(B50,B)
    f = B50./(B50 + B);
end

function p = poisspdf(x,lambda)
    p = lambda.^x./factorial(x).*exp(-lambda); 
end

function TupelM = Idx2Tupel(n)
    TupelM = nan((n+1)^2,2);
    %This Matrix contains at row i the number of short-lived bulbs (first
    %column) and long lived bulbs (second column)
    counter = 0;
    for i = 0:n
        for j = 0:n
            counter = counter +1;
            TupelM(counter,1) = i;
            TupelM(counter,2) = j;
        end
    end
end

function enhance = enhancing(t,thalf,h)
scale = 2^(-h); 
enhance = scale.*(1+tanh(3./thalf.*(t-thalf))).^(h);
end

function d = dampening(p,t)
    d = polyval(p,t);
    d = max(1e-4,d);
end
