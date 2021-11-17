clear all
close all

load('AllData');
%load('AllParameters');
time = {'P60'};

Mutants = {'T18C';'T25C';'T29C'};
%mutant = 'T18C';%'Atg6_mutant';%'Control';

%       time prob number
N = zeros(3,1);
T = [];
AllData = [];
for i = 1:length(Mutants)
    mutant = char(Mutants(i));
    Data = Output.(mutant).P60.Bulb.SurvProb;
    edges = Data(:,1)';
    T = [T,edges];
    AllData = [AllData;Data(:,2)];
    N(i) = length(Data);
    figure(i)
    hold on
    stairs(edges,Data(:,2),'Linewidth',3)
end
N = cumsum(N);

options = optimset('TolFun',rand,'TolX',1e-10,'Display','iter');
LB = zeros(7,1);
UB = ones(7,1);
% 2 compartments
%decline slow
UB(1) = 0.01;
LB(2) = 0.01;


r = rand;%*0.0001;
[y,res] = fminsearch(@fit2,[4 10*r rand rand rand],options,AllData',T,LB,UB,N);%,VarData);

c5 = LB(1)+ (UB(1)-LB(1)) * 1./(1+exp(y(1)));

c2 = LB(2)+ (UB(2)-LB(2)) * 1./(1+exp(y(2)));
prob_synBulb = zeros(3,1);
prob_synBulb(1) = LB(5)+ (UB(5)-LB(5)) * 1./(1+exp(y(3)));
prob_synBulb(2) = LB(6)+ (UB(6)-LB(6)) * 1./(1+exp(y(4)));
prob_synBulb(3) = LB(7)+ (UB(7)-LB(7)) * 1./(1+exp(y(5)));

c5
c2
prob_synBulb

counter = 1;
for i = 1:length(Mutants)
    t = T(counter:N(i));
    pred1 = targetfun(c5,t);
    pred2 = targetfun(c2,t);
    pred = prob_synBulb(i).*pred1 + (1-prob_synBulb(i)).*pred2; 
    figure(i) 
    h(1) = plot(t,prob_synBulb(i).*pred1,'k--','LineWidth',2);
    h(2) = plot(t,(1-prob_synBulb(i)).*pred2,'k-.','LineWidth',2);
    h(3) = plot(t,pred,'b:','LineWidth',3);
    xlim([0 60])
    ylim([0 1])
    xlabel('lifetime (in min)','FontSize',18)
    ylabel('Probability','FontSize',18)
    title(strcat('Temperature = ',char(Mutants(i)),''),'FontSize',18)
    legend(h,'long-lived bulbs','short-lived bulbs','sum','FontSize',14)
    set(gca,'Fontsize', 14)
    if i < 3
        counter = N(i)+1;
    end
end


function pred = targetfun(k,t)
    pred = exp(-t.*k);
end

function res = fit2(x,Data,T,LB,UB,N)%,VarData);

c5 = LB(1)+ (UB(1)-LB(1)) * 1./(1+exp(x(1)));
c2 = LB(2)+ (UB(2)-LB(2)) * 1./(1+exp(x(2)));
prob_synBulb = zeros(3,1);
prob_synBulb(1) = LB(5)+ (UB(5)-LB(5)) * 1./(1+exp(x(3)));
prob_synBulb(2) = LB(6)+ (UB(6)-LB(6)) * 1./(1+exp(x(4)));
prob_synBulb(3) = LB(7)+ (UB(7)-LB(7)) * 1./(1+exp(x(5)));

AllPred = zeros(1,length(T));
counter = 1;
for i = 1:3
    t = T(counter:N(i));
    pred1 = targetfun(c5,t);
    pred2 = targetfun(c2,t);
    pred = prob_synBulb(i).*pred1 + (1-prob_synBulb(i)).*pred2;  
    AllPred(counter:N(i)) = pred;
    if i < 3
        counter = N(i)+1;
    end
end

res = sum(((AllPred-Data).^2));%./VarData);
end

