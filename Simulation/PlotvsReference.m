%function plotGrowthConeSimSimple_BulbLifeTimeThreshold(Time,XStore,time_points,ens_Data,param_set,AutophagyType)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all
clear all

path = './Figures/';
try cd (path)
    cd ..
catch
    mkdir(path)
end

FigureName = 'AdjustedRate_BulbToSynapse';

Mutants = {'T18C';'T25C';'T29C'};
MutantNames = {'T18C';'T25C';'T29C'};
Scales.T18C = 2.046454768;
Scales.T25C = 1;
Scales.T29C = 0.897310513;

Colors = [.5 .5 .5;0 0 1;1 0 0;1 0 1];



for i = 1:length(Mutants)
    mutant = char(Mutants(i));
    MutantName = char(MutantNames(i))
    Filename = strcat('./EnsembleData/DataSimple_',mutant,'.mat');
    S = load(Filename);
    PlotColor = Colors(i,:);
    
    T = length(S.ens_Data);
    XPositions = (1:100:T);
    Scale = Scales.(mutant);
    XPositionsPlot = XPositions./(60*Scale)+40;

    
    %Bulbs
    Ref1 = S.ens_Data(:,XPositions,1);
    Ref2 = S.ens_Data(:,XPositions,2);
    MRef = mean(Ref1+Ref2);
    StdRef = std(Ref1+Ref2);
    figure(8)
    hold on
    plot(XPositionsPlot,MRef,'Color',PlotColor,'LineWidth',3)
    %plot(XPositionsPlot,MRef+StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')
    %plot(XPositionsPlot,MRef-StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')

    %Synapses
    Ref1 = S.ens_Data(:,XPositions,3);
    MRef = mean(Ref1);
    StdRef = std(Ref1);
    figure(9)
    hold on
    plot(XPositionsPlot,MRef,'Color',PlotColor,'LineWidth',3)
    %plot(XPositionsPlot,MRef+StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')
    %plot(XPositionsPlot,MRef-StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')
%---------

    %All Filopodia


end

figure(10)
%xlabel('time (hours)')
%title(strcat('Number of filopodia (',AutophagyTypeName,')'))
title('Filopodia number (simulated)')
ylabel('Number/terminal')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'XTicklabel',cellstr(strcat('P',num2str((40:10:100)'))))
set(gca,'FontSize',14);
ylim([0 20])
print(10,'-dtiff',strcat(path,FigureName,'Filopodia.tiff'))


figure(8)
%xlabel('time (hours)')
title('Bulbous tip number (simulated)')
%title(strcat('Number of bulbous tips (',AutophagyTypeName,')'))
ylabel('Number/terminal')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'XTicklabel',cellstr(strcat('P',num2str((40:10:100)'))))
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 4])
print(8,'-dtiff',strcat(path,FigureName,'Bulb.tiff'))


figure(9)
%xlabel('time (hours)')
%title(strcat('Number of synapses (',AutophagyTypeName,')'))
title('Synapses number (simulated)')
ylabel('Number/terminal')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'XTicklabel',cellstr(strcat('P',num2str((40:10:100)'))))
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 30])
print(9,'-dtiff',strcat(path,FigureName,'Synapse.tiff'))


