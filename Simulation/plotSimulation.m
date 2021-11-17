function plotSimulation(time_points,ens_Data,Type,Scale)


TypeName = strrep(Type,'_','-');

Ref_yellow = [0.9 0.8 0.2];%[0.9100 0.4100 0.1700];


ens_sbul = ens_Data(:,:,1);
ens_Lbul = ens_Data(:,:,2);
ens_syn = ens_Data(:,:,3);


XPositions = (1:time_points)./(60*Scale)+40;

% %short lived Filopodia
% figure(4)
% hold on

% TMP_5_50_95 = prctileMat(ens_filo,[5 50 95 25 75],1);
% TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
% area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
% area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
% area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
% area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
%  
% plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
% line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
% xlabel('time (hours)')
% title(strcat('Short lived filopodia (',TypeName,')'))
% ylabel('Number per growth cone')
% set(get(gca,'xlabel'),'Fontsize',16);
% set(get(gca,'ylabel'),'Fontsize',16);
% set(get(gca,'title'),'Fontsize',16);
% set(gca,'FontSize',14);
% %set(gca,'FontWeight','b');
% ylim([0 20])

% %long-lived Filopodia
% figure(5)
% hold on
% TMP_5_50_95 = prctileMat(ens_sFil,[5 50 95 25 75],1);
% TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
% area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
% area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
% area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
% area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
%  
% plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
% line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
% xlabel('time (hours)')
% title(strcat('Long-lived filopodia (',TypeName,')'))
% ylabel('Number per Growthcone')
% set(get(gca,'xlabel'),'Fontsize',16);
% set(get(gca,'ylabel'),'Fontsize',16);
% set(get(gca,'title'),'Fontsize',16);
% set(gca,'FontSize',14);
% %set(gca,'FontWeight','b');
% ylim([0 25])


%% NEW Block regarding the Data ------


% %All Filopodia
% figure(10)
% hold on
% %XPositions = (1:time_points)./(60*Scale)+40;
% TMP_5_50_95 = prctileMat(ens_filo+ens_sFil,[5 50 95 25 75],1);
% TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
% area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
% area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
% area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
% area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
%  
% plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
% line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
% if ~strcmp(Type,'T25C')
%     S = load('./EnsembleData/DataSimple_T25C.mat');
%     Ref1 = S.ens_Data(:,:,1);
%     Ref2 = S.ens_Data(:,:,2);
%     MRef = mean(Ref1+Ref2);
%     plot(S.time_out,MRef,'Color',Ref_yellow,'LineWidth',3)
% end
% xlabel('time (hours)')
% title(strcat('Number of filopodia (',TypeName,')'))
% ylabel('Number per growth cone')
% set(get(gca,'xlabel'),'Fontsize',16);
% set(get(gca,'ylabel'),'Fontsize',16);
% set(get(gca,'title'),'Fontsize',16);
% set(gca,'FontSize',14);
% %set(gca,'FontWeight','b');
% 
% % %plot Data
% % errorbar(t./60+40,Data,StdData,'Linewidth',2)
% % plot(t./60+40, Data,'bo','MarkerSize',14,'MarkerFaceColor','b')
% 
% ylim([0 20])
% %print(10,'-dtiff',strcat('./Figures/SimulationFilopodia',Type,'.tiff'))

% Short-lived bulbs
figure(6)
hold on
TMP_5_50_95 = prctileMat(ens_sbul,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
plot(XPositions,mean(ens_sbul),'k:','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
if ~strcmp(Type,'T25C')
    S = load('./EnsembleData/DataSimple_T25C.mat');
    Ref = S.ens_Data(:,:,1);
    MRef = mean(Ref);
    plot(S.time_out,MRef,'Color',Ref_yellow,'LineWidth',3)
end

xlabel('time (hours)')
title(strcat('Number of short-lived bulbous tips (',TypeName,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 5])


% Long-lived bulbs
figure(7)
hold on
TMP_5_50_95 = prctileMat(ens_Lbul,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
plot(XPositions,mean(ens_Lbul),'k:','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
if ~strcmp(Type,'T25C')
    S = load('./EnsembleData/DataSimple_T25C.mat');
    Ref = S.ens_Data(:,:,2);
    MRef = mean(Ref);
    plot(S.time_out,MRef,'Color',Ref_yellow,'LineWidth',3)
end

xlabel('time (hours)')
title(strcat('Number of synaptogenic bulbous tips (',TypeName,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 5])

% All bulbs
figure(8)
hold on
TMP_5_50_95 = prctileMat(ens_Lbul+ens_sbul,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
plot(XPositions,mean((ens_sbul+ens_Lbul)),'k:','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
if ~strcmp(Type,'T25C')
    S = load('./EnsembleData/DataSimple_T25C.mat');
    Ref = S.ens_Data(:,:,2)+S.ens_Data(:,:,1);
    MRef = mean(Ref);
    plot(S.time_out,MRef,'Color',Ref_yellow,'LineWidth',3)
end

xlabel('time (hours)')
title(strcat('Number of bulbous tips (',TypeName,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 5])
%print(8,'-dtiff',strcat('./Figures/SimulationBulb',Type,'.tiff'))


%Synapses
figure(9)
hold on
TMP_5_50_95 = prctileMat(ens_syn,[5 50 95 25 75],1);
TMP_5_50_95(isnan(TMP_5_50_95)) = 0;
area(XPositions,TMP_5_50_95(3,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(5,:),-1.1,'EdgeColor','none','FaceColor',[.6 .6 .6]);%-> Any mutated
area(XPositions,TMP_5_50_95(4,:),-1.1,'EdgeColor','none','FaceColor',[.8 .8 .8]);%-> Any mutated
area(XPositions,TMP_5_50_95(1,:),-1.1,'EdgeColor','none','FaceColor',[1 1 1]);%-> Any mutated
 
plot(XPositions,TMP_5_50_95(2,:),'r','LineWidth',3)
line([XPositions(1) XPositions(end)],[0 0],'Color','k','LineStyle','--')
plot(XPositions,mean(ens_syn),'k:','LineWidth',3)
if ~strcmp(Type,'T25C')
    S = load('./EnsembleData/DataSimple_T25C.mat');
    Ref = S.ens_Data(:,:,3);
    MRef = mean(Ref);
    plot(S.time_out,MRef,'Color',Ref_yellow,'LineWidth',3)
end

xlabel('time (hours)')
title(strcat('Number of synapses (',TypeName,')'))
ylabel('Number per growth cone')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 50])
%print(9,'-dtiff',strcat('./Figures/SimulationSynapse',Type,'.tiff'))


end

