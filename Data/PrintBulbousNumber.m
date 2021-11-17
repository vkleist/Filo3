clear all
close all

load('AllData.mat');

Mutants = {'T18C';'T25C'; 'T29C'};
%NrBootstrap = 1000;

edges = (0:10:60);
Labels = cellstr(num2str(edges'));
Labels{end} = strcat('>',Labels{end});
%Colors = ['b';'r';'k';'m';'y';'c';'g'];
Colors = [.5 .5 .5;0 0 1;1 0 0;1 0 1];
AllBulbTips = [];
Groups = {};
EventVars = {};
counter = 1;
for i = 1:length(Mutants)
    mutant = char(Mutants(i))
    BulbTips = Output.(mutant).P60.Bulb.LTimes;
    GrowthCones = Output.(mutant).P60.Bulb.GC;
    Censored = Output.(mutant).P60.Bulb.Censored;
    AllBulbTips = [AllBulbTips;BulbTips];
    N = length(BulbTips);
    Groups(counter:counter+N-1) = {strcat(mutant(2:3),'°C')};
    EventVars(counter:counter+N-1) = {'Event'};
    idx = find(Censored);
    EventVars(idx+counter-1) = {'Censored'};
    counter = counter + N;
        
    figure(i)
    [N,edges] = histcounts(BulbTips,edges);%
    heights = N./sum(N);
    xplot = edges(1:end-1)+0.5*(edges(2)-edges(1));
    bar(xplot,heights,1,'FaceColor',Colors(i,:))
    %[N,edges] = histogram(BulbTips,6,'Normalization','probability');
    MutantName = strrep(mutant,'_','-');
    title(MutantName)
    xlabel('Lifetime (min)')
    ylabel('Frequency')
    fs = 18;
    set(gca,'FontSize',fs);
    set(get(gca,'title'),'Fontsize',fs);  
    set(get(gca,'xlabel'),'Fontsize',fs); 
    set(get(gca,'ylabel'),'Fontsize',fs);   
    set(gca,'XTick',edges,'XTickLabels',Labels)
    xlim([0 60])

    NrGC = length(unique(GrowthCones));
    NrInCategory = nan(length(N),NrGC);
    NrBulbsPerGC = nan(1,NrGC);
    NrDisappeared = zeros(1,NrGC);
    NrEmerged = zeros(1,NrGC);
    NrCensoredInCategory = nan(length(N),NrGC);
    for j = 1:NrGC
        idx = find(GrowthCones==j);
        NrBulbsPerGC = length(idx);
        
        %measurement of short lived bulbs, assuming all disappear
        %eventually
        BulbsInGC = BulbTips(idx);
        CensoredInGC = Censored(idx);
        lgx = BulbsInGC<40;
        NrDisappeared(j) = sum(lgx);
        % measurement of all bulb tips that eventually disappeared
        lgx = BulbsInGC<60;
        NrEmerged(j) = sum(lgx);
         
        for z = 1:length(N)-1
            IdxGC = BulbsInGC>=edges(z) & BulbsInGC<edges(z+1);
            NrInCategory(z,j) = sum(IdxGC);
            NrCensoredInCategory(z,j) = sum(CensoredInGC(IdxGC));
        end
        %Interval right-open
        z = length(N);
        IdxGC = BulbsInGC>=edges(z) & BulbsInGC<=edges(z+1);
        NrInCategory(z,j) = sum(IdxGC);
        NrCensoredInCategory(z,j) = sum(CensoredInGC(IdxGC));
        
        txt = strcat(num2str(NrInCategory(:,j)),'/',num2str(NrBulbsPerGC));
        txt2 = strcat('(',num2str(NrCensoredInCategory(:,j)),'*)');
        cs = [.4 .4 .4];
        if i == 1
            cs = [1 0 0];
        end
%         for z = 1:length(N)
%             %text(edges(z)+(j-1)*(edges(2)-edges(1))/NrGC,heights(z),txt(z,:),'Color',Colors(j),'VerticalAlignment','bottom')
%             text(xplot(z),0.85-0.05*(j-1) ,txt(z,:),'Color',cs,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',10)
%             if NrCensoredInCategory(z,j) > 0
%                 text(xplot(z),0.85-0.05*(j-1)-0.025 ,txt2(z,:),'Color',cs,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',6)
%         
%             end
%         end
    end
    ylim([0 0.9])
    Path = './Figures/';
    try cd (Path)
        cd ..
    catch
        mkdir(Path)
    end
    Name = strcat('Fig_BulbousLifeTime',mutant);

    %print(i-1,'-dpdf',strcat(Path,Name,'.pdf'))
    print(i,'-dtiff',strcat(Path,Name,'.tiff'))
    print(i,'-depsc',strcat(Path,Name,'.eps'))
    
    %lgx = isoutlier(BulbTips);
    lgx = BulbTips>=40;
    if mean(BulbTips(lgx)) > mean(BulbTips(~lgx))
        Long = lgx;
        Short = ~lgx;
    else
        Long = ~lgx;
        Short = lgx;
    end
    
    ShortLived = BulbTips(Short);  
    
end
[p, fh, stats,data] = MatSurv(AllBulbTips, EventVars, Groups,'Xlabel','minutes','YLabel','Bulb survival probability','CensorInRT',true);%,'PairWiseP',true);%,'Print',true);%,'RT_KMplot',true)
for i = 1:length(Mutants)
    Output.(char(Mutants{i})).P60.Bulb.SurvProb = data.GROUPS(i).KM_ALL;
end
save('./AllData.mat','Output');