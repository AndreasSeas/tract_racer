% TO MAKE HISTOGRAMS
close all;clear;clc
%% get and prepare the data
cd 'results_210419_213352(subgenual, 1vox ROIs, n6)'/
load outputvars.mat
outputs1=outputs; roipts1=roipts;savepts1=savepts;
cd ..
cd 'results_210420_043946(subgenual, 6vox ROIs, n1)'/
load outputvars.mat
outputs3=outputs; roipts3=roipts;savepts3=savepts;
cd ..
alldata1=cellfun(@str2num,outputs1(:,7:end));alldata3=cellfun(@str2num,outputs3(:,7:end));

% %% do ANOVA
% % number
% idx=3;
% [p,tbl,stats] = anovan(alldata1(:,idx),{outputs1(:,1),outputs1(:,2),outputs1(:,3),outputs1(:,4),outputs1(:,5),outputs1(:,6)});
% [p,tbl,stats] = anovan(alldata3(:,idx),{outputs3(:,1),outputs3(:,2),outputs3(:,3),outputs3(:,4),outputs3(:,5),outputs3(:,6)});
% 
% % vol
% idx=9;
% [p,tbl,stats] = anovan(alldata1(:,idx),{outputs1(:,1),outputs1(:,2),outputs1(:,3),outputs1(:,4),outputs1(:,5),outputs1(:,6)});
% [p,tbl,stats] = anovan(alldata3(:,idx),{outputs3(:,1),outputs3(:,2),outputs3(:,3),outputs3(:,4),outputs3(:,5),outputs3(:,6)});
% 
% % trunkvol
% idx=10;
% [p,tbl,stats] = anovan(alldata1(:,idx),{outputs1(:,1),outputs1(:,2),outputs1(:,3),outputs1(:,4),outputs1(:,5),outputs1(:,6)});
% [p,tbl,stats] = anovan(alldata3(:,idx),{outputs3(:,1),outputs3(:,2),outputs3(:,3),outputs3(:,4),outputs3(:,5),outputs3(:,6)});
% 
% % recruitment
% idx=numel(alldata1(1,:));
% [p,tbl,stats] = anovan(alldata1(:,idx),{outputs1(:,1),outputs1(:,2),outputs1(:,3),outputs1(:,4),outputs1(:,5),outputs1(:,6)});
% [p,tbl,stats] = anovan(alldata3(:,idx),{outputs3(:,1),outputs3(:,2),outputs3(:,3),outputs3(:,4),outputs3(:,5),outputs3(:,6)});


%% number of tracts
maketheplot(3, "fig_tracts_squat", "# tracts",alldata1,outputs1);


%% volume
maketheplot(9, "fig_volume_squat", "volume",alldata1,outputs1);

%% trunk volume
maketheplot(10, "fig_trunkvolume", "trunk volume",alldata1,outputs1);

%% recruitment
maketheplot(65, "fig_recruitment_squat", "recruitment",alldata1,outputs1);


%% function build 
function maketheplot(idx, filename, xlabelname,alldata1,outputs1)

xyz={'x','y','z','x','y','z'};
braincoord={'right','left';'anterior','posterior';'inferior','superior'};
braincoord=[braincoord;braincoord];
roii={'A','A','A','B','B','B'};

figure
set(gcf,'visible','off');
h=histogram(alldata1(:,idx));
binedges=h.BinEdges;
xvals=binedges(2:end)-0.5*h.BinWidth;
y0=h.Values;

for i=1:6
    g1=outputs1(:,i)=="shiftn"+xyz{i};
    g2=outputs1(:,i)=="noshift";
    g3=outputs1(:,i)=="shift"+xyz{i};
    h=histogram(alldata1(g1,idx),binedges);
    y1=h.Values;
    h=histogram(alldata1(g2,idx),binedges);
    y2=h.Values;
    h=histogram(alldata1(g3,idx),binedges);
    y3=h.Values;

    f=figure;
    set(gcf,'visible','off');
    hold on;
    b=bar(xvals,y1);
    b.FaceAlpha=0.4;
    b=bar(xvals,y2);
    b.FaceAlpha=0.4;
    b=bar(xvals,y3);
    b.FaceAlpha=0.4;b.FaceAlpha=0.4;
    bar(xvals,y0,"FaceColor","none");
    xlabel(xlabelname);
    ylabel("frequency");
    
    legend("ROI"+roii{i}+" "+ braincoord{i,1},...
        "ROI"+roii{i}+" neutral",...
        "ROI"+roii{i}+" "+ braincoord{i,2},"total","location","best")
    grid on
    set(gcf,'position',[10,10,600,400]);
    set(gca,'FontSize',20)
    set(gca, 'YScale', 'log')
    hold off;
    saveas(f,filename + "_" + roii{i} + "_" + xyz{i},"jpg");
end

end
% 
% figure
% set(gcf,'visible','off');
% h=histogram(alldata1(:,3));
% binedges=h.BinEdges;
% xvals=binedges(2:end)-0.5*h.BinWidth;
% y0=h.Values;
% 
% 
% % segment by ROIA_y
% idx=3;
% g1=outputs1(:,2)=="shiftny";
% g2=outputs1(:,2)=="noshift";
% g3=outputs1(:,2)=="shifty";
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("number of tracts");
% ylabel("frequency");
% legend("ROIA posterior","ROIA neutral","ROIA anterior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"number_A_y","jpg");
% 
% 
% % segment by ROIA_z
% figure
% set(gcf,'visible','off');
% idx=3;
% g1=outputs1(:,3)=="shiftnz";
% g2=outputs1(:,3)=="noshift";
% g3=outputs1(:,3)=="shiftz";
% 
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("number of tracts");
% ylabel("frequency");
% 
% legend("ROIA inferior","ROIA neutral","ROIA superior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"number_A_z","jpg");
% 
% figure
% set(gcf,'visible','off');
% h=histogram(alldata1(:,3));
% binedges=h.BinEdges;
% xvals=binedges(2:end)-0.5*h.BinWidth;
% y0=h.Values;
% 
% % segment by ROIB_y
% idx=3;
% g1=outputs1(:,5)=="shiftny";
% g2=outputs1(:,5)=="noshift";
% g3=outputs1(:,5)=="shifty";
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("number of tracts");
% ylabel("frequency");
% legend("ROIB posterior","ROIB neutral","ROIB anterior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"number_B_y","jpg");
% 
% % TO MAKE HISTOGRAMS
% close all;clc
% %% volume
% idx=9;
% figure
% set(gcf,'visible','off');
% h=histogram(alldata1(:,idx));
% binedges=h.BinEdges;
% xvals=binedges(2:end)-0.5*h.BinWidth;
% y0=h.Values;
% 
% % segment by ROIA_y
% 
% g1=outputs1(:,2)=="shiftny";
% g2=outputs1(:,2)=="noshift";
% g3=outputs1(:,2)=="shifty";
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("region volume");
% ylabel("frequency");
% legend("ROIA posterior","ROIA neutral","ROIA anterior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"volume_A_y","jpg");
% 
% 
% % segment by ROIA_z
% figure
% set(gcf,'visible','off');
% g1=outputs1(:,3)=="shiftnz";
% g2=outputs1(:,3)=="noshift";
% g3=outputs1(:,3)=="shiftz";
% 
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("region volume");
% ylabel("frequency");
% 
% legend("ROIA inferior","ROIA neutral","ROIA superior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"volume_A_z","jpg");
% 
% 
% % segment by ROIB_y
% g1=outputs1(:,5)=="shiftny";
% g2=outputs1(:,5)=="noshift";
% g3=outputs1(:,5)=="shifty";
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("region volume");
% ylabel("frequency");
% legend("ROIB posterior","ROIB neutral","ROIB anterior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"volume_B_y","jpg");
% 
% %% trunk volume
% idx=10;
% figure
% set(gcf,'visible','off');
% h=histogram(alldata1(:,idx));
% binedges=h.BinEdges;
% xvals=binedges(2:end)-0.5*h.BinWidth;
% y0=h.Values;
% 
% % segment by ROIA_y
% 
% g1=outputs1(:,2)=="shiftny";
% g2=outputs1(:,2)=="noshift";
% g3=outputs1(:,2)=="shifty";
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("trunk volume");
% ylabel("frequency");
% legend("ROIA posterior","ROIA neutral","ROIA anterior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"trunkvolume_A_y","jpg");
% 
% 
% % segment by ROIA_z
% figure
% set(gcf,'visible','off');
% g1=outputs1(:,3)=="shiftnz";
% g2=outputs1(:,3)=="noshift";
% g3=outputs1(:,3)=="shiftz";
% 
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("trunk volume");
% ylabel("frequency");
% 
% legend("ROIA inferior","ROIA neutral","ROIA superior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"trunkvolume_A_z","jpg");
% 
% 
% % segment by ROIB_y
% g1=outputs1(:,5)=="shiftny";
% g2=outputs1(:,5)=="noshift";
% g3=outputs1(:,5)=="shifty";
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("trunk volume");
% ylabel("frequency");
% legend("ROIB posterior","ROIB neutral","ROIB anterior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"trunkvolume_B_y","jpg");
% 
% 
% %% trunk volume
% idx=65;
% figure
% set(gcf,'visible','off');
% h=histogram(alldata1(:,idx));
% binedges=h.BinEdges;
% xvals=binedges(2:end)-0.5*h.BinWidth;
% y0=h.Values;
% 
% % segment by ROIA_y
% 
% g1=outputs1(:,2)=="shiftny";
% g2=outputs1(:,2)=="noshift";
% g3=outputs1(:,2)=="shifty";
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("recruitment");
% ylabel("frequency");
% legend("ROIA posterior","ROIA neutral","ROIA anterior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"recruitment_A_y","jpg");
% 
% 
% % segment by ROIA_z
% figure
% set(gcf,'visible','off');
% g1=outputs1(:,3)=="shiftnz";
% g2=outputs1(:,3)=="noshift";
% g3=outputs1(:,3)=="shiftz";
% 
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("recruitment");
% ylabel("frequency");
% 
% legend("ROIA inferior","ROIA neutral","ROIA superior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"recruitment_A_z","jpg");
% 
% 
% % segment by ROIB_y
% g1=outputs1(:,5)=="shiftny";
% g2=outputs1(:,5)=="noshift";
% g3=outputs1(:,5)=="shifty";
% h=histogram(alldata1(g1,idx),binedges);
% y1=h.Values;
% h=histogram(alldata1(g2,idx),binedges);
% y2=h.Values;
% h=histogram(alldata1(g3,idx),binedges);
% y3=h.Values;
% 
% f=figure;
% set(gcf,'visible','off');
% hold on;
% b=bar(xvals,y1);
% b.FaceAlpha=0.4;
% b=bar(xvals,y2);
% b.FaceAlpha=0.4;
% b=bar(xvals,y3);
% b.FaceAlpha=0.4;b.FaceAlpha=0.4;
% bar(xvals,y0,"FaceColor","none");
% xlabel("recruitment");
% ylabel("frequency");
% legend("ROIB posterior","ROIB neutral","ROIB anterior","sum","location","best")
% grid on
% set(gcf,'position',[10,10,600,800]);
% set(gca,'FontSize',20)
% hold off;
% saveas(f,"recruitment_B_y","jpg");