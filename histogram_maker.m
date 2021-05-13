% code  :: histogram_maker.m
% descr :: make histograms of entire dataset
% auth  :: Andreas Seas
% edits :: May 13, 2021

%% clean slate
close all;clear;clc

%% decide what to save
savefigs=1;
savetables=0;

%% get and prepare the data
foldername='sampleresults_210419_213352(subgenual, 1vox ROIs, n6)';
cd(foldername)
load outputvars.mat% load data
alldata=cellfun(@str2num,outputs(:,7:end));

%% do ANOVA for dataset
% number
idx=3; %the column number
[p,tbl_num,stats] = anovan(alldata(:,idx),{outputs(:,1),outputs(:,2),outputs(:,3),outputs(:,4),outputs(:,5),outputs(:,6)});

% vol
idx=9;
[p,tbl_vol,stats] = anovan(alldata(:,idx),{outputs(:,1),outputs(:,2),outputs(:,3),outputs(:,4),outputs(:,5),outputs(:,6)});

% trunkvol
idx=10;
[p,tbl_trunkvol,stats] = anovan(alldata(:,idx),{outputs(:,1),outputs(:,2),outputs(:,3),outputs(:,4),outputs(:,5),outputs(:,6)});

% recruitment
idx=numel(alldata(1,:));
[p,tbl_rec,stats] = anovan(alldata(:,idx),{outputs(:,1),outputs(:,2),outputs(:,3),outputs(:,4),outputs(:,5),outputs(:,6)});

if savetables==1
    writecell(tbl_num,"anova_number.txt")
    writecell(tbl_vol,"anova_volume.txt")
    writecell(tbl_trunkvol,"anova_trunkvolume.txt")
    writecell(tbl_rec,"anova_recruitment.txt")
end

%% make histograms for number of tracts
maketheplot(3, "fig_tracts_squat", "# tracts",alldata,outputs,savefigs);

%% make histograms for volume
maketheplot(9, "fig_volume_squat", "volume",alldata,outputs,savefigs);

%% make histograms for trunk volume
maketheplot(10, "fig_trunkvolume", "trunk volume",alldata,outputs,savefigs);

%% make histograms for recruitment
maketheplot(65, "fig_recruitment_squat", "recruitment",alldata,outputs,savefigs);

cd .. % go back to OG folder

%% function build 
function maketheplot(idx, filename, xlabelname,alldata1,outputs1,savemeorno)
% get set of xyz values
xyz={'x','y','z','x','y','z'};

% get -x, x, -y, y, -z, z values
braincoord={'right','left';'anterior','posterior';'inferior','superior'};
braincoord=[braincoord;braincoord];

% corresponding roi values... in this case ROIA and ROIB... all is function
% of how you did loops for the run
roii={'A','A','A','B','B','B'};

figure
set(gcf,'visible','off');% set figure to not display for purposes of saving memory
h=histogram(alldata1(:,idx));% make baseline histogram to get optimal bin edges
binedges=h.BinEdges;% get bin edges
xvals=binedges(2:end)-0.5*h.BinWidth;% get bin centers
y0=h.Values;% get values

for i=1:6% go through ROIA and ROIB histograms... first 3 are ROIA, next are ROIB
    
    %segment out the data based on if has shiftnx, shifty, noshift, etc...
    g1=outputs1(:,i)=="shiftn"+xyz{i};
    g2=outputs1(:,i)=="noshift";
    g3=outputs1(:,i)=="shift"+xyz{i};
    
    % get base histograms for each subset in the x, y, or z dimension
    h=histogram(alldata1(g1,idx),binedges);
    y1=h.Values;
    h=histogram(alldata1(g2,idx),binedges);
    y2=h.Values;
    h=histogram(alldata1(g3,idx),binedges);
    y3=h.Values;

    % make final figure
    f=figure;
    set(gcf,'visible','off');% turns figure off again
    hold on;
    b=bar(xvals,y1);% makes barplot for shiftn
    b.FaceAlpha=0.4;% sets face to 0.4
    b=bar(xvals,y2);% makes barplot for noshift
    b.FaceAlpha=0.4;
    b=bar(xvals,y3);% makes barplot for shift
    b.FaceAlpha=0.4;b.FaceAlpha=0.4;
    bar(xvals,y0,"FaceColor","none");% barplot for all data (==y1+y2+y3)
    xlabel(xlabelname);%dynamic label name
    ylabel("frequency");
    % make legend 
    legend("ROI"+roii{i}+" "+ braincoord{i,1},...
        "ROI"+roii{i}+" neutral",...
        "ROI"+roii{i}+" "+ braincoord{i,2},"total","location","best")
    
    grid on
    
    set(gca,'FontSize',20)
    
    %%% below parameters can be modulated through an if statement based on
    %%% what you are graphing (recruitment or other) for better vis
    set(gcf,'position',[10,10,600,800]);% set fig size
    % set(gca, 'YScale', 'log') % set figure yscale to log if needed
    hold off;
    if savemeorno==1
        saveas(f,filename + "_" + roii{i} + "_" + xyz{i},"jpg");
    end
end

end

