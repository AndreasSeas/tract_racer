% code  :: Perturbations.m
% descr :: Perform perturbations on JONES ROIs
% auth  :: Andreas Seas
% edits :: April 29, 2021 

%% clean slate
close all;clear;clc

%% set path
% this is the unique path to the various dsi_studio apps, each path
% delimited by a ":", make sure to change it to your own! Also note that on
% pc, the main folder delimiter is "\" v.s. "/" on mac and linux

setenv('PATH',"/Applications/dsi_studio.app/Contents/MacOS:/Users/andreasseas/Camino/bin:/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/Library/Apple/usr/bin:/Users/andreasseas/Developer/flutter/bin")

% echoes the path so you know you've actually changed it
!echo $PATH

%% get datetime
% for use to save the unique data files - record the parameters of each
% initial run at the time of run, that way you can go back and check 
dt=datestr(datetime('now'),"yymmdd_HHMMSS");

%% get unique transformation combinations for each ROI in 3x3x3 box
% these three strings for x, y, and z are all appended in the cmd line
% script to perturb the ROI

xcombos={"shiftnx", "noshift","shiftx"};
ycombos={"shiftny", "noshift","shifty"};
zcombos={"shiftnz", "noshift","shiftz"};
% %% prepare stuff for permutation
xyzXform=[];
addstr=[];
for x=1:3
    for y=1:3
        for z=1:3
            xyzXform=[xyzXform;[xcombos{x},ycombos{y},zcombos{z}]];
            temp="";
            for i=1:3
               if ~strcmp(xyzXform{end,i},"noshift")
                  temp=temp+","+ xyzXform{end,i};
               end
            end
            addstr=[addstr;temp];
        end
    end
end
% 

%% define ROI/As
% all ROIs were made in itk-SNAP with the original dwi image and parameters
% described in Jones, et al (2013)
% ... = name of roi; ... = niftiread of the ROI; ... = voxel number in ROI,
% to be used to calculate seeding density in order to initialize with 1000
% seeds per voxel
A="ROI_A_1vox.nii.gz";roiA=niftiread("ROI_A_1vox.nii.gz");...
    voxA=sum(sum(sum(roiA)));
B="ROI_B_1vox.nii.gz";roiB=niftiread("ROI_B_1vox.nii.gz");...
    voxB=sum(sum(sum(roiB)));
C="ROI_C_1vox.nii.gz";roiC=niftiread("ROI_C_1vox.nii.gz");...
    voxC=sum(sum(sum(roiC)));
D="ROI_D_1vox.nii.gz";roiD=niftiread("ROI_D_1vox.nii.gz");...
    voxD=sum(sum(sum(roiD)));
E="ROI_E_1vox.nii.gz";roiE=niftiread("ROI_E_1vox.nii.gz");...
    voxE=sum(sum(sum(roiE)));

%% define ROI/As, in this case 3 vox thick
% % this is just an import of the same voxel "sheets" except 3 voxels thick
A="ROI_A_3vox.nii.gz";roiA=niftiread("ROI_A_3vox.nii.gz");...
    voxA=sum(sum(sum(roiA)));
B="ROI_B_3vox.nii.gz";roiB=niftiread("ROI_B_3vox.nii.gz");...
    voxB=sum(sum(sum(roiB)));
C="ROI_C_3vox.nii.gz";roiC=niftiread("ROI_C_3vox.nii.gz");...
    voxC=sum(sum(sum(roiC)));
D="ROI_D_3vox.nii.gz";roiD=niftiread("ROI_D_3vox.nii.gz");...
    voxD=sum(sum(sum(roiD)));
E="ROI_E_3vox.nii.gz";roiE=niftiread("ROI_E_3vox.nii.gz");...
    voxE=sum(sum(sum(roiE)));


%% Subgenual Region, n=N per ROI loc comb

N=1;% number of iterations of tractography for each combination of ROIA + shift and ROIB + shift
numseed=(voxA+voxB)*1000;% calculation of the number of seeds per voxel
outputs=[];% data on ROI locations as well as output data for each ROI run
% designed with the intent of ANOVAn
% contents: ROIA transformation, ROIB transformation, volume of tracts,
%      number of tracts, statistical values from stat file from all
%      fibers, statistical values from stat file from fibers intersecting
%      ROI, ratio of number of fibers activated by stimulation
clc % clear the slate

fa=waitbar(0,"for ROIA");
fb=waitbar(0,"for ROIB");
fn=waitbar(0,"for combination");

savepts=[];
roipts=[];

for a=1:numel(addstr)
    waitbar(a/numel(addstr),fa);
    for b=1:numel(addstr)
        waitbar(b/numel(addstr),fb);
        for n=1:N
            waitbar(n/N,fn);
            disp("a = " + a + "/" + numel(addstr)+...
                ", b = " + b + "/" + numel(addstr)+...
                ", n = " + n + "/" + N);
            [temp1,temp2]=system("dsi_studio --action=trk"+...
                " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + A +addstr{a}+ ...
                " --roi2=" + B +addstr{b}+ ...
                " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+...
                " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track.tt.gz");
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=track.tt.gz --output=track.txt");
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=track.tt.gz --export=stat");
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=track.tt.gz "+...
                "--roi=elecmultisphereroi_rad10mm_subj1DW.nii.gz "+...
                "--output=trackstim.tt.gz");
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=trackstim.tt.gz --output=trackstim.txt");
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=trackstim.tt.gz --export=stat");
            % read tracks
            la=readcell("track.txt");
            
            [ro,co]=size(la);
            
            if co>1
                pts=cell(ro,1);
                allpts=[];
                la=readmatrix("track.txt");% read it in as matrix instead
                for i=1:ro
                    la2=rmmissing(la(i,:))';
                    pts{i}=[la2(1:3:end-2),la2(2:3:end-1),la2(3:3:end)];
                    allpts=[allpts;pts{i}];
                end
                
            else %if all comes in string
                pts=cell(size(la));
                allpts=[];

                for i=1:numel(la)
                    la2=strsplit(la{i}," ");
                    la2=cellfun(@str2num,la2)';
                    pts{i}=[la2(1:3:end-2),la2(2:3:end-1),la2(3:3:end)];
                    allpts=[allpts;pts{i}];
                end
            end
            
            shp=alphaShape(allpts(:,1),allpts(:,2),allpts(:,3));

            savepts=[savepts;{pts}];
            
            % read stimtracts
            la=readcell("trackstim.txt");
            
            [ro,co]=size(la);
            
            if co>1
                pts=cell(ro,1);
                allpts=[];
                la=readmatrix("trackstim.txt");% read it in as matrix instead
                for i=1:ro
                    la2=rmmissing(la(i,:))';
                    pts{i}=[la2(1:3:end-2),la2(2:3:end-1),la2(3:3:end)];
                    allpts=[allpts;pts{i}];
                end
                
            else %if all comes in string
                pts=cell(size(la));
                allpts=[];

                for i=1:numel(la)
                    la2=strsplit(la{i}," ");
                    la2=cellfun(@str2num,la2)';
                    pts{i}=[la2(1:3:end-2),la2(2:3:end-1),la2(3:3:end)];
                    allpts=[allpts;pts{i}];
                end
            end
            
            shp2=alphaShape(allpts(:,1),allpts(:,2),allpts(:,3));

            roipts=[roipts;{pts}];
            
            % load stats
            statvals=readcell('track.tt.gz.stat.txt');
            
            % get stats for intersect only
            
            intersectstat=readcell('trackstim.tt.gz.stat.txt');
            
            ratiointersect=intersectstat{1,2}/statvals{1,2};
            
            % save outputs
            outputs=[outputs;...
                [xyzXform(a,:),xyzXform(b,:),volume(shp),numel(la),...
                statvals(:,2)',intersectstat(:,2)'...
                ratiointersect]];
            
        end
        
    end
end

mkdir("results_"+dt);
cd("results_"+dt);
save outputvars.mat outputs savepts roipts
cd ..


%% individual track maker

numseed=(voxA+voxB)*1000;
[temp1,temp2]=system("dsi_studio --action=trk"+...
                " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + A +...
                " --roi2=" + B +...
                " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
                " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_AB.tt.gz");
            
numseed=(voxC+voxB)*1000;
[temp1,temp2]=system("dsi_studio --action=trk"+...
                " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + C +...
                " --roi2=" + B +...
                " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
                " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_BC.tt.gz");
            
numseed=(voxC+voxD)*1000;
[temp1,temp2]=system("dsi_studio --action=trk"+...
                " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + C +...
                " --roi2=" + D +...
                " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
                " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_CD.tt.gz");
            
numseed=(voxD+voxE)*1000;
[temp1,temp2]=system("dsi_studio --action=trk"+...
                " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + D +...
                " --roi2=" + E +...
                " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
                " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_DE.tt.gz");
            
[temp1,temp2]=system("dsi_studio --action=trk"+...
                " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + D +...
                " --roi2=" + E +" --roa=" + C + ...
                " --roa2=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
                " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_DEmC.tt.gz");