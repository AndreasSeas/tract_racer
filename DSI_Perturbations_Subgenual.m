% code  :: DSI_Perturbations_Subgenual.m
% descr :: Perform perturbations on JONES ROIs for subgenual aspect alone
% auth  :: Andreas Seas
% edits :: May 13, 2021

% relevant links:
%    - http://dsi-studio.labsolver.org/Manual/command-line-for-dsi-studio
%      (DSI studio command line interface manual)
%    - http://dsi-studio.labsolver.org/Manual/Fiber-Tracking
%      (DSI studio manual tracking guide, helpful for visualizing)

% file requirements within code folder
%    - ROI_*_#vox.nii.gz 
%        for all ROIs being used
%    - data.nii.gz.src.gz.gqi.1.25.fib.gz 
%        this is a modified version of the original dti that has tensor
%        data calculated... performed manually using dsi defaults
%    - mni1mm_rhem-csf_subj1DW.nii.gz
%        used as a region of avoidance here, already extant within patient
%        file. may need to regenerate to only include CSF and not exclude
%        rhem, but left this way for now to minimize run times
%    - elecmultisphereroi_rad10mm_subj1DW.nii.gz
%        roi showing region of activation thru stimulation... can be any
%        ROI outside of the ROIs used for gating
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
% individual runs
dt=datestr(datetime('now'),"yymmdd_HHMMSS");

%% get unique transformation combinations for each ROI in 3x3x3 box
% these three strings for x, y, and z are all appended in the cmd line 
% script to shift the ROI.

xcombos={"shiftnx", "noshift","shiftx"};% set of potential x shifts
ycombos={"shiftny", "noshift","shifty"};% set of potential x shifts
zcombos={"shiftnz", "noshift","shiftz"};% set of potential x shifts
% %% prepare stuff for permutation
xyzXform=[];% a summary of the different transforms... 27x3 cell
addstr=[];% the string to add to each ROI... 27x1 cell
% like this: ",shiftnx,shiftny,shiftnz"
for x=1:3% shift thru x combos
    for y=1:3% shift thru x combos
        for z=1:3% shift thru x combos
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
% A="ROI_A_3vox.nii.gz";roiA=niftiread("ROI_A_3vox.nii.gz");...
%     voxA=sum(sum(sum(roiA)));
% B="ROI_B_3vox.nii.gz";roiB=niftiread("ROI_B_3vox.nii.gz");...
%     voxB=sum(sum(sum(roiB)));
% C="ROI_C_3vox.nii.gz";roiC=niftiread("ROI_C_3vox.nii.gz");...
%     voxC=sum(sum(sum(roiC)));
% D="ROI_D_3vox.nii.gz";roiD=niftiread("ROI_D_3vox.nii.gz");...
%     voxD=sum(sum(sum(roiD)));
% E="ROI_E_3vox.nii.gz";roiE=niftiread("ROI_E_3vox.nii.gz");...
%     voxE=sum(sum(sum(roiE)));

%% Subgenual Region, n=N per ROI location combination
% run on the subgenual aspect!

N=6;% number of iterations of tractography for each combination of 
% ROIA + shift and ROIB + shift

numseed=(voxA+voxB)*1000;% calculation of the number of seeds per voxel

outputs=[];% data on ROI locations as well as output data for each ROI run
% designed with the intent of ANOVAn
% contents: ROIA transformation, ROIB transformation, volume of tracts,
%      number of tracts, statistical values from stat file from all
%      fibers, statistical values from stat file from fibers intersecting
%      ROI, ratio of number of fibers activated by stimulation
clc % clear the slate

% initialize waitbars to help track progress
fa=waitbar(0,"for ROIA");
fb=waitbar(0,"for ROIB");
fn=waitbar(0,"for combination (n = "+string(N)+")");

savepts=[];
% this initializes a cell of cells... let me break that down:
% savepts{i} = cell of cells,size is k x 1, with each cell element
%      containing the xyz matrix on the next line
% savepts{i}{j} = m x 3 matrix of xyz values for specific streamline within
%      specific run of dsi studio tractography

roipts=[];
% this initializes a similar cell of cells, this time to capture only the
% streamlines passing through a specific ROI (on top of ROIA and ROIB)...
% this is here used to find the streamlines also recruited by stimulation
% of SCC25

for a=1:numel(addstr)% go through each permutation of ROIA center (27 with 3x3x3 matrix)
    waitbar(a/numel(addstr),fa);% update wait bar
    for b=1:numel(addstr)% go through each permutation of ROIB center (27 with 3x3x3 matrix)
        waitbar(b/numel(addstr),fb);% update wait bar
        for n=1:N% go through above permutations N times (for ANOVAn)
            waitbar(n/N,fn);% update wait bar
            
            disp("a = " + a + "/" + numel(addstr)+...
                ", b = " + b + "/" + numel(addstr)+...
                ", n = " + n + "/" + N);% display the stage of this run
            
            % start sending system commands to dsi_studio... for reference,
            % the [temp1,temp2] at the beginning of each call is meant only
            % to decrease the verbosity of the system function, otherwise
            % you would see each of the dsi studio outputs for every one of
            % the (in this case) 27x27x6 runs
            
            % initial tractography run with ROIA+ROIB permutations
            [temp1,temp2]=system("dsi_studio --action=trk"+...
                " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + A +addstr{a}+ ...
                " --roi2=" + B +addstr{b}+ ...
                " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+...
                " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track.tt.gz");
            
            % convert track.tt.gz file to track.txt thru dsi... for import
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=track.tt.gz --output=track.txt");
            
            % get stats from dsi studio data
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=track.tt.gz --export=stat");
            
            
            % look for intersection of track.tt.gz dataset with the roi
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=track.tt.gz "+...
                "--roi=elecmultisphereroi_rad10mm_subj1DW.nii.gz "+...
                "--output=trackstim.tt.gz");
            %... and repeat similar method...
            
            % similarly export trackstim.tt.gz to trackstim.txt
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=trackstim.tt.gz --output=trackstim.txt");
            
            % similarly get stats of the trackstim data
            [temp1, temp2]=system("dsi_studio --action=ana "+...
                "--source=data.nii.gz.src.gz.gqi.1.25.fib.gz "+...
                "--tract=trackstim.tt.gz --export=stat");
            
            % read all tracts of this particular run into MATLAB
            la=readcell("track.txt");
            
            [ro,co]=size(la);% get the size of the import
            
            % this upcoming loop is there to avoid errors based on the
            % import... sometimes the import of the .txt comes in as a cell
            % of k x 1 size with each row representing a tract, but
            % sometimes there is a spacing error, or something that shifts
            % around a delimiter that confuses readcell, requiring a
            % readmatrix upload.
            
            if co>1% checks for case where more than 1 column
                % this indicates that it did not import correctly through
                % readcell and opts for the readmatrix approach
                pts=cell(ro,1);% rows still correspond to line # (\n is hard to miss)
                allpts=[];% initializes cell to hold all streamlines
                la=readmatrix("track.txt");% read it in as matrix instead
                for i=1:ro% go through each row
                    la2=rmmissing(la(i,:))';% removes any missing elements
                    pts{i}=[la2(1:3:end-2),la2(2:3:end-1),la2(3:3:end)];
                    % gets every 3rd element for x, and shifts for y and z
                    allpts=[allpts;pts{i}];% saves streamline i 
                end
                
            else %if all comes in k x 1 cell as expected
                pts=cell(size(la));% gets # lines
                allpts=[];% initializes cell to hold all lines

                for i=1:numel(la)% go thru each row
                    la2=strsplit(la{i}," ");% splits string to get #s
                    la2=cellfun(@str2num,la2)';% converts cell of #s to matrix
                    pts{i}=[la2(1:3:end-2),la2(2:3:end-1),la2(3:3:end)];
                    % gets every 3rd element for x, and shifts for y and z
                    allpts=[allpts;pts{i}];% saves streamline i
                end
            end
            
            shp=alphaShape(allpts(:,1),allpts(:,2),allpts(:,3));
            % computes a shape covering all streamline points

            savepts=[savepts;{pts}];% append to savepts cell of cells
            
            
            %%% below perform the same operations as above, just read data
            %%% from stimtracts, the tracts within the stimulated ROI
            % read stimtracts
            la=readcell("trackstim.txt");
            
            % 
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
            
            shp2=alphaShape(allpts(:,1),allpts(:,2),allpts(:,3));% save 
            % get second shape

            roipts=[roipts;{pts}];% save the points passing through ROI
            
            % load stats for all tracts
            statvals=readcell('track.tt.gz.stat.txt');
            
            % get stats for intersect only
            intersectstat=readcell('trackstim.tt.gz.stat.txt');
            
            % obtain the ratio of fibers that intersect to total number
            ratiointersect=intersectstat{1,2}/statvals{1,2};
            
            % save all outputs
            outputs=[outputs;...
                [xyzXform(a,:),xyzXform(b,:),volume(shp),numel(la),...
                statvals(:,2)',intersectstat(:,2)'...
                ratiointersect]];
            % all output values:
            %   - the x, y, and z transforms for that run for ROIA
            %   - the x, y, and z transforms for that run for ROIB
            %   - the volume of the shape encapsulating all fibers
            %   - the number of all fibers 
            %   - the stat values for all fibers (some redundancy with
            %     above)
            %   - the stat values for fibers intersecting ROI
            %   - # fibers intersecting ROI/# total fibers
            % note again there is some repetition of values, and also the
            % stat values contain more accurate measures of volume for the
            % entire bundle
            
        end
        
    end
end

mkdir("results_"+dt);% make a new directory with the datetime label from 
% the initiation of this code
cd("results_"+dt);% cd into that location
save outputvars.mat outputs savepts roipts %save results
cd .. % cd back to base directory


%% individual track maker (for visualization purposes)

% numseed=(voxA+voxB)*1000;
% [temp1,temp2]=system("dsi_studio --action=trk"+...
%                 " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + A +...
%                 " --roi2=" + B +...
%                 " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
%                 " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_AB.tt.gz");
%             
% numseed=(voxC+voxB)*1000;
% [temp1,temp2]=system("dsi_studio --action=trk"+...
%                 " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + C +...
%                 " --roi2=" + B +...
%                 " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
%                 " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_BC.tt.gz");
%             
% numseed=(voxC+voxD)*1000;
% [temp1,temp2]=system("dsi_studio --action=trk"+...
%                 " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + C +...
%                 " --roi2=" + D +...
%                 " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
%                 " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_CD.tt.gz");
%             
% numseed=(voxD+voxE)*1000;
% [temp1,temp2]=system("dsi_studio --action=trk"+...
%                 " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + D +...
%                 " --roi2=" + E +...
%                 " --roa=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
%                 " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_DE.tt.gz");
%             
% [temp1,temp2]=system("dsi_studio --action=trk"+...
%                 " --source=data.nii.gz.src.gz.gqi.1.25.fib.gz --roi=" + D +...
%                 " --roi2=" + E +" --roa=" + C + ...
%                 " --roa2=mni1mm_rhem-csf_subj1DW.nii.gz --seed_count="+numseed+" --max_length=120"+...
%                 " --fa_threshold=0.2 --turning_angle=60 --min_length=40 --output=track_DEmC.tt.gz");