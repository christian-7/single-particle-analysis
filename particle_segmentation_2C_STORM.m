%% Particle segmentation from 2C data (TS format)

% Input: 
% 
%         individual localization files   --> loc path and name
%         individual WF images            --> WF path and name


% Workflow: 

% Load data and WF images
% Adjust WF image contrast
% Binarize and segment
% Produce overlay images and save the extraced particles

%           Other options:

%                           DBSCAN filter
%                           overview plotting
%                           save images

% Output: 
% 
%       Variable Cent

%     Cent{i,1} = locs_Ch1(target_Ch1,1:9);
%     Cent{i,2} = length(target_Ch1);
%     Cent{i,3} = intI(i);
%     Cent{i,4} = Particles_WF{i,1};
%     Cent{i,5} = locs_Ch2(target_Ch2,1:9);
%     Cent{i,6} = length(target_Ch2);
%     Cent{i,7} = intI2(i);
%     Cent{i,8} = Particles_WF2{i,1};


%% Read Data
clear, clc, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pxl       = 107.99;                                                                 % Pixel size in nm
number    = 2;                                                                      % image number
filetype  = 2;                                                                      % 1 for TS, 2 for bstore
wfnumber  = 13;

%%%%%%%%%%%%%%%%% Manual Input %%%%%%%%%%%%%%%%%%%%%%%

WFpath          = ['Z:\Christian-Sieben\data_HTP\2016-04-01_humanCentriole_aTubNB_Sas6'];
WF_name         = ['FOV2_2.tif'];         

WFpath2         = ['Z:\Christian-Sieben\data_HTP\2016-04-01_humanCentriole_aTubNB_Sas6'];
WF_name2        = ['FOV2.tif'];      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Locpath1         = ['Z:\Christian-Sieben\data_HTP\2016-04-01_humanCentriole_aTubNB_Sas6\locResults_A647\humCent_aTubNB_Sas6_Pos' num2str(number) '_1'];
locName1         = ['humCent_aTubNB_Sas6_Pos' num2str(number) '_1_MMStack_locResults_DC.dat'];

Locpath2         = ['Z:\Christian-Sieben\data_HTP\2016-04-01_humanCentriole_aTubNB_Sas6\locResults_A750\humCent_aTubNB_Sas6_Pos' num2str(number) '_1'];
locName2         = ['humCent_aTubNB_Sas6_Pos' num2str(number) '_1_MMStack_locResults_processed.dat'];

savename         = ['humCent_aTubNB_Sas6_aTub' num2str(number) '_extractedParticles'];
savepath         =  'Z:\Christian-Sieben\data_HTP\2016-04-01_humanCentriole_aTubNB_Sas6\2C STORM analysis';

savepath_Images = 'Z:\Christian-Sieben\data_HTP\2016-04-01_humanCentriole_aTubNB_Sas6\2C STORM analysis';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n -- Path and File information --\n')
cd(Locpath1);

%% Load image

% GFP Channel

cd(WFpath);
ICh1   =   imread(WF_name);

% Storm Channel

cd(WFpath2);
ICh2  =  imread(WF_name2);

% Load localization data set

cd(Locpath1);
locs_Ch1=dlmread(locName1,',',1,0);

cd(Locpath2);
locs_Ch2=dlmread(locName2,',',1,0);

% Load header

cd(Locpath1);
file    = fopen(locName1);
line    = fgetl(file);
h       = regexp( line, ',', 'split' );

if filetype == 1;
    
x       = strmatch('"x [nm]"',h);
y       = strmatch('"y [nm]"',h);
LL      = strmatch('"loglikelihood"',h);

else

x       = strmatch('x [nm]',h);
y       = strmatch('y [nm]',h);
LL      = strmatch('loglikelihood',h);

end

fprintf('\n -- Data loaded --\n')

%% Adjust segmentation parameters

% Select Channel for Segmentation

Seg_Channel = 1;

if Seg_Channel==1;
    I=ICh1;
else I=ICh2;
end

close all

minWF = 70;
maxWF = 422;

% blur the image
figure('Position',[10 600 500 500],'name','Raw GFP Image'), imshow(I,[minWF maxWF],'InitialMagnification','fit');

% adjust the contrast of the raw image
I2 = imadjust(I,[0.001 0.01],[]);
figure('Position',[600 600 500 500],'name','Image after adjusted contrast'), imshow(I2,'InitialMagnification','fit');

% lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
G = fspecial('gaussian',[3 3],60); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imG = imfilter(I2,G,'same');
figure('Position',[1200 600 500 500],'name','Image after adjusted contrast, blurring'), imshow(imG,'InitialMagnification','fit');

% adjust the background
I3 = imadjust(imG,[0.3 0.5],[]);
figure('Position',[10 10 500 500],'name','Image after adjusted contrast, blurring, adjusted contrast'), imshow(I3,'InitialMagnification','fit');

% Make binary image
bin = im2bw(I3,0.3);
figure('Position',[600 10 500 500],'name','Binary image result'),imshow(bin,'InitialMagnification','fit')
[B,L,N,A] = bwboundaries(bin); % B - connectivity

%% Extract particles

close all

% Extract the integrated intensity of the Ch1 WF (A647) image for each ROI
% Extract the integrated intensity of the Ch2 WF (A750) image for each ROI

intI = [];
intI2 = [];

for i = 1:length(B);
    
    intI(i,1) = sum(sum(ICh1(min(B{i,1}(:,1)):max(B{i,1}(:,1)),min(B{i,1}(:,2)):max(B{i,1}(:,2)))));
    intI2(i,1) = sum(sum(ICh2(min(B{i,1}(:,1)):max(B{i,1}(:,1)),min(B{i,1}(:,2)):max(B{i,1}(:,2)))));
    
end

% Extract subimages from both WF channels

Particles_WF  ={};
Particles_WF2 ={};

box_size = 20;

for k=1:length(B)
    
%   the border particles WF ROIs are saved without the box
    
    
       if      min(B{k,1}(:,1))<box_size+1 | min(B{k,1}(:,2)) < box_size+1 | max(B{k,1}(:,1))+box_size>length(ICh1) | max(B{k,1}(:,2))+box_size>length(ICh1)
           
               Particles_WF{k,1}  = ICh1((min(B{k,1}(:,1))):(max(B{k,1}(:,1))),(min(B{k,1}(:,2))):(max(B{k,1}(:,2))));
               Particles_WF2{k,1} = ICh2((min(B{k,1}(:,1))):(max(B{k,1}(:,1))),(min(B{k,1}(:,2))):(max(B{k,1}(:,2))));
           
       else
               Particles_WF{k,1}  = ICh1((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));  
               Particles_WF2{k,1} = ICh2((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));   
            
       end
    
           
end

%Find the center of each particle and transform into an X,Y coordinate

Center=[];

for k=1:length(B)
        
            boundary    = B{k};
            Center(k,1) = (((max(B{k,1}(:,1))-min(B{k,1}(:,1)))/2)+min(B{k,1}(:,1)))*(pxl);             % Center of the segemented spot in nm
            Center(k,2) = (((max(B{k,1}(:,2))-min(B{k,1}(:,2)))/2)+min(B{k,1}(:,2)))*(pxl);             % Center of the segemented spot in nm
            Center(k,3) = max(pdist(B{k,1}))/1;                                                         % Size of the box, measure the max distance as input for the box size
            
end

CFX = (max(locs_Ch1(:,x)/pxl))./size(I);
CFY = (max(locs_Ch1(:,y)/pxl))./size(I);

center2(:,1)    =   Center(:,2)*CFX(:,1); % Center of the segmented spot in nm
center2(:,2)    =   Center(:,1)*CFY(:,1); % Center of the segmented spot in nm
center2(:,3)    =   Center(:,3);
Center          =   center2;

fprintf('\n -- Particles extracted --\n')

%% Build box around each Center and copy locs into separate variable -> structure Cent

Cent={}; 
    
for i=1:length(Center);
    
    vx1 = Center(i,1)+Center(i,3)*pxl;
    vx2 = Center(i,1)-Center(i,3)*pxl;
    vy1 = Center(i,2)+Center(i,3)*pxl;
    vy2 = Center(i,2)-Center(i,3)*pxl;
    
    target_Ch1 = find(locs_Ch1(:,x)>vx2 & locs_Ch1(:,x)<vx1 & locs_Ch1(:,y)>vy2 & locs_Ch1(:,y)<vy1);
    target_Ch2 = find(locs_Ch2(:,x)>vx2 & locs_Ch2(:,x)<vx1 & locs_Ch2(:,y)>vy2 & locs_Ch2(:,y)<vy1);
    
    Cent{i,1} = locs_Ch1(target_Ch1,1:9);
    Cent{i,2} = length(target_Ch1);
    Cent{i,3} = intI(i);
    Cent{i,4} = Particles_WF{i,1};
    Cent{i,5} = locs_Ch2(target_Ch2,1:9);
    Cent{i,6} = length(target_Ch2);
    Cent{i,7} = intI2(i);
    Cent{i,8} = Particles_WF2{i,1};
       
end 
    
fprintf('\n -- %f Particles selected from localitaion dataset --\n',length(Cent))

%% Overlay with Widefield image and plot correlations

% the coordinates need be backtransformed into pxl coordinates
% 1. divide by correction factor
% 2. divide trough pxl size

cd(savepath)

figure('Position',[150 150 400 400],'name',['Extracted Particles from Ch1 on WF of Ch1']);
imshow(ICh1,[minWF maxWF]); hold on;

for i = 1:length(Cent);
    
CentO=[];
CentO(:,1) = Cent{i,1}(:,x)/CFX(:,1);
CentO(:,2) = Cent{i,1}(:,y)/CFY(:,1);
    
scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');

hold on;
end

savefig('Overlay_extracted_particles_on_WF_Ch1.fig');

figure('Position',[150 150 400 400],'name',['Extracted Particles from Ch2 on WF of Ch2']);
imshow(ICh2,[minWF maxWF]); hold on;

for i = 1:length(Cent);
    
CentO=[];
CentO(:,1) = Cent{i,5}(:,x)/CFX(:,1);
CentO(:,2) = Cent{i,5}(:,y)/CFY(:,1);
    
scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');

hold on;
end

savefig('Overlay_extracted_particles_on_WF_Ch2.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% figure('Position',[1000 600 500 500],'name','All localizations');
% imshow(I,[minWF_GFP maxWF_GFP]); hold on;
% CentO = [];
% CentO(:,1) = locs_Ch1(:,x)/CFX(:,1);
% CentO(:,2) = locs_Ch1(:,y)/CFX(:,1);;
% scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');
% 
% savefig('Overlay_all_Locs_onGFP.fig');

% Plot for each Particle (1) the integrate intensity vs (2) the nbr of locs 

figure('Position',[400 100 900 300],'name','# of Locs vs. GFP Intensity');

subplot(1,3,1)
scatter(cell2mat(Cent(:,2)),cell2mat(Cent(:,3)),5,'filled','r');hold on;
xlabel('Nbr of locs Ch1');
ylabel('Ch1 WF int intensity');
box on;
axis square;

subplot(1,3,2);
scatter(cell2mat(Cent(:,6)),cell2mat(Cent(:,7)),5,'filled','r');hold on;
xlabel('Nbr of locs Ch2');
ylabel('Ch2 WF int intensity');
box on;
axis square;

subplot(1,3,3);
scatter(cell2mat(Cent(:,3)),cell2mat(Cent(:,7)),5,'filled','r');hold on;
xlabel('Ch1 intensity');
ylabel('Ch2 intensity');
box on;
axis square;


savefig('WF_correlation.fig');
save(savename,'Cent');


fprintf('\n -- Extracted centrioles saved --\n')


