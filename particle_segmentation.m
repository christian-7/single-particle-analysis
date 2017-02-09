%% Particle segmentation from single color data

% Load data and WF images
% Adjust WF image contrast
% Binarize and segment
% Produce overlay images and save the extraced particles

% Other options:

% DBSCAN filter
% overview plotting
% save images


%% Read Data
clear, clc, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pxl = 107.99;                                                                 % Pixel size in nm

%%%%%%%%%%%%%%%%% Manual Input %%%%%%%%%%%%%%%%%%%%%%%

WFpath      = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent\humanCent_aTub_NB_A647_WF13';
WF_name     = 'humanCent_aTub_NB_A647_WF13_MMStack_Pos0.ome.tif';         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

Locpath     = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent\locResults_aTub\humanCent_aTub_NB_A647_13';
locName     = 'humanCent_aTub_NB_A647_13_MMStack_Pos0_locResults_DC.dat';

savename =  'humanCent_aTub_NB_A647_Pos_13_extracted';
savepath =  'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent\analysed';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Data

% Load WF image
cd(WFpath);
I=imread(WF_name);

% Load data set
cd(Locpath);
locs=dlmread(locName,',',1,0);

% Load header
file = fopen(locName);
line = fgetl(file);
h = regexp( line, ',', 'split' );

x = strmatch('x [nm]',h);
y = strmatch('y [nm]',h);
LL = strmatch('loglikelihood',h);

fprintf('\n -- Data loaded --\n')

%% Adjust segmentation parameters
close all

minWF = 10;
maxWF = 200;

% #1: the original image
figure('Position',[10 600 500 500],'name','Raw Image'), imshow(I,[minWF maxWF],'InitialMagnification','fit');

% adjust the contrast of the raw image
I2 = imadjust(I,[0.1 0.5],[]);
figure('Position',[600 600 500 500],'name','Image after adjusted contrast'), imshow(I2,'InitialMagnification','fit');

G = fspecial('gaussian',[7 7],50); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imG = imfilter(I2,G,'same');
figure('Position',[1200 600 500 500],'name','Image after adjusted contrast, blurring'), imshow(imG,'InitialMagnification','fit');

% adjust the background
I3 = imadjust(imG,[0.01 0.3],[]);
figure('Position',[10 10 500 500],'name','Image after adjusted contrast, blurring, adjusted contrast'), imshow(I3,'InitialMagnification','fit');

% Make binary image
bin = im2bw(I3,0.3);
figure('Position',[600 10 500 500],'name','Binary image result'),imshow(bin,'InitialMagnification','fit')
[B,L,N,A] = bwboundaries(bin); % B - connectivity

%% Extract particles

close all

% Extract the integrated intensity of the GFP WF image for each ROI

intI = [];

for i = 1:length(B);
    
    intI(i,1) = sum(sum(I(min(B{i,1}(:,1)):max(B{i,1}(:,1)),min(B{i,1}(:,2)):max(B{i,1}(:,2)))));
    
end

% Extract subimages from GFP channel

Particles_WF={};
box_size = 1;

for k = 1:length(B)
    
    if      min(B{k,1}(:,1))<box_size+1; % if the x minimum overlaps with the border
        
            Particles_WF{k,1} = I((1):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));    
            
    elseif  min(B{k,1}(:,2))<box_size+1 % if the y minimum overlaps with the border
        
            Particles_WF{k,1} = I((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(1):(max(B{k,1}(:,2))+box_size));    
    
    elseif  max(B{k,1}(:,1))+box_size>length(I) | max(B{k,1}(:,2))+box_size>length(I)
        
            Particles_WF{k,1} = I((min(B{k,1}(:,1))-box_size):(length(I)),(min(B{k,1}(:,2))-box_size):(length(I)));
        
    else
        
            Particles_WF{k,1} = I((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));   
    end
           
end


%Find the center of each particle and transform into an X,Y coordinate

Center=[];

for k=1:length(B)
        
            boundary    = B{k};
            Center(k,1) = (((max(B{k,1}(:,1))-min(B{k,1}(:,1)))/2)+min(B{k,1}(:,1)))*(pxl);             % Center of the segemneted spot in nm
            Center(k,2) = (((max(B{k,1}(:,2))-min(B{k,1}(:,2)))/2)+min(B{k,1}(:,2)))*(pxl);             % Center of the segemneted spot in nm
            Center(k,3) = max(pdist(B{k,1}))/2;                                                         % Ssize of the box, measure the max distance as input for the box size
            
end

CFX = (max(locs(:,x)/pxl))./size(I);
CFY = (max(locs(:,y)/pxl))./size(I);

center2(:,1)    =   Center(:,2)*CFX(:,1); % Center of the segmented spot in nm
center2(:,2)    =   Center(:,1)*CFY(:,1); % Center of the segmented spot in nm
center2(:,3)    =   Center(:,3);
Center          =   center2;

fprintf('\n -- Particles extracted --\n')

%% Build box around each Center and copy locs into separate variable -> structure Cent

Cent={}; 
    
for i=1:length(Center);
    
    vx1=Center(i,1)+Center(i,3)*pxl;
    vx2=Center(i,1)-Center(i,3)*pxl;
    vy1=Center(i,2)+Center(i,3)*pxl;
    vy2=Center(i,2)-Center(i,3)*pxl;
    
    vx=find(locs(:,x)>vx2 & locs(:,x)<vx1);
    subset=locs(vx,1:9);

    vy=find(subset(:,y)>vy2 & subset(:,y)<vy1);

    Cent{i,1} = subset(vy,1:end);
    Cent{i,2} = intI(i);
    Cent{i,3} = Particles_WF{i,1};
       
end 
    
fprintf('\n -- %f Particles selected from localitaion dataset --\n',length(Cent))

%% Overlay with Widefield image

close all

% the coordinates need be backtransformed into pxl coordinates
% 1. divide by correction factor
% 2. divide trough pxl size

figure('Position',[100 400 600 600])
imshow(I,[minWF maxWF]); hold on;

for i=1:length(Cent);
    
CentO=[];
CentO(:,1) = Cent{i,1}(:,x)/CFX(:,1);
CentO(:,2) = Cent{i,1}(:,y)/CFY(:,1);
    
scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');

hold on;
end

cd(savepath);
savefig(['Overlay_extractedParticles_' savename '.fig']);


figure('Position',[700 100 300 300])
for i=1:length(Cent);
scatter(length(Cent{i}),Cent{i,2},5,'filled','r');hold on;
xlabel('Nbr of localizations');
ylabel('WF integrated intensity');
box on;
axis square;
end
% 
% cd(savepath);
% savefig(['Locs_vs_intInt_' savename '.fig']);

figure('Position',[200 100 300 300])
for i=1:length(Cent);
scatter(length(Cent{i}),Cent{i,2}/8400,5,'filled','r');hold on;
xlabel('Nbr of localizations');
ylabel('A647 molecules');
box on;
axis square;
end

cd(savepath);
savefig(['Locs_vs_molecules_' savename '.fig']);


% Save variable Particles

% 1. localization data
% 2. WF integrated intensity
% 3. WF Image subset

cd(savepath);
save(savename,'Cent');

fprintf('\n -- File Saved --\n')

%% DBSCAN particle size Filter

close all;

tic
fprintf('\n -- DBSCAN started --\n')

for m = 1:length(Cent);

dataDBS      = [];
dataDBS(:,1) = Cent{m,1}(:,1); % x in mum
dataDBS(:,2) = Cent{m,1}(:,2); % y in mum

if isempty(dataDBS)
else

% Run DBSCAN on each particle 

k   = 10;                                                 % minimum number of neighbors within Eps
Eps = 15;                                                 % minimum distance between points, nm

% fprintf('\n -- DBSCAN input and parameters selected --\n')

[class,type]=DBSCAN(dataDBS,k,Eps);     % uses parameters specified at input
class2=transpose(class);                % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);                  % (core: 1, border: 0, outlier: -1)

% fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

coreBorder = [];
coreBorder = find(type2 >= 0);

subset          = [];
subset          = Cent{m,1}(coreBorder,1:end);
subset(:,end+1) = class2(coreBorder);

% 
% subsetP = [];
% 
% subsetP(:,1)=dataDBS(coreBorder,1);
% subsetP(:,2)=dataDBS(coreBorder,2);
% subsetP(:,3)=class2(coreBorder);
% 
% figure('Position',[700 600 900 400])
% subplot(1,2,1)
% scatter(dataDBS(:,1),dataDBS(:,2),1);
% title('Raw Data')
% axis on
% axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
% box on
% 
% subplot(1,2,2)
% scatter(subsetP(:,1),subsetP(:,2),1,mod(subsetP(:,3),10))
% title('identified Clusters')
% axis on
% axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
% box on

% Select only large cluster(s)

for i=1:max(subset(:,end));
    
    vx = find(subset(:,end)==i);
    
    if length(vx)>100;
        
         Cent{m,4}{i,1} = subset(vx,1:end);
         
    end
end

end

end

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

% fprintf(' -- Particles DBSCAN filtered -- \n',toc)

%% Plot all particles in subplots

subplotSize = [];

for i=1:length(Cent);
    
 subplotSize(i,1) = length(Cent{i,4});
          
end

NbrSubplots = ceil(sqrt(sum(subplotSize)));
count = 1;
figure('Position',[100 100 1500 1000])

for i = 1:length(Cent);
    
    for j = 1:length(Cent{i,4});
    
            if isempty(Cent{i,4}{j,1})
            else
            
        subplot(NbrSubplots,NbrSubplots,count);
        scatter(Cent{i,4}{j,1}(:,1),Cent{i,4}{j,1}(:,2),1,'.b');
        title(['ROI :', num2str(i) ' Particle :', num2str(j)],'FontSize',7);
        count = count+1;
        box on
        axis off
        
            end             
    end
   
end

%% Plot each Isolated particles with the WF ROI

NbrSubplots = round(sqrt(length(Cent)));

figure('Position',[100 100 1000 1000])
count = 1;

for i=1:length(Cent);
    
    subplot(NbrSubplots+1,NbrSubplots,count);
    imshow(Cent{i,3});hold on;
    
for j = 1:length(Cent{i,4});
    
            if isempty(Cent{i,4}{j,1})
            else
                
CentO=[];
CentO(:,1) = Cent{i,1}(:,1)/CFX(:,1);
CentO(:,2) = Cent{i,1}(:,2)/CFY(:,1);

CentO(:,1) = CentO(:,1)-min(CentO(:,1))+box_size*pxl; % this normalizes to 0
CentO(:,2) = CentO(:,2)-min(CentO(:,2))+box_size*pxl;

scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');hold on; % Position over full WF
title(['ROI :', num2str(i) ' Particle :', num2str(j)],'FontSize',9);

            end
end
count = count+1;
end




%% Save variable Particles

% 1. localization data
% 2. WF integrated intensity
% 3. WF Image subset
% 4. DBSCAN filtered particles

savepath = 'Z:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\segmented_centrioles\extracted_DBSCAN_filtered'
savename =  'FOV10_Sas6_1000mW_10ms_A647_1_extractedParticles';

cd(savepath);
assignin('base',savename,Cent);
save(savename,savename);

close all
fprintf('\n -- File Saved --\n')



% end
