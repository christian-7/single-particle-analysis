%% Particle segmentation from single color data in H5 format

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

pxl = 106; % 107.99                                                                 % Pixel size in nm

savepath =  'Z:\Christian-Sieben\data_HTP\2017-01-17_humanCent_Cep152-aTub\analysis';

H5Folder =  'Z:\Christian-Sieben\data_HTP\2017-01-17_humanCent_Cep152-aTub';

filename =  '2016-01-19_humanCent_aTub-Cep152_A647';

grouppath = '//humanCent_aTub-Cep152-ML/humanCent_aTub-Cep152-ML_'; % ['/' filename '/' filename '_']

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1;
    
cd(H5Folder); 
fileinfo = h5info([filename '.h5']);

% Load WF image

I = h5read([filename '.h5'],[[grouppath num2str(i)] '/WidefieldImage_ChannelA647_Pos0/image_data']);

I = (fliplr(imrotate(I,-90)));

% Load data set and generate locs variable

dataset = h5read([filename '.h5'],[[grouppath num2str(i)] '/Localizations_ChannelA647_Pos0/table']);

xCol = 1; yCol = 2; frameCol = 3; uncCol = 4; photonCol = 5; LLCol = 6; sigmaCol = 7;

locs(:,xCol)        = dataset.x0x5Bnm0x5D;
locs(:,yCol)        = dataset.y0x5Bnm0x5D;
locs(:,frameCol)    = dataset.frame;
locs(:,uncCol)      = dataset.uncertainty0x5Bnm0x5D;
locs(:,photonCol)   = dataset.intensity0x5Bphoton0x5D;
locs(:,LLCol)       = dataset.loglikelihood;
locs(:,sigmaCol)    = dataset.sigma0x5Bnm0x5D;

savename =  [filename '-ML_extractedCent_' num2str(i)];

end

fprintf('\n -- Data loaded --\n')

%% Find segmentation parameters

close all

% Miji;
MIJ.createImage('result', I, true);
cd(H5Folder); 

%% Adjust segmentation parameters

close all
MIJ.run('Close All');

minWF = 100;
maxWF = 19000;

% #1: the original image
figure('Position',[10 600 500 500],'name','Raw Image'), imshow(I,[minWF maxWF],'InitialMagnification','fit');

% #2: adjust the contrast of the raw image
I2 = imadjust(I,[0.04 0.08],[]);
figure('Position',[600 600 500 500],'name','Image after adjusted contrast'), imshow(I2,'InitialMagnification','fit');

G = fspecial('gaussian',[7 7],50); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imG = imfilter(I2,G,'same');
figure('Position',[1200 600 500 500],'name','Image after adjusted contrast, blurring'), imshow(imG,'InitialMagnification','fit');

% adjust the background
I3 = imadjust(imG,[0.01 0.3],[]);
figure('Position',[10 10 500 500],'name','Image after adjusted contrast, blurring, adjusted contrast'), imshow(I3,'InitialMagnification','fit');

% Make binary image
bin = im2bw(I3,0.2);
figure('Position',[600 10 500 500],'name','Binary image result'),imshow(bin,'InitialMagnification','fit')
[B,L,N,A] = bwboundaries(bin); % B - connectivity

% Otsu automatic Thresholding

% level = graythresh(I);
% bin = im2bw(I,level);
% figure('Position',[600 600 500 500],'name','Image after adjusted contrast'), imshow(bin,'InitialMagnification','fit');
% [B,L,N,A] = bwboundaries(bin); % B - connectivity

%% Extract particles

close all

% Extract the integrated intensity of the WF image for each ROI

intI = [];

for i = 1:length(B);
    
    intI(i,1) = sum(sum(I(min(B{i,1}(:,1)):max(B{i,1}(:,1)),min(B{i,1}(:,2)):max(B{i,1}(:,2)))));
    
end

% Extract subimages from GFP channel
% Exlcude border particles

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
            Center(k,3) = max(pdist(B{k,1}))/1;                                                         % Size of the box, measure the max distance as input for the box size
            
end

CFX = (max(locs(:,xCol)/pxl))./size(I);
CFY = (max(locs(:,yCol)/pxl))./size(I);

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
    
    vx = find(locs(:,xCol)>vx2 & locs(:,xCol)<vx1 & locs(:,yCol)>vy2 & locs(:,yCol)<vy1);
    
    subset = locs(vx,1:end);

    Cent{i,1} = subset;
    Cent{i,2} = length(vx);
    Cent{i,3} = intI(i);
    Cent{i,4} = Particles_WF{i,1};
           
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
CentO(:,1) = Cent{i,1}(:,1)/CFX(:,1);
CentO(:,2) = Cent{i,1}(:,2)/CFY(:,1);
    
scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');

hold on;
end

cd(savepath);
savefig(['Overlay_extractedParticles_' savename '.fig']);

% figure('Position',[1000 400 600 600])
% imshow(I,[minWF maxWF]); hold on;
% CentO = [];
% CentO(:,1) = locs(:,xCol)/CFX(:,1);
% CentO(:,2) = locs(:,yCol)/CFX(:,1);;
% scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');

% 
% cd(savepath);
% savefig(['Overlay_allLocs' savename '.fig']);

% Plot for each Particle (1) the integrate intensitz vs (2) the nbr of locs 

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
%% 
close all
cd(H5Folder);
