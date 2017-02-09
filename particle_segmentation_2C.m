%% Particle segmentation from 2C data

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

pxl       = 107.99;                                                                 % Pixel size in nm
number    = 12;                                                                      % image number
filetype  = 1;                                                                      % 1 for TS, 2 for bstore
wfnumber  = 13;

%%%%%%%%%%%%%%%%% Manual Input %%%%%%%%%%%%%%%%%%%%%%%

WFpath     = ['Z:\Christian-Sieben\data_HTP\2016-07-28_Yeast\Kog1_GFP_30C_WF' num2str(wfnumber)];
WF_name    = ['Kog1_GFP_30C_WF' num2str(wfnumber) '_MMStack_Pos0.ome.tif'];         

WFpath2      = ['Z:\Christian-Sieben\data_HTP\2016-07-28_Yeast\Kog1_GFP_30C_NB_A647_WF' num2str(wfnumber)];
WF_name2     = ['Kog1_GFP_30C_NB_A647_WF' num2str(wfnumber) '_MMStack_Pos0.ome.tif'];      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Locpath     = ['Z:\Christian-Sieben\data_HTP\2016-07-28_Yeast\locResults\Yeast_Kog1_GFP_30C_NB_A647_' num2str(number)];
locName     = ['Yeast_Kog1_GFP_30C_NB_A647_' num2str(number) '_MMStack_Pos0_locResults_DC_TS.csv'];

savename    = ['Yeast_Kog1_GFP_30C_NB_A647' num2str(number) '_extractedParticles_2nd'];
savepath    =  'Z:\Christian-Sieben\data_HTP\2016-07-28_Yeast\locResults\newAnalysis';

savepath_Imgages = 'Z:\Christian-Sieben\data_HTP\2016-07-28_Yeast\locResults\newAnalysis';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(Locpath);

%% Load image

% GFP Channel

cd(WFpath);
I   =   imread(WF_name);

% Storm Channel

cd(WFpath2);
ICh2  =  imread(WF_name2);

% Load localization data set

cd(Locpath);
locs=dlmread(locName,',',1,0);

% Load header

file    = fopen(locName);
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

close all

minWF_GFP = 389;
maxWF_GFP = 8000;

% blur the image
figure('Position',[10 600 500 500],'name','Raw GFP Image'), imshow(I,[minWF_GFP maxWF_GFP],'InitialMagnification','fit');

% adjust the contrast of the raw image
I2 = imadjust(I,[0.06 0.07],[]);
figure('Position',[600 600 500 500],'name','Image after adjusted contrast'), imshow(I2,'InitialMagnification','fit');

% lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
G = fspecial('gaussian',[3 3],60); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imG = imfilter(I2,G,'same');
figure('Position',[1200 600 500 500],'name','Image after adjusted contrast, blurring'), imshow(imG,'InitialMagnification','fit');

% adjust the background
I3 = imadjust(imG,[0.08 0.1],[]);
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

% Extract the integrated intensity of the A647 WF image for each ROI

intI2 = [];

for i = 1:length(B);
    
    intI2(i,1) = sum(sum(ICh2(min(B{i,1}(:,1)):max(B{i,1}(:,1)),min(B{i,1}(:,2)):max(B{i,1}(:,2)))));
    
end

% Extract subimages from GFP channel

Particles_WF  ={};
Particles_WF2 ={};

box_size = 20;

for k=1:length(B)
    
%   the border particles WF ROIs are saved without the box
    
    
       if      min(B{k,1}(:,1))<box_size+1 | min(B{k,1}(:,2)) < box_size+1 | max(B{k,1}(:,1))+box_size>length(I) | max(B{k,1}(:,2))+box_size>length(I)
           
               Particles_WF{k,1} = I((min(B{k,1}(:,1))):(max(B{k,1}(:,1))),(min(B{k,1}(:,2))):(max(B{k,1}(:,2))));
               Particles_WF2{k,1} = ICh2((min(B{k,1}(:,1))):(max(B{k,1}(:,1))),(min(B{k,1}(:,2))):(max(B{k,1}(:,2))));
           
       else
               Particles_WF{k,1}  = I((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));  
               Particles_WF2{k,1} = ICh2((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));   
            
       end
    
    
%     if      min(B{k,1}(:,1))<box_size+1 % lower side
%         
%             Particles_WF{k,1} = I((1):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));    
%             
%     elseif  min(B{k,1}(:,2))<box_size+1 % left side
%         
%             Particles_WF{k,1} = I((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(1):(max(B{k,1}(:,2))+box_size));    
%     
%     elseif  max(B{k,1}(:,1))+box_size>length(I) | max(B{k,1}(:,2))+box_size>length(I) % upper or right side
%         
%             Particles_WF{k,1} = I((min(B{k,1}(:,1))-box_size):(length(I)),(min(B{k,1}(:,2))-box_size):(length(I)));
%         
%     else
%             Particles_WF{k,1} = I((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));   
%     end
           
end

% Extract subimages from A647 channel

% Particles_WF2={};
% box_size = 20;
% 
% for k=1:length(B)
%     
%     if      min(B{k,1}(:,1))<box_size+1 
%         
%             Particles_WF2{k,1} = ICh2((1):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));    
%             
%     elseif  min(B{k,1}(:,2))<box_size+1 
%         
%             Particles_WF2{k,1} = ICh2((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(1):(max(B{k,1}(:,2))+box_size));    
%     
%     elseif  max(B{k,1}(:,1))+box_size>length(ICh2) | max(B{k,1}(:,2))+box_size>length(ICh2)
%         
%             Particles_WF2{k,1} = ICh2((min(B{k,1}(:,1))-box_size):(length(ICh2)),(min(B{k,1}(:,2))-box_size):(length(ICh2)));
%         
%     else
%         
%             Particles_WF2{k,1} = ICh2((min(B{k,1}(:,1))-box_size):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));   
%     end
%            
% end

%Find the center of each particle and transform into an X,Y coordinate

Center=[];

for k=1:length(B)
        
            boundary    = B{k};
            Center(k,1) = (((max(B{k,1}(:,1))-min(B{k,1}(:,1)))/2)+min(B{k,1}(:,1)))*(pxl);             % Center of the segemneted spot in nm
            Center(k,2) = (((max(B{k,1}(:,2))-min(B{k,1}(:,2)))/2)+min(B{k,1}(:,2)))*(pxl);             % Center of the segemneted spot in nm
            Center(k,3) = max(pdist(B{k,1}))/1;                                                         % Ssize of the box, measure the max distance as input for the box size
            
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
    Cent{i,4} = intI2(i);
    Cent{i,5} = Particles_WF2{i,1};
       
end 
    
fprintf('\n -- %f Particles selected from localitaion dataset --\n',length(Cent))

%% Overlay with Widefield image and plot correlations

% the coordinates need be backtransformed into pxl coordinates
% 1. divide by correction factor
% 2. divide trough pxl size


figure('Position',[10 600 500 500],'name','Extracted Particles');
imshow(I,[minWF_GFP maxWF_GFP]); hold on;

for i=1:length(Cent);
    
CentO=[];
CentO(:,1) = Cent{i,1}(:,x)/CFX(:,1);
CentO(:,2) = Cent{i,1}(:,y)/CFY(:,1);
    
scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');

hold on;
end

savefig('Overlay_extracted_particles_onGFP.fig');


figure('Position',[1000 600 500 500],'name','All localizations');
imshow(I,[minWF_GFP maxWF_GFP]); hold on;
CentO = [];
CentO(:,1) = locs(:,x)/CFX(:,1);
CentO(:,2) = locs(:,y)/CFX(:,1);;
scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');


savefig('Overlay_all_Locs_onGFP.fig');

% Plot for each Particle (1) the integrate intensity vs (2) the nbr of locs 

figure('Position',[400 100 900 300],'name','# of Locs vs. GFP Intensity');

subplot(1,3,1);
for i=1:length(Cent);
scatter(length(Cent{i}),Cent{i,2},5,'filled','r');hold on;
xlabel('Nbr of localizations');
ylabel('GFP WF integrated intensity');
box on;
axis square;
end

subplot(1,3,2);
for i=1:length(Cent);
scatter((Cent{i,2}),Cent{i,4},5,'filled','r');hold on;
xlabel('GFP intensity');
ylabel('A647 intensity');
box on;
axis square;
end

subplot(1,3,3);
for i=1:length(Cent);
scatter((Cent{i,2}),Cent{i,4}/8300,5,'filled','r');hold on;
xlabel('GFP intensity');
ylabel('A647 molecules');
box on;
axis square;
end


savefig('WF_correlation.fig');
save(savename,'Cent');

%% DBSCAN particle size Filter

close all;

tic
fprintf('\n -- DBSCAN started --\n')

Cent{length(Cent),6} = [];

for m = 1:length(Cent);

dataDBS      = [];
dataDBS(:,1) = Cent{m,1}(:,x); % x in mum
dataDBS(:,2) = Cent{m,1}(:,y); % y in mum

if isempty(dataDBS)
   Cent{m,6} = [];
else

% Run DBSCAN on each particle 

k   = 10;                                                 % minimum number of neighbors within Eps
Eps = 20;                                                 % minimum distance between points, nm

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


% subsetP = [];
% 
% subsetP(:,1)    = dataDBS(coreBorder,1);
% subsetP(:,2)    = dataDBS(coreBorder,2);
% subsetP(:,3)    = class2(coreBorder);
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
    
    if length(vx)>50;

    Cent{m,6}{i,1} = subset(vx,1:end);
         
    else    end
    
end

end

end

% figure('Position',[10 600 400 400],'name','# of Locs vs. GFP Intensity after DBSCAN');
% 
% for i=1:length(Cent);
%     
%     for j = 1:length(Cent{i,6});
%     
%     if isempty(Cent{i,6}{j,1})
%     else
%     
%         scatter(length(Cent{i,6}{j,1}),Cent{i,2},5,'filled','r');hold on;
%         xlabel('Nbr of localizations');
%         ylabel('WF GFP integrated intensity');
%         box on;
%         axis square;
%     end
%     end
%     
% end

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

%% Plot all particles in subplots

NbrSubplots = [];

for i = 1:length(Cent);
    
    NbrSubplots(i) = length(Cent{i,6});
end

NbrSubplots = round(sqrt(sum(NbrSubplots)))+1;

count = 1;
figure('Position',[100 100 1500 1000])

for i = 1:length(Cent);
    
    for j = 1:length(Cent{i,6});
    
            if isempty(Cent{i,6}{j,1})
            else
            
        subplot(NbrSubplots,NbrSubplots,count);
        scatter(Cent{i,6}{j,1}(:,x),Cent{i,6}{j,1}(:,y),10,'filled','black');
        title(['ROI :', num2str(i) ' Particle :', num2str(j)],'FontSize',9);
        count = count+1;
        box on
        axis off
        
            end             
    end
   
end

%% Plot each Isolated particles with the WF ROI

NbrSubplots = [];

for i = 1:length(Cent);
    
    NbrSubplots(i) = length(Cent{i,6});
end

NbrSubplots = round(sqrt(sum(NbrSubplots)))+1;

figure('Position',[100 100 1000 1000])
count = 1;

for i=1:length(Cent);
    
       
for j = 1:length(Cent{i,6});
    
    subplot(NbrSubplots+1,NbrSubplots,count);
    imshow(Cent{i,3},[minWF_GFP maxWF_GFP]);hold on;
    
            if isempty(Cent{i,6}{j,1})
            else
                
CentO=[];
CentO(:,1) = Cent{i,1}(:,x)/CFX(:,1);
CentO(:,2) = Cent{i,1}(:,y)/CFY(:,1);

CentO(:,1) = CentO(:,1)-min(CentO(:,1)) + box_size*pxl; % this normalizes to 0
CentO(:,2) = CentO(:,2)-min(CentO(:,2)) + box_size*pxl;

scatter(CentO(:,1)/pxl,CentO(:,2)/pxl,1,'r');hold on; % Position over full WF
title(['ROI :', num2str(i) ' Particle :', num2str(j)],'FontSize',9);

            end
            
count = count+1;            
end

end




%% Save variable Particles

% 1. localization data
% 2. WF integrated intensity GFP
% 3. WF Image subset GFP
% 4. WF integrated intensity A657
% 5. WF Image subset A647
% 6. DBSCAN filtered particles

close all

cd(savepath);
% assignin('base',savename,Cent);
% save(savename,savename);
save(savename,'Cent');

close all
fprintf('\n -- File Saved --\n')
cd(Locpath);
%% Save all Particles as Rendered images

pxlsize = 10;
cd(savepath_Imgages);

% Determine the box size form the largest particle

im_size = [];
count = 1;

for i=1:length(Cent);
    
    for j = 1:length(Cent{i,6});
        
        if isempty(Cent{i,6}{j,1})
        else

im_size(count,1) = round((max(Cent{i,6}{j,1}(:,y))-min(Cent{i,6}{j,1}(:,y)))/pxlsize);
im_size(count,2) = round((max(Cent{i,6}{j,1}(:,x))-min(Cent{i,6}{j,1}(:,x)))/pxlsize);
        
count = count + 1;

        end
    end
end

count=1;

for i = 1:length(Cent);
    
    for j = 1:length(Cent{i,6});
    
            if isempty(Cent{i,6}{j,1})
                
            else
            
        heigth =round((max(Cent{i,6}{j,1}(:,x)) - min(Cent{i,6}{j,1}(:,x)))/pxlsize);
        width = round((max(Cent{i,6}{j,1}(:,y)) - min(Cent{i,6}{j,1}(:,y)))/pxlsize);
        
        rendered = hist3([Cent{i,6}{j,1}(:,y),Cent{i,6}{j,1}(:,x)],[heigth width]);
        
        empty  = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));
        center = round(length(empty)/2);

        empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = rendered;
 

name = [savename 'image_10nm_' num2str(i),'_',num2str(j) '.tiff'];  

I32=[];
I32=uint32(empty);

t = Tiff(name,'w');

tagstruct.ImageLength     = size(I32,1);
tagstruct.ImageWidth      = size(I32,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)

t.write(I32);
t.close()

count=count+1;
            end
    end
end

fprintf('\n -- Images Saved --\n')
cd(Locpath);
