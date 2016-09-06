clear, clc, close all

%% Read Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pxl=107.99;                                                                         % Pixel size in nm
WF_name  = 'Kog1_GFP_30C_NB_A647_WF11_MMStack_Pos0.tiff';                              % name of the WF image
locName  = 'Yeast_Kog1_GFP_30C_NB_A647_10_MMStack_Pos0_locResults_processed_DC_TS.csv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load image

% cd(WFpath);
I=imread(WF_name);

% Load data set
% cd(Locpath);
% name2=[base,'_MMStack_Pos0_locResults_DC.dat'];
locs=dlmread(locName,',',1,0);

% Load header

file = fopen(locName);
line = fgetl(file);
h = regexp( line, ',', 'split' );

x = strmatch('"x [nm]"',h);
y = strmatch('"y [nm]"',h);
LL = strmatch('"loglikelihood"',h);

fprintf('\n -- Data loaded --\n')


%% Adjust segmentation parameters

close all

% blur the image
figure('name','Raw Image'), imshow(I);
% adjust the contrast of the raw image
I2 = imadjust(I,[0.3 0.5],[]);
figure('name','Image after adjusted contrast'), imshow(I2);

G = fspecial('gaussian',[7 7],50); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imG = imfilter(I2,G,'same');
figure('name','Image after adjusted contrast, blurring'), imshow(imG);

% adjust the background
I3 = imadjust(imG,[0.01 0.3],[]);
figure('name','Image after adjusted contrast, blurring, adjusted contrast'), imshow(I3);

% Make binary image
bin = im2bw(I3,0.3);
figure('name','Binary image result'),imshow(bin)
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
box_size = 20;

for k=1:length(B)
    
    if      min(B{k,1}(:,1))<box_size+1 
        
            Particles_WF{k,1} = I((1):(max(B{k,1}(:,1))+box_size),(min(B{k,1}(:,2))-box_size):(max(B{k,1}(:,2))+box_size));    
            
    elseif  min(B{k,1}(:,2))<box_size+1 
        
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
            Center(k,1) = (((max(B{k,1}(:,1))-min(B{k,1}(:,1)))/2)+min(B{k,1}(:,1)))*(pxl);%*0.9723); Center of the segemneted spot in nm
            Center(k,2) = (((max(B{k,1}(:,2))-min(B{k,1}(:,2)))/2)+min(B{k,1}(:,2)))*(pxl);%*0.9751); Center of the segemneted spot in nm
            Center(k,3) = max(pdist(B{k,1}))/2;                                                         %  size of the box, measure the max distance as input for the box size
            
end

CFX = (max(locs(:,x)/pxl))./size(I);
CFY = (max(locs(:,y)/pxl))./size(I);

center2(:,1)    =   Center(:,2)*CFX(:,1);
center2(:,2)    =   Center(:,1)*CFY(:,1);
center2(:,3)    =   Center(:,3);
Center=center2;

fprintf('\n -- Particles extracted --\n')

%% Build box around each Center and copy locs into separate variable -> structure Cent

Cent={}; 
box_size = 200; 
    
for i=1:length(Center);
    
    vx1=Center(i,1)+Center(i,3)*pxl;
    vx2=Center(i,1)-Center(i,3)*pxl;
    vy1=Center(i,2)+Center(i,3)*pxl;
    vy2=Center(i,2)-Center(i,3)*pxl;
    
%     vx1=Center(i,1) + box_size;
%     vx2=Center(i,1) - box_size;
%     vy1=Center(i,2) + box_size;
%     vy2=Center(i,2) - box_size;
    
    vx=find(locs(:,x) > vx2 & locs(:,x)<vx1);
    subset=locs(vx,1:9);

    vy=find(subset(:,y)>vy2 & subset(:,y)<vy1);

    Cent{i,1} = subset(vy,1:end);  
    Cent{i,2} = intI(i);
    Cent{i,3} = Particles_WF{i,1};
    
end 
    
fprintf('\n -- %f Centrioles selected from localitaion dataset --\n', length(Cent))

%% Overlay with Widefield image

figure
imshow(I); hold on;
for i=1:length(Cent);
scatter(Cent{i}(:,1)/pxl,Cent{i}(:,2)/pxl,1,'r');
hold on;
end

figure
for i=1:length(Cent);
scatter(length(Cent{i}),Cent{i,2},5,'filled','r');hold on;
xlabel('Nbr of localizations');
ylabel('GFP integrated intensity');
box on;
axis square;
end

% Use this to plot all localizations

% figure
% imshow(I); hold on;
% scatter(locs(:,x)/pxl,locs(:,y)/pxl,1,'r');

% savefig('Overlay_NB_GFP_12.fig');

save('Extracted_Particles_4.mat','Cent');

%% Render extracted Particles, Save as 32-bit tiff

pxlsize=10;
cd('Z:\Christian-Sieben\data_HTP\2016-07-28_Yeast\locResults\Yeast_Kog1_GFP_30C_NB_A647_12\segmented');

for i=1:length(Cent);

    if length(Cent{i,1})>200;
        
        
        heigth = round((max(Cent{i,1}(:,2))-min(Cent{i,1}(:,2)))/pxlsize);
        width  = round((max(Cent{i,1}(:,1))-min(Cent{i,1}(:,1)))/pxlsize);
        
%         width =  400/pxlsize;

rendered = hist3([Cent{i,1}(:,2),Cent{i,1}(:,1)],[heigth width]);

name = ['image_10nm' num2str(i) '.tiff'];

I32=[];
I32=uint32(rendered);

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

    else end

end

%% Render extracted Particles, Save as 32-bit tiff

pxlsize=10;
cd('Z:\Christian-Sieben\data_HTP\2016-07-28_Yeast\locResults\Yeast_Kog1_GFP_30C_NB_A647_7\segmented');

% Determine the box size form the largest particle

im_size = [];

for i=1:length(Cent);

im_size(i,1) = round((max(Cent{i,1}(:,2))-min(Cent{i,1}(:,2)))/pxlsize);
im_size(i,2) = round((max(Cent{i,1}(:,1))-min(Cent{i,1}(:,1)))/pxlsize);

end

% Build 2Dhistogram and embed it into a larger black matrix

for i=1:length(Cent);

    if length(Cent{i,1})>200;
        
        empty  = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));
        center = round(length(empty)/2);

heigth = round((max(Cent{i,1}(:,2))-min(Cent{i,1}(:,2)))/pxlsize);
width  = round((max(Cent{i,1}(:,1))-min(Cent{i,1}(:,1)))/pxlsize);

rendered = hist3([Cent{i,1}(:,2),Cent{i,1}(:,1)],[heigth width]);

empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = rendered;

name = ['image_10nm_' num2str(i) '.tiff'];

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

    else end

end



