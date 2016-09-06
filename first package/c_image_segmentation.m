function c_image_segmentation(WFpath,WF_name,Locpath,locName,savepath,savename) 

%% Read Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pxl=107.99;                                                                 % Pixel size in nm

%%%%%%%%%%%%%%%%% Manual Input %%%%%%%%%%%%%%%%%%%%%%%%%

WFpath      = 'Y:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647';
WF_name     = 'SAS6_FOV1.tif';                    
base        = 'humanCent_aTub_NB_A647_3';   
Locpath     = 'Y:\Christian-Sieben\data_HTP\2016-06-03_humanCent_Sas6_A647\locsResults_Sas6\FOV1_Sas6_1000mW_10ms_A647_1';
locName     = 'FOV1_Sas6_1000mW_10ms_A647_1_MMStack_Pos0_locResults_DC.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load image

cd(WFpath);
I=imread(WF_name);

% Load data set
cd(Locpath);
% name2=[base,'_MMStack_Pos0_locResults_DC.dat'];
locs=dlmread(locName,',',1,0);

% Load header

file = fopen(locName);
line = fgetl(file);
h = regexp( line, ',', 'split' );

x = strmatch('x [nm]',h);
y = strmatch('y [nm]',h);
LL = strmatch('loglikelihood',h);

fprintf('\n -- Data loaded --\n')


%% Find and Plot Center of mass for each object

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

% blur the image

figure;imshow(I);

G = fspecial('gaussian',[7 7],50); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imG = imfilter(I,G,'same');

figure;imshow(imG);

% adjust the background

% I2 = imadjust(imG,[0.1 0.15],[]); % (I,[low_in; high_in],[low_out; high_out])
I2 = imadjust(I,[0.4 0.5],[]);
figure;imshow(I2);


% Make binary image

bin = im2bw(I2,0.3);
[B,L,N,A] = bwboundaries(bin);

% Plot result

%Find the center of each particle an dtransform into a x,y coordinate

Center=[];

for k=1:length(B)
        
            boundary    = B{k};
            Center(k,1) = (((max(B{k,1}(:,1))-min(B{k,1}(:,1)))/2)+min(B{k,1}(:,1)))*(pxl);%*0.9723);
            Center(k,2) = (((max(B{k,1}(:,2))-min(B{k,1}(:,2)))/2)+min(B{k,1}(:,2)))*(pxl);%*0.9751);
            Center(k,3) = max(pdist(B{k,1}))/2;                                                         %  measure the max distance as input for the box size
            
end

CFX = (max(locs(:,x)/pxl))./size(I);
CFY = (max(locs(:,y)/pxl))./size(I);

center2(:,1)    =   Center(:,2)*CFX(:,1);
center2(:,2)    =   Center(:,1)*CFY(:,1);
center2(:,3)    =   Center(:,3);
Center=center2;

fprintf('\n -- Centrioles identified --\n')

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

    Cent{i,1}=subset(vy,1:end);
       
end 
    
fprintf('\n -- %f Centrioles selected from localitaion dataset --\n',length(Cent))

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

%% Save variable Cent
cd(savepath);

save(savename,'Cent');

fprintf('\n -- File Saved --\n')

end
