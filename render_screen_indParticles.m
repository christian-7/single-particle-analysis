%% Use this script to choose the best filtering parameters for a dataset

% Plot histogram of all relevant parameters
% Choose filtering range
% Render result before/after
% option to save the image as 32-bit tiff

% INPUT: Cent

%% Select a particle and define the columns

i = 200;

% all_Cent = Cent_filt;

locs = [];
locs = all_Particles{i,1};

xCol               = 1;
yCol               = 2;
framesCol          = 3;
uncertaintyCol     = 4;
photonsCol         = 5;
LLCol              = 6;
sigmaCol           = 7;

%% PLot filter parameters to set the filter 

close all;

figure('Position',[1000 300 700 700]); 
subplot(2,2,1);
hist(locs(:,photonsCol),50);
title('Photons');
subplot(2,2,2);
hist(locs(:,sigmaCol),50);
title('Sigma');
subplot(2,2,3);
hist(locs(:,LLCol),50);
title('LL');
subplot(2,2,4);
hist(locs(:,uncertaintyCol),50);
title('Uncertainty');


%% Plot locs per frame

locs_p_frame =[];
count = 1;

for i=min(locs(:,framesCol)):max(locs(:,framesCol));
    
    vx = find(locs(:,framesCol)==i);
    
    locs_p_frame(count,1) = length(vx);
    locs_p_frame(count,2) = i;
    count=count+1;
    
end

figure('Position',[100 500 1200 300])
subplot(1,3,1)
scatter(locs_p_frame(:,2),locs_p_frame(:,1),1,'.k')
xlabel('frames');
ylabel('locs per frame');
box on;
axis([0 3e4 0 5])

subplot(1,3,2)
bins = 0 : 1 : 3;
h = hist(locs_p_frame(:,1),bins)
bar(bins, h/sum(h),'k')
xlabel('locs per frame');
ylabel('frequency');
box on;

subplot(1,3,3)
scatter(locs(:,framesCol),locs(:,uncertaintyCol),1,'.k')
xlabel('frames');
ylabel('uncertainty');
box on;

%% Apply the filter

close all

% Set Filter parameters 

minFrame            = 2000;
maxLL               = 150;
MinPhotons          = 1500;
Minsigma            = 120; 
Maxsigma            = 180;
Maxuncertainty      = 15;

filter   = find(locs(:,sigmaCol) < Maxsigma & locs(:,sigmaCol) > Minsigma & locs(:,photonsCol) > MinPhotons & locs(:,uncertaintyCol) < Maxuncertainty & locs(:,framesCol) > minFrame );

subsetLL = locs(filter,1:end);

fprintf('\n -- Data Filtered %f of the locs are left --\n', ((length(subsetLL)/length(locs))));

%% Calculate 2D histogram 

close all

pxlsize = 10;

heigth  = round((max(subsetLL(:,yCol))-min(subsetLL(:,yCol)))/pxlsize);
width   = round((max(subsetLL(:,xCol))-min(subsetLL(:,xCol)))/pxlsize);

im = hist3([locs(:,yCol),locs(:,xCol)],[heigth width]);
imfiltered = hist3([subsetLL(:,yCol),subsetLL(:,xCol)],[heigth width]);


figure('Position',[500 300 1200 600])
subplot(1,2,1)
imshow(imgaussfilt(im, 1),[0.3 50]);
title('Not Filtered')
colormap('hot');

subplot(1,2,2)
imshow(imgaussfilt(imfiltered, 1),[0.3 20]);
title('Filtered')
colormap('hot');


%% Save 32-bit image

cd(savepath)

name = [savename '_rendered_' num2str(pxlsize) 'nm_per_pxl.tiff'];  

I32=[];
I32=uint32(im);

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


cd(locpath);

fprintf('\n -- Image saved --\n');
