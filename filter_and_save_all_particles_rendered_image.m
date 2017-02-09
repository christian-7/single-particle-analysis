%% Use this script to filter, render and save all particles in a dataset

% INPUT: Cent
% Choose filter parameters with render_screem_indParticles.m
clear, clc
%% 
% Cent = all_Particles;

pxlsize = 10;

savepath = 'Z:\Christian-Sieben\data_HTP\2016-04-01_humanCentriole_aTubNB_Sas6\2C STORM analysis\rendered_images';

%% Set the filter  
 
minFrame            = 2000;
maxLL               = 150;
MinPhotons          = 1500;
Minsigma            = 130; 
Maxsigma            = 150;
Maxuncertainty      = 15;


xCol               = 1;
yCol               = 2;
framesCol          = 4;
uncertaintyCol     = 5;
photonsCol         = 6;
LLCol              = 8;
sigmaCol           = 9;

%% Determine the box size form the largest particle

im_size = [];

for i=1:length(Cent);
    
if isempty(Cent{i,1})==1;
    
else
    
im_size(i,1) = round((max(Cent{i,1}(:,2))-min(Cent{i,1}(:,2)))/pxlsize);
im_size(i,2) = round((max(Cent{i,1}(:,1))-min(Cent{i,1}(:,1)))/pxlsize);

end
        
end

fprintf(['Image Size ' num2str(max(im_size(:,1))) ' by ' num2str(max(im_size(:,2)))])

%% Render each particle into 32-bit and save in savepath
tic
cd(savepath);
count=1;
Channel = 5;

for i = 1:length(Cent);
    
    if isempty(Cent{i,1})==1;
    
    else
    
        filter   = find(Cent{i,Channel}(:,sigmaCol) < Maxsigma & Cent{i,Channel}(:,sigmaCol) > Minsigma & Cent{i,Channel}(:,photonsCol) > MinPhotons & Cent{i,Channel}(:,uncertaintyCol) < Maxuncertainty & Cent{i,Channel}(:,framesCol) > minFrame );
        subsetLL = Cent{i,Channel}(filter,1:end);
        
        if isempty(subsetLL)==1;
            
        else
      
        heigth =round((max(subsetLL(:,2)) - min(subsetLL(:,2)))/pxlsize);
        width = round((max(subsetLL(:,1)) - min(subsetLL(:,1)))/pxlsize);
        
        rendered = hist3([subsetLL(:,2),subsetLL(:,1)],[heigth width]);
        
        empty  = zeros(round(max(max(im_size))*1.1),round(max(max(im_size))*1.1));
        
        size_DC    = size(empty);
        center_X   = round(size_DC(2)/2);
        center_Y   = round(size_DC(1)/2);     

        empty(round(center_Y-heigth/2):(round(center_Y-heigth/2))+heigth-1,round(center_X-width/2):(round(center_X-width/2))+width-1) = rendered;        
 
        name32rendered  = ['image_10nm_32bit_rendered_Ch2_Cent_' num2str(i) '.tiff'];

I32=[];
I32=uint32(empty);

t = Tiff(name32rendered,'w');
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

X = [' Finished ',num2str(i),' of ',num2str(length(Cent)),];
clc;disp(X); 

        end
    end       
end
 toc  