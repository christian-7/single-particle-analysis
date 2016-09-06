function c_render_image(Locpath,locName,Segpath,segName,savename,impath,outputFileName1);

cd(Locpath);

file = fopen(locName);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol = strmatch('x [nm]',h);
yCol = strmatch('y [nm]',h);
frameCol = strmatch('frame',h);
photonsCol = strmatch('intensity [photon]',h);
sigmaCol = strmatch('sigma [nm]',h);
uncertaintyCol = strmatch('uncertainty [nm]',h);

% Load the segemented centrioles

cd(Segpath);
load(segName);

fprintf('-- 1/4 Header and Data Loaded --\n')



%% Concatenate all Centrioles in one variable

c1 = [];
c2 = [];
c3 = [];
c4 = [];
c5 = [];
c6 = [];
all_locs = [];

for i=1:length(Cent);
    
   c1 = vertcat(c1, Cent{i,2}(:,1));
   c2 = vertcat(c2, Cent{i,2}(:,2));
   c3 = vertcat(c3, Cent{i,2}(:,3));
   c4 = vertcat(c4, Cent{i,2}(:,4));
   c5 = vertcat(c5, Cent{i,2}(:,5));
   c6 = vertcat(c6, Cent{i,2}(:,6));
    
end

all_locs(:,1) = c1;
all_locs(:,2) = c2;
all_locs(:,3) = c3;
all_locs(:,4) = c4;
all_locs(:,5) = c5;
all_locs(:,6) = c6;

clear c1 c2 c3 c4 c5 c6 

save(savename,'all_locs');

fprintf('-- 2/4 Saved smart merged data --\n')

%% Determine the bin for each localization and wheight the pixel acc to the photon count

pxlsize=10;

heigth=round((max(all_locs(:,2))-min(all_locs(:,2)))/pxlsize);
width=round((max(all_locs(:,1))-min(all_locs(:,1)))/pxlsize);

imWMerging = hist3([all_locs(:,1),all_locs(:,2)],[width heigth]); 

% Generate a matrix where each pxl contains the sum of the photons from all



% locs that fall into that pxl

tic
binMerged2 = zeros(width,heigth);       % matrix with the size of the 2D histogram image

for i=1:length(all_locs);               % for all molecules, find x and y pixel
    
    x = round((all_locs(i,1)/pxlsize)-(min(all_locs(:,1))/pxlsize));
    y = round((all_locs(i,2)/pxlsize)-(min(all_locs(:,2))/pxlsize));
    
    if      x == 0;
            x = 1;
    elseif  y == 0;      
            y = 1;
    end

    binMerged2(x,y) = binMerged2(x,y)+all_locs(i,3);
    
    clear x y
    
end
toc

imWMerging2 = times(imWMerging, binMerged2);

fprintf('-- 3/4 Generated photon wheighted 2D Histogram --\n')

%% Write 32-bit file

cd(impath);

I32=[];
I32=uint32(imWMerging2);

t = Tiff(outputFileName1,'w');

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

fprintf('-- 4/4 Saved 32-bit image --\n')

end
