function c_render_image_wo_merging(Locpath,locName,pxlsize,impath,outputFileName1) 
% Load file

cd(Locpath);
locs=dlmread(locName,',',1,0);

file = fopen(locName);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol = strmatch('x [nm]',h);
yCol = strmatch('y [nm]',h);
frameCol = strmatch('frame',h);
photonsCol = strmatch('intensity [photon]',h);
sigmaCol = strmatch('sigma [nm]',h);
uncertaintyCol = strmatch('uncertainty [nm]',h);


fprintf('-- 1/2 Data Loaded --\n')

%% Filter, build historgam, blur -> show

% Find width and heigth
heigth=round((max(locs(:,yCol))-min(locs(:,yCol)))/pxlsize);
width=round((max(locs(:,xCol))-min(locs(:,xCol)))/pxlsize);

% Calculate 2D histogram --> 10 nm/pxl

im=hist3([locs(:,xCol),locs(:,yCol)],[width heigth]); % heigth x width
I16 = uint16(round(im*65535));
I32 = uint32(im);

% Apply Gaussian Filter

G = fspecial('gaussian',[4 4],1.5); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imG = imfilter(I32,G,'same');

%% Write 32-bit Files

cd(impath);

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
t.close();

fprintf('-- 2/2 Image saved -- \n')

end

