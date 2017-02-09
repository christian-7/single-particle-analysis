%% 

Cent = all_Cent;

pxlsize = 10;

cd('Z:\Christian-Sieben\data_HTP\2016-10-13_humanCentriole_Cep63\analysis\rendered');

% Determine the box size form the largest particle

im_size = [];

for i=1:length(Cent);
    
im_size(i,1) = round((max(Cent{i,1}(:,2))-min(Cent{i,1}(:,2)))/pxlsize);
im_size(i,2) = round((max(Cent{i,1}(:,1))-min(Cent{i,1}(:,1)))/pxlsize);
        
end

count=1;

for i = 1:length(Cent);
    
            
        heigth =round((max(Cent{i,1}(:,2)) - min(Cent{i,1}(:,2)))/pxlsize);
        width = round((max(Cent{i,1}(:,1)) - min(Cent{i,1}(:,1)))/pxlsize);
        
        rendered = hist3([Cent{i,1}(:,2),Cent{i,1}(:,1)],[heigth width]);
        
        empty  = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));
        center = round(length(empty)/2);

        empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = rendered;        
 
name32rendered  = ['image_10nm_32bit_rendered_FOV_' num2str(i) '.tiff'];
% name16          = ['image_10nm_16bit_' num2str(i),'_',num2str(j) '.tiff'];    
% name32          = ['image_10nm_32bit_' num2str(i),'_',num2str(j) '.tiff'];

I32=[];
I32=uint32(imgaussfilt(empty,1));

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
            
end
   