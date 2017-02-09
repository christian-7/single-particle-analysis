%% First Classification - TopOrSide / Unclassified

% Start this after the density filtering

Cent = Cent_filt;

TopOrSide       = [];
unclassified    = [];
          
% Top = {};
% Side= {};

clc, close all

for i = 1:length(Cent);

for j = 1:length(Cent{i,6});
    
            if isempty(Cent{i,6}{j,1})
            else
            
        figure('Position',[500 300 300 300]);
        
        scatter(Cent{i,6}{j,1}(:,1)-min(Cent{i,6}{j,1}(:,1)),Cent{i,6}{j,1}(:,2)-min(Cent{i,6}{j,1}(:,2)),10,'filled','black');
        title(['ROI :', num2str(i) ' Particle :', num2str(j)],'FontSize',9);
        axis([0 max(Cent{i,6}{j,1}(:,1))-min(Cent{i,6}{j,1}(:,1)) 0 max(Cent{i,6}{j,1}(:,2))-min(Cent{i,6}{j,1}(:,2))])
        box on
        axis on  
        
        fprintf('Mouse --> Top or Side');
        fprintf('\n');
        fprintf('\n');
        fprintf('Key --> unclassified');
                
        w = waitforbuttonpress;                                     % 0 if it detects a mouse button click / 1 if it detects a key press
        
        if w == 0;  
            
        TopOrSide = vertcat(TopOrSide,Cent{i,6}(j,1));              % Put on the shown cluster into the new variable  
        disp('Top')
        
        else
            
        unclassified = vertcat(unclassified,Cent{i,6}(j,1));        % Put on the shown cluster into the new variable 
        disp('Side')
        end

            end  
    close all  , clc     
    end
   
    
end

%% Second Classification - Top / Side
       
Top = {};
Side= {};

for j = 1:length(TopOrSide);

    
            if isempty(TopOrSide{j,1})
            else
            
        figure('Position',[500 300 300 300]);
        scatter(TopOrSide{j,1}(:,1)-min(TopOrSide{j,1}(:,1)),TopOrSide{j,1}(:,2)-min(TopOrSide{j,1}(:,2)),10,'filled','black');
        title([' Particle :', num2str(j)],'FontSize',9);
        axis([0 max(TopOrSide{j,1}(:,1))-min(TopOrSide{j,1}(:,1)) 0 max(TopOrSide{j,1}(:,2))-min(TopOrSide{j,1}(:,2))])
        box on
        axis on  
        
        fprintf('Mouse --> Top');
        fprintf('\n');
        fprintf('\n');
        fprintf('Key --> Side');
             
        w = waitforbuttonpress;                     % 0 if it detects a mouse button click / 1 if it detects a key press
        
        if w == 0;  
            
        Side = vertcat(Side,TopOrSide{j,1});        % Put on the shown cluster into the new variable  
        disp('Top')
        
        else
            
        Top = vertcat(Top,TopOrSide{j,1});          % Put on the shown cluster into the new variable 
        disp('Side')
        end

            end  
    close all  , clc     
   
   
    
end


%% Render and save all particles

pxlsize = 10;

cd('Z:\Christian-Sieben\data_HTP\2016-09-29_Yeast_Kog1_GFP\analysis\rendered');

% Top or Side to render

Cent = Side;

% Determine the box size form the largest particle

im_size = [];

for j=1:length(Cent);
       
im_size(j,1) = round((max(Cent{j,1}(:,2))-min(Cent{j,1}(:,2)))/pxlsize);
im_size(j,2) = round((max(Cent{j,1}(:,1))-min(Cent{j,1}(:,1)))/pxlsize);

end

count=1;

for j = 1:length(Cent);
    
        heigth = round((max(Cent{j,1}(:,1)) - min(Cent{j,1}(:,1)))/pxlsize);
        width  = round((max(Cent{j,1}(:,2)) - min(Cent{j,1}(:,2)))/pxlsize);
        
        rendered = hist3([Cent{j,1}(:,2),Cent{j,1}(:,1)],[heigth width]);
        
        empty  = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));
        center = round(length(empty)/2);

        empty(round(center-heigth/2):(round(center-heigth/2)) + heigth-1,round(center-width/2):(round(center-width/2)) + width-1) = rendered;        
 
name32rendered  = ['image_10nm_32bit_rendered_' num2str(j) '.tiff'];

% name16          = ['image_10nm_16bit_' num2str(i),'_',num2str(j) '.tiff'];    
% name32          = ['image_10nm_32bit_' num2str(i),'_',num2str(j) '.tiff'];

% I16 = uint16(empty);
% imwrite(I16, name16);


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

end



