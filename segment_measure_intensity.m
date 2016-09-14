clear, close all, clc
%% Load image

WF_name = 'NB_A647_200mW_filterUp_200ms.tif'

I=imread(WF_name);

% I2 = im2uint8(I);
% I=I2;

fprintf('\n -- Data loaded --\n')

%% Adjust segmentation parameters
close all

% blur the image
figure('Position',[10 600 500 500],'name','Raw Image'), imshow(I,'InitialMagnification','fit');

% adjust the contrast of the raw image
I2 = imadjust(I,[0.02 0.09],[]);
figure('Position',[600 600 500 500],'name','Image after adjusted contrast'), imshow(I2,'InitialMagnification','fit');

G = fspecial('gaussian',[7 7],10); % lowpass filter of size and gaussian blur sigma, [lowpass filter] sigma
imG = imfilter(I2,G,'same');
figure('Position',[1200 600 500 500],'name','Image after adjusted contrast, blurring'), imshow(imG,'InitialMagnification','fit');

% adjust the background
I3 = imadjust(imG,[0.2 0.9],[]);
figure('Position',[10 10 500 500],'name','Image after adjusted contrast, blurring, adjusted contrast'), imshow(I3,'InitialMagnification','fit');

% Make binary image
bin = im2bw(I2,0.3);
figure('Position',[600 10 500 500],'name','Binary image result'),imshow(bin,'InitialMagnification','fit')
[B,L,N,A] = bwboundaries(bin); % B - connectivity

%% Extract particles

close all

% Extract the integrated intensity of the GFP WF image for each ROI

intI = [];

for i = 1:length(B);
    
    intI(i,1) = sum(sum(I(min(B{i,1}(:,1)):max(B{i,1}(:,1)),min(B{i,1}(:,2)):max(B{i,1}(:,2)))));
    
end

% Intensity Histogram 

bins = 1:50:1e3;
h = hist(intI,bins);
bar(bins, h);

save('Sas6_5.mat','intI');

