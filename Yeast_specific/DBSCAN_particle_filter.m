%% 
clear, clc, close all

%% 

for m=1:length(Cent);

dataDBS = [];

dataDBS(:,1) = Cent{m,1}(:,1); % x in mum
dataDBS(:,2) = Cent{m,1}(:,2); % y in mum

scatter(dataDBS(:,1),dataDBS(:,2),1)

% Run DBSCAN on each centriole 

k   = 20;                                                   % minimum number of neighbors within Eps
Eps = 30;                                                 % minimum distance between points, nm

fprintf('\n -- Parameters selected --\n')

tic
[class,type]=DBSCAN(dataDBS,k,Eps);     % uses parameters specified at input
class2=transpose(class);                % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);                  % (core: 1, border: 0, outlier: -1)

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

coreBorder = [];
coreBorder = find(type2 >= 0);

subset          = [];
subset          = Cent{m,1}(coreBorder,1:end);
subset(:,end+1) = class2(coreBorder);


% figure('Position',[700 600 900 400])
% subplot(1,2,1)
% scatter(dataDBS(:,1),dataDBS(:,2),1);
% title('Raw Data')
% axis on
% axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
% box on
% 
% subplot(1,2,2)
% scatter(subset(:,1),subset(:,2),1,mod(subset(:,3),10))
% title('identified Clusters')
% axis on
% axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
% box on

% Select only large cluster(s)

for i=1:max(subset(:,end));
    
    vx = find(subset(:,end)==i);
    
    if length(vx)>800;
        
         Cent{m,4}{i,1} = subset(vx,1:end);
         
    end
end

end

save('Extracted_Particles_5.mat','Cent');

close all;
%% 

Cent = Cent8.Cent;

% figure
for i=1:length(Cent);
    
    if isempty(Cent{i,4});
    
    elseif Cent{i,4}>10000;
    
    else
        
scatter(length(Cent{i}),Cent{i,4},5,'filled','r');hold on;
xlabel('Nbr of localizations');
ylabel('GFP integrated intensity');
box on;
axis square;

    end
end

