m =30;

dataDBS      = [];
dataDBS(:,1) = Cent{m,1}(:,1); % x in mum
dataDBS(:,2) = Cent{m,1}(:,2); % y in mum

% Run DBSCAN on each particle 

k   = 10;                                                 % minimum number of neighbors within Eps
Eps = 15;                                                 % minimum distance between points, nm

fprintf('\n -- DBSCAN input and parameters selected --\n')

tic
[class,type]=DBSCAN(dataDBS,k,Eps);     % uses parameters specified at input
class2=transpose(class);                % class - vector specifying assignment of the i-th object to certain cluster (m,1)
type2=transpose(type);                  % (core: 1, border: 0, outlier: -1)

fprintf(' -- DBSCAN computed in %f sec -- \n',toc)

coreBorder = [];
coreBorder = find(type2 >= 0);

subsetP = [];
subsetP(:,1)=dataDBS(coreBorder,1);
subsetP(:,2)=dataDBS(coreBorder,2);
subsetP(:,3)=class2(coreBorder);

figure('Position',[700 600 900 400])
subplot(1,2,1)
scatter(dataDBS(:,1),dataDBS(:,2),1);
title('Raw Data')
axis on
axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
box on

subplot(1,2,2)
scatter(subsetP(:,1),subsetP(:,2),1,mod(subsetP(:,3),10))
title('identified Clusters')
axis on
axis([min(dataDBS(:,1)) max(dataDBS(:,1)) min(dataDBS(:,2)) max(dataDBS(:,2))])
box on

