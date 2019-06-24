%% --- initialize workspace -----------------------------------------------------

clear all; close all; clc
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- metadata setup -----------------------------------------------------

% data directories
dataDir = 'extracted_data';  % where to find the necessary data

%% --- load data ----------------------------------------------------------

% Load water mask
IGBP = load('NLDAS_IGBPpredomveg.asc');


km2 = IGBP(:,5);
watermask = IGBP(:,22);
watermask(watermask > km2/2) = 0/0;
watermask(watermask <= km2/2) = 1;

% screen report
tic;
fprintf('Loading data ...');

% load the data into a table format
fname = strcat(dataDir,'/all_data.txt');
rawTable = readtable(fname);

% remove catchments with missing values
iMissing = find(all(~isnan(rawTable{:,:}')));
rawTable = rawTable(iMissing,:);

%% --- train model with Aridity---------------------------------------

% number of k-folds
Nkfold = 5;
fprintf('.............................................................. \r')
fprintf('Training Spatial Model plus Aridity ... \r');tic

predictorNames = {'gauge_lat', 'gauge_lon', 'aridity'};
responseNames = 'q_mean';

trainingData = rawTable;
trainingData{:,1} = log(trainingData{:,1});

[trainedArid, AridRMSE] = ...
    jf_trainRegressionModel(trainingData, predictorNames, responseNames,...
    Nkfold);

fprintf('Finished: time = %f[s] \n \r');toc
fprintf('Root mean squared error: \r');
disp(AridRMSE)
%% --- train model just spatial  --------------------------------------

% number of k-folds
Nkfold = 5;

fprintf('Training Spatial Model ... \r');tic

predictorNames = {'gauge_lat', 'gauge_lon'};
responseNames = 'q_mean';

trainingData = rawTable;
trainingData{:,1} = log(trainingData{:,1});

[trainedModel, validationRMSE] = ...
    jf_trainRegressionModel(trainingData, predictorNames, responseNames,...
    Nkfold);

fprintf('Finished: time =%f[s] \n');toc
fprintf('Root mean squared error: \r');
disp(validationRMSE);
%% -- make predictions

% make a spatial grid
Wlon = -124.9375;
Elon = -67.0625;
Nlat = 52.9375;
Slat = 25.0625;

%grid resolution
delta = 1/8;

% define the lat/lon coordinates
latdex = Slat:delta:Nlat;
londex = Wlon:delta:Elon;

%dimensions
Dlat = length(latdex);
Dlon = length(londex);

%lat/lon lists --> lat/lon grid
[latgrid, longrid] = meshgrid (latdex, londex); 

% make gpr predictions

predictTable = table(latgrid(:), longrid(:),'VariableNames',{'gauge_lat','gauge_lon'});
Ylist = trainedModel.predictFcn(predictTable);
Ylist = Ylist .* watermask;

% turn into spatial field
Ygrid = reshape(Ylist, [Dlon, Dlat]);

%% ---plot spatial predictions

%initialze figure
fig = 1; figure(fig); close(fig); figure(fig);

%plot
% Ygrid = log(Ygrid);
surface(longrid, latgrid, Ygrid, 'edgecolor','none');
title('Kriging (Gaussian process) interpolation of 643 Catchements)',...
    'fontsize', 18);
xlabel('longitude [deg]','fontsize', 18);
ylabel('latitude [deg]','fontsize', 18);
hold on;

% plot gauge locations
[X,Y] = meshgrid(rawTable{:,'gauge_lon'}, rawTable{:,'gauge_lat'});
Z = repmat(max(Ylist), length(X));
plot3(rawTable{:,'gauge_lon'}, rawTable{:, 'gauge_lat'},Z,'ok')

cb = colorbar;
cbtl = cb.TickLabels;
for i = 1:length(cbtl)
    adjval= exp(str2num(cbtl{i}));
    cbtl{i}= num2str(round(adjval*100)/100);
end
cb.TickLabels = cbtl;
cb.Label.String = 'mean surface runoff log(mm/d)';
cb.Label.FontSize=18;

colormap('bone')

%% end script
