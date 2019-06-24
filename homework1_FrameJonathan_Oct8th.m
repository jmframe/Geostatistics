% The paper that describes this data set is on EarthArXiv:
%   https://eartharxiv.org/2em53/

%% --- init workspace -----------------------------------------------------

clear all; close all; clc
restoredefaultpath; 
fignum = 0;

%% --- metadata setup -----------------------------------------------------

% data directories
% where to find the necessary data
dataDir = 'extracted_data/';
%   Figure directory
fileDir = '/Users/jf/box sync/class/geostat/hw1/write/';
% Save figure files if == 1
saveFile = 0;

%% --- load data ----------------------------------------------------------

% screen report
tic;
fprintf('Loading data ...');

% load the data into a table format
fname = strcat(dataDir,'/all_data.txt');
rawTable = readtable(fname);

% remove catchments with missing values
iMissing = find(all(~isnan(rawTable{:,:}')));
rawTable = rawTable(iMissing,:);

% extract some key data columns - you may need more 
qmean    = rawTable{:,'q_mean'};            % mean flow in each catchment
runRat   = rawTable{:,'runoff_ratio'};      % mean flow in each catchment

% separate attributes from signatures
attTable = rawTable(:,14:end);
sigTable = rawTable(:,1:13);

% data dimensions
[N,D] = size(rawTable);

% extract data from tables
rawData = rawTable{:,:};
sigData = sigTable{:,:};
attData = attTable{:,:};

% screen report
t = toc;
fprintf('. finished - time = %f [seconds] \n',t);
disp('------------------------------------------------------');
% screen report
fprintf('number of catchments = %d \n',N);
fprintf('number of data types = %d \n',D);

%% Homework Question: b)  Calculate the mean, standard deviation,......
%  skewness and kurtosis of all variables involved in these three corrs.
%   Generating Statistics for all variables involved in these three 
%   correlations: "p_mean", "q_mean", "low_prec_freq",...
%                 "high_prec_freq", "runoff_ratio"
hiCorVars = {'q mean','p mean', 'runoff ratio', 'low prec freq', 'q 95'};
%   Create table with only the variables in the three correlations.
hiCorMat = [rawTable.q_mean, rawTable.p_mean, rawTable.runoff_ratio,...
    rawTable.low_prec_freq, rawTable.q95];
%   Calculate the stats for the high correlation variables.
hiCorStats = [mean(hiCorMat); std(hiCorMat); skewness(hiCorMat);...
    kurtosis(hiCorMat)];
%   Name the columns and rows in the statistics table
varNames = ["qmean", "pmean",  "runoff_ratio", "low_prec_freq", "q95"];
disp('------------------------------------------------------');
disp(varNames);
disp(hiCorStats);
V = size(hiCorVars,2);
%% HOMEWORK Question c)  Plot the histograms of all variables involved...

% in these correlations.
%   create synthetic data with statistical values
%    just to compare with the empiracal data
%    not part of the aissignment
hiCorMatSynth = zeros(1000000,5);
for i = 1:V
    hiCorMatSynth(:,i) = pearsrnd(hiCorStats(1,i), hiCorStats(2,i),...
        hiCorStats(3,i), hiCorStats(4,i),1000000,1);
end
%   Now plot Normalized histograms
%    for fun add theoretical distributions on top of histogram.
for i = 1:V
    fignum = fignum+1;
    figure(fignum); close(fignum); figure(fignum);
   
    %Define the edges of the histogram bins.
    plotEdges = floor(min(hiCorMat(:,i))):(ceil(max(hiCorMat(:,i))) - ...
        floor(min(hiCorMat(:,i)))) / 25:max(hiCorMat(:,i));
      
    %Create and plot the theoretical normal distribuion
     %    just for fun
    histogram(hiCorMatSynth(:,i),plotEdges,...
        'normalization','probability', 'LineStyle', 'none',...
        'FaceColor','r','facealpha',0.4);
    set(gcf,'units','centimeters','position',[10,10,7,7])
    hold on; %   plot synthetics data behind, and the real data in front
    
    %   Add normalize tags to these histograms.
    histogram(hiCorMat(:,i),plotEdges,'normalization','probability',...
        'FaceColor','b','facealpha',0.8);
    xlabel(hiCorVars(i),'fontsize',11);
    ylabel('probability of occurance','fontsize',11);    
    %legend('Synthetic data','Empirical','fontsize',11,'Location',...
    %    'southoutside');
    hold off;
    if saveFile == 1
        figFile = strcat(fileDir,'_','Hist_',...
            varNames(i));
        saveas(figure(fignum),figFile, 'png');
    end
end

%% Homework question Plot heat maps or contour plots of the joint dist.----
%in all three of these correlations.
%   Code adapted from https://stackoverflow.com/questions/16313949
%   Correlation 5:
%   high_prec_freq vs runoff_ratio, cor=0.69141
%   Correlation 4: 
%   p_mean vs Runoff_ratio, cor=0.7175
%   Correlation 3: 
%   low_prec_freq = hiCorMat(4) vs runoff_ratio = hiCorMat(3), cor:0.73791
%   Correlation 2:
%   q_95  = hiCorMat(5)    vs p_mean = hiCorMat(2), cor: 0.845
%   Correlation 1: 
%   p_mean hiCorMat(2), q_mean hiCorMat(1), cor: 0.88689
%   Correlation Pairs
corPairs = [2,1;2,5;4,3];
%   loop through and plot contour maps for each correlation pair
for i = 1:3
    i1 = corPairs(i,1);
    i2 = corPairs(i,2);
    xi = linspace(min(hiCorMat(:,i1)), max(hiCorMat(:,i1)), 20);
    yi = linspace(min(hiCorMat(:,i2)), max(hiCorMat(:,i2)), 20);
    hst = hist3([hiCorMat(:,i1),hiCorMat(:,i2)],{xi yi});
    % normalize the histogram data
    dx = xi(2)-xi(1);
    dy = yi(2)-yi(1);
    area = dx*dy;
    pdfData = hst/sum(sum(hst))/area;
    
    % plot pdf
    fignum = fignum+1;
    figure(fignum); close(fignum); figure(fignum);
    contour(xi,yi,pdfData);
    set(gcf,'units','centimeters','position',[10,10,8,8])
    xlabel(hiCorVars(i1),'fontsize',11);
    ylabel(hiCorVars(i2),'fontsize',11);
    colorbar;
    if saveFile == 1
        figFile = strcat(fileDir,'_Joint_',...
            varNames(i1),'_vs_',hiCorVars(i2));
        saveas(figure(fignum),figFile, 'png');
    end
    
end

%% Use random sampling to sample 50 catchments from the list and 
% calculate the sample mean and 90% confidence intervals of the sample mean 
% for the variable <runoff_ratio>. Plot the distribution of the sample mean 
% against the actual mean from all 671 catchments.
% indexes for random sample
d = runRat;
K = 50;
N = length(d);
%INPUTS: DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K: NUMBER OF RANDOM SAMPLES
iRandom = randsample(N,K);
%   get the statistics for all the runoff ratios
dStats = [mean(d), std(d), skewness(d), kurtosis(d)];
%   grab the runoff ratio of the 50/100 random catchments
dRan = d(iRandom);
dRanStats = [mean(dRan), std(dRan),skewness(dRan),...
    kurtosis(dRan), var(dRan)];
disp('Statistics of random sample...');
disp('mean, standard deviation, skewness, kurtosis');
format short
disp(dRanStats);

%   Plot the probability distributions of both the sample and total 
%   code adapted from https://stackoverflow.com/questions/16133703/
%   Data Plotting range
x_d = (-5 * dStats(2):0.01:5 * dStats(2)) + dStats(1);  
%   Data Probability distribution function
y_d = exp(- 0.5 * ((x_d - dStats(1)) / dStats(2)) .^ 2)/...
    (dStats(2) * sqrt(2 * pi));
y_d_n = y_d/sum(y_d);

%   Random sampled Data Plotting range
x_dRan = (-5 * dRanStats(2):0.01:5 * dRanStats(2)) +...
    dRanStats(1);  
%   Random Sample Probability distribution function
y_dRan = exp(- 0.5 * ((x_dRan - dRanStats(1)) /...
    dRanStats(2)) .^ 2) / (dRanStats(2) * sqrt(2 * pi));
y_dRan_n = y_dRan/sum(y_dRan);

%   plot pdf
fignum = fignum+1;
figure(fignum); close(fignum); figure(fignum);
%   plot random sample
%   plot total record
plot(x_d, y_d_n, 'k', 'LineWidth', 3);
hold on;
plot(x_dRan, y_dRan_n, 'r', 'LineWidth', 2);
grid on;
set(gcf,'units','centimeters','position',[10,10,15,15])
xlabel('Runoff Ratio','fontsize',11);
ylabel(strcat('p(Runoff Ratio)'),'fontsize',11);
legend('All catchments',strcat(num2str(K),' random catchments'),...
    'fontsize',11,'Location','southoutside');
hold off;
if saveFile == 1
    figFile = strcat(fileDir,'_',signature,num2str(K),'_RandSampleDist');
    saveas(figure(fignum), figFile, 'png');
end

SEM = std(dRan)/sqrt(length(dRan));               % Standard Error
z = tinv([0.05  0.95],length(dRan)-1);          % T-Score
CIr = mean(dRan) + z * SEM;                      % Confidence Intervals
disp(strcat('The 90 percent confidence interval of the mean'));
disp(strcat('with K =',{' '},num2str(K),' samples is =',{' '},...
    num2str(CIr)));
disp('------------------------------------------------------');   


%%  stratified sampling by geo or veg class
% CHOOSE THE ATTRIBUTE OF WHICH TO STRATIFY THE RUNOFF RATIO

for stratAtt = [1 35] % Geology = 1, Vegetation = 35
    stratData = attData(:,stratAtt);
    if stratAtt == 1
        Attribute = 'Geological';
    elseif stratAtt == 35
        Attribute = 'Vegetation';
    end

    N = length(runRat);
    u = unique(stratData);  Nstrata = length(u);
    %   define stratified statistics matrix
    %Stats: stratum, n, mean, stdev, skewness, kurtosis
    stratStat = zeros(Nstrata,5)./0;
    % initialize a vector to store the length of each stratum
    Nstratum = zeros(Nstrata,1) /.0;
    % initialize a vector of standard deviation and mean of each stratum
    stdStrat = zeros(Nstrata,1) /.0;
    meanStrat = zeros(Nstrata,1) /.0;    
    % product of stratum sample mean and total stratum size, for use in
    % estimation of the mean
    Nkzhat = zeros(Nstrata,1) /.0; 

    %  Calculate length and the standard deviation for each stratum
    for i = 1:Nstrata
        % index of the stratification values, don't store value, since
        % it will change with each strata. Or could have used {}
        iStrat = find(stratData == u(i));
        Nstratum(i) = length(iStrat);
        stdStrat(i) = std(runRat(iStrat));
        meanStrat(i) = mean(runRat(iStrat));        
    end
    % Calculate the sum of the product between the stratum N & std
    sumSknk = dot(Nstratum, stdStrat);

    %  calulate strat-sampling statististics
    %  this loop generates a table of statistics for the strata
    for i = 1:Nstrata
        %  index of the stratification values, don't store value, since
        %   it will change with each strata
        iStrat = find(stratData == u(i));
        %  Choose optimal stratum sample
        nOPT = round(K * (Nstratum(i) * stdStrat(i)) / sumSknk);
        % THRESHOLD for using a sample size
        if nOPT < 4
            nOPT = 0;
        end
        %  generating a random sample within the strata
        %  CAN STORE THIS WITH THE {}
        iStrat = randsample(iStrat,nOPT); 
        %  collect the values from our stratum sample
        stratSample = runRat(iStrat);
        %   CALCULATE THE SAMPLE MEAN
        stratSampleMean =  mean(stratSample);
        %   Calculate the sum-product of the sample mean & stratum size
        Nkzhat(i) = stratSampleMean * Nstratum(i);
        %Calculate the within-stratum sample variance
        sk2 = 1/(nOPT - 1) * sum( ( stratSample - meanStrat(i) ).^2 );
        %stratified sample statistics
        stratStat(i,:) = [u(i), length(stratSample),...
            mean(stratSample), sk2, stratSampleMean];
    end
    totalStratSamples = sum(stratStat(:,2));
    disp('------------------------------------------------------');
    disp(strcat(Attribute, ' stratification, samples = ', num2str(K)));
    disp({'stratum','n Samples','Sample mean','Sample Var'});
    disp(stratStat);
    disp('Total samples through all strata =');
    disp(totalStratSamples);    

    % Index the true value variances
    iNotNAN = [];
    for i = 1:Nstrata
        if stratStat(i,4) > 0
            iNotNAN = [iNotNAN,i];
        end    
    end
    disp('Strata with true value variances');
    disp(iNotNAN)

    %   Stratified Sample Mean
    sumSkzk = sum(Nkzhat(iNotNAN));
    stratMean = (1/N) * sumSkzk;
    disp('Estimate of the mean =');
    disp(stratMean);

    stratVar = 0;
    for i = iNotNAN
        Nk = Nstratum(i,1);
        nk = stratStat(i,2);
        sk2 = stratStat(i,4);
        stratVar = stratVar + (sk2/nk) * (Nk/N)^2 * (1-nk/Nk);
    end
    stratStd = sqrt(stratVar);
    disp('Variance of the mean =');
    format short
    disp(round(stratVar,6));

    % calcualte the relative precisions of strat sampling
    %   Relative precision = S^2(z_hat) / (S_s)^2(z_hat)
    relativePrecision = dRanStats(5) / stratVar;

    disp('------------------------------------------------------');
    disp(strcat(Attribute, {' '} ,'stratified relative precision',...
        ' with',{' '},num2str(totalStratSamples),' samples is =',...
        num2str(relativePrecision)));


    % plot stratified mean distributions
    %   stratified Plotting range
    x_strat = (-5 * stratStd:0.01:5 * stratStd ) +  stratMean;  
    %   stratified Probability distribution function
    y_strat =  exp(- 0.5 * ((x_strat - stratMean) /...
        stratStd) .^ 2) / (stratStd * sqrt(2 * pi));
    y_nStrat = y_strat / sum(y_strat);

    % HW Question----------------------------------------------------------
    % Calculate the 90% confidence intervals of the sample mean
    % Code adapted from https://www.mathworks.com/matlabcentral/answers/...
    % 159417-how-to-calculate-the-confidence-interval
    SEM = stratStd/sqrt(totalStratSamples);               % Standard Error
    z = tinv([0.05  0.95],totalStratSamples-1);      % T-Score
    CIs = stratMean + z * SEM;    % Confidence Intervals
    disp(strcat('The 90 percent confidence intervals of the',{' '},...
        Attribute,' stratified sample mean'));
    disp(strcat('with ',{' '}, num2str(totalStratSamples),{' '},...
        'samples are',{' '},num2str(CIs)));

    if strcmp(Attribute, 'Geological')
        GEO_CIs = CIs;
        GEO_x_strat = x_strat;
        GEO_y_nStrat = y_nStrat;
    elseif strcmp(Attribute, 'Vegetation')
        VEG_CIs = CIs;
        VEG_x_strat = x_strat;
        VEG_y_nStrat = y_nStrat;
    end
end

%% PLOTTING THE MEAN FROM THE STRATIFIED and Random SAMPLES
%   plot pdf
fignum = fignum+1;
figure(fignum); close(fignum); figure(fignum);
%   Plot the probability distributions of the random sample, 
%    stratified samples and total 
%   plot total record first
plot(x_d,y_d_n, 'k', 'LineWidth',3);
hold on;
plot(GEO_x_strat, GEO_y_nStrat, 'b', 'LineWidth',3);
plot(VEG_x_strat, VEG_y_nStrat, 'r', 'LineWidth',3);
set(gcf,'units','centimeters','position',[10,10,15,15])
xlabel('Runoff Ratio','fontsize',11);
ylabel('p(Runoff Ratio)','fontsize',11);
%   plot stratified samples
xlim([floor(min(min(x_strat))), ceil(max(max(x_strat)))]);
legend('All catchments','Geo stratified samples',...
    'Veg stratified samples',...
    'fontsize',11,'Location','southoutside');
if saveFile == 1
    figFile = strcat(fileDir,'_','RunoffRatio_',...
        num2str(K),'StratSampleDist');
    saveas(figure(fignum), figFile, 'png');
end

hold off;
disp('------------------------------------------------------');


%% -- test for normality -----------------------------------------------

% test raw runoff ratio data with k-s test
%   Test runoff ratios for normality.
for Trans = ["None", "BoxCox", "Log", "Log_Alpha"]

    if strcmp(Trans, 'None')
        disp('----------------------------------------------------------');
        d = runRat;
        T = 'Runoff ratios with no transformation';
        disp(T);
    elseif strcmp(Trans, 'BoxCox')
        disp('----------------------------------------------------------');
        d = boxcox(runRat);
        T = 'Runoff ratios with Box Cox transformation';
        disp(T);
    elseif strcmp(Trans, 'Log')
        disp('----------------------------------------------------------');
        T = 'Runoff ratios with Log transformation';
        d = log(runRat);
        disp(T);
    elseif strcmp(Trans, 'Log_Alpha')
        disp('----------------------------------------------------------');
        T = 'Runoff ratios with Log transformation and Alpha shift';
        d = log(runRat + 1.2);
        disp(T);        
    end

    m = mean(d);
    s = std(d);
    ds = (d-m)/s;
    [h,p] = kstest(d);      %Kolmagor-Schmirnof Test
    [hs,ps] = kstest(ds); % WITH STANDARDIZATION
    if h == 1
        hyp = ' rejected';
    else
        hyp = ' not rejected';
    end
    disp('------------------------------------------------------');    
    disp(strcat('Testing for normal distribution hypothesis is ',...
        hyp, ', p value = ',{' '}, num2str(p)));
    if hs == 1
        hyp_s = ' rejected';
    else
        hyp_s = ' not rejected';
    end
    disp('------------------------------------------------------');    
    disp(strcat('Testing for normal distribution hypothesis is ',...
        hyp_s, ', p value = ',{' '}, num2str(ps)));

    % plot the emirical and theoretical cdfs
    y = randn(10000,1);
    [cdf1a,cdf1b] = ecdf(d);
    [cdf2a,cdf2b] = ecdf(ds);
    [cdf3a,cdf3b] = ecdf(y);
    %   plot CDFs with subtracted mean and divided by std
    fignum = fignum+1;
    figure(fignum); close(fignum); figure(fignum);
    plot(cdf1b,cdf1a, 'k', 'LineWidth',3);
    title(T);
    grid on;
    hold on;
    plot(cdf2b,cdf2a, 'b', 'LineWidth',3);
    plot(cdf3b,cdf3a, ':r', 'LineWidth',3);
    xlabel('Runoff Ratio','fontsize',11);
    ylabel('Cumulative probability','fontsize',11);
    legend('Empirical CDF','Empirical CDF Standardized',...
        'Theoretical CDF','Location','southoutside','fontsize',11);
    set(gcf,'units','centimeters','position',[10,10,15,15])    
    hold off;
end
%% --- hypothesis testing -------------------------------------------------

% test whether evergreen forests have different runoff ratios...
%than other land cover types
%FROM HOMEWORK: a)  Is there a statistically significant difference...
%   between runoff ratios at forested sites vs. sites with other...
vegClass = attData(:,35);
%   land cover types? Why do you think this might be?
%FOREST CLASSES = 1, 2, 11 & 12
iForest = zeros(1,1) ./0;
iNotForest = zeros(1,1) ./0;
for i = 1:length(vegClass)
    if vegClass(i) == 1 || vegClass(i) == 2 ||...
            vegClass(i) == 11 || vegClass(i) == 12 
        iForest = [iForest; i];
    else
        iNotForest = [iNotForest; i];
    end
end
iForest(1,:) = [];
iNotForest(1,:) = [];
iForest = sortrows(iForest);
iNotForest = sortrows(iNotForest);
rrForest = runRat(iForest);
rrNotForest = runRat(iNotForest);
% test whether the collecation of all forest-type land covers have a
% different mean
[h_Forest,p_Forest] = ttest2(rrForest,rrNotForest,'Vartype','unequal');
if h_Forest == 1
    hyp_For = ' rejected';
else
    hyp_For = ' not rejected';
end
disp('------------------------------------------------------');
fprintf(strcat('Forest Stratification hypothesis is ',...
    hyp_For, ', p value = ',...
    num2str(p_Forest),'\n'));


% test whether all sedimentary geo types have a different mean
%FROM HOMEWORK: b)  Is there a statistically significant difference...
%   between runoff ratios at sites where the primary geological class...
%   is sedimentary vs. other classes? What about sites with...
%   unconsolidated sediments? Why do you think this might be?
%SEDIMENTARY CLASSES = 1, 4, 6, 7 & 8
geoClass = attData(:,1);
iSed = zeros(1,1) ./0;
iNotSed = zeros(1,1) ./0;
for i = 1:length(geoClass)
    if geoClass(i) == 1 || geoClass(i) == 4 ||...
            geoClass(i) == 6 || geoClass(i) == 7 || geoClass(i) == 8
        iSed = [iSed; i];
    else
        iNotSed = [iNotSed; i];
    end
end
iSed(1,:) = [];
iNotSed(1,:) = [];
iSed = sortrows(iSed);
iNotSed = sortrows(iNotSed);
rrSed = runRat(iSed);
rrNotSed = runRat(iNotSed);
[h_Sed,p_Sed] = ttest2(rrSed,rrNotSed,'Vartype','unequal');
if h_Sed == 1
    hyp_Sed = ' rejected';
else
    hyp_Sed = ' not rejected';
end
disp('------------------------------------------------------');
fprintf(strcat('Geological stratification hypothesis is ',...
    hyp_Sed, ', p value = ',...
    num2str(p_Sed),'\n'));


%   Unconsolidated sedimentary only (Geological class = 7)
iSed_u = find(geoClass() == 7);
iNotSed_u = find( geoClass() ~= 7);

iSed_u(1,:) = [];
iNotSed_u(1,:) = [];
iSed_u = sortrows(iSed_u);
iNotSed_u = sortrows(iNotSed_u);
rrSed_u = runRat(iSed_u);
rrNotSed_u = runRat(iNotSed_u);
[h_Sed_u,p_Sed_u] = ttest2(rrSed_u,rrNotSed_u,'Vartype','unequal');
if h_Sed_u == 1
    hyp_Sed_u = ' rejected';
else
    hyp_Sed_u = ' not rejected';
end
disp('------------------------------------------------------');
fprintf(strcat('Unconsolidated sedimentary hypothesis is ',...
    hyp_Sed_u, ', p value = ',...
    num2str(p_Sed_u),'\n'));

%% --- END ---------------------------------------------------------
