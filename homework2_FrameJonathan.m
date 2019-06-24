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
% separate attributes from signatures
nSig = 13; % Can be changed later to test code perfornance
sigTable = rawTable(:,1:nSig);
attTable = rawTable(:,(nSig+1):end);
% data dimensions
[N,D] = size(rawTable);
% extract data from tables
rawData = rawTable{:,:};
sigData = sigTable{:,:};
sigNames = string(sigTable.Properties.VariableNames);
attData = attTable{:,:};
% Normalize data from the beginning, I
for di = 1:size(attData,2)
    attData(:,di)=(attData(:,di)-mean(attData(:,di)))./std(attData(:,di));
end
attNames = string(attTable.Properties.VariableNames);
% screen report
disp('------------------------------------------------------');
% screen report
fprintf('number of catchments = %d \n',N);
fprintf('number of data types = %d \n',D);

%% Homework 2. Show best MLR model combination for the 13 signatures
% SET ALPHA VALUE, FOR THE F-STATISTIC.
alpha = 0.05;
% HYPOTHESIS: Predictor significant improve our model P < alpha
% NULL HYPOTHESIS: Predictor variable doesn't significantly improve our MLR
% P >= alpha
disp('-----------------------------------------------------')
disp('The alpha value has been set for')
disp(alpha)
%------------------------------------------------------------------    
% SCRAMBLE DATA, TO GET RID OF ANY GROUPING OF CORRELATIONS
% Add one column of random numbers. Sort by those random numbers, then rm
R = [rand(N,1) rawData];
R = sortrows(R);
%Remove sorting vector
R = R(:,2:(D+1));
% Seperate the 37 Attributes
A = R(:,(nSig+1):D);
[~,nAtt] = size(A); % should be 38 = 37 + intercept
% From the 13 signatures
S = R(:,1:nSig);
%% ---------------------------------------------------------------------
% SET UP K-FOLD
% Set number of basins in each training group
n = 30;
% Split into training groups
K = floor(N/n);
disp('-----------------------------------------------------')
disp('The number of K-fold groups (K) is:'); disp(K)
%Left our due to rounding. Will be put back in later.
remainder = N - (n*K);
Nhat = N-remainder;
disp('-----------------------------------------------------')
disp('The Remainder data to be added back into K-fold groups is:')
disp(remainder)
if remainder > K
    disp('-----------------------------------------------------')
    disp('WARNING: not all data is being used.')
    disp('Adjust n value such that K >= remainder')
    disp('with 643 data points, 30 is a good number for n')    
end
% create K-fold group indecies
iKf = zeros(K,n)./0;
iKF = {}; %Cell array for intecies of different sizes (i.e., with leaveOuts
for iK = 1:K
   iKf(iK,:) = 1 + (n * (iK-1)) :  (n * (iK));
   % Add the left out data to the end of the first #leaveOuts rows.
   if iK <= remainder
       iKF{iK} = [iKf(iK,:) (Nhat+iK)];
   else
       iKF{iK} = iKf(iK,:);
   end
end
%------------------------------------------------------------------  
% Set constants used in the main loop.
regMethod = string(["All","Stepwise","PCA"]);
M = length(regMethod);
% SET ATTRIBUTES TO USE IN ANALYSIS
% 38 total for 'use-all' approach, will be narrowed if stepwise is on
% 13 & 28 are troublesome. Think about eliminating them... 
iA = [1:12, 14:37];
nB = length(iA);   
%% ------------------------------------------------------------------
% SET MEMORY FOR THE RESULTS/ANALYSIS
B = zeros(K, nB+1, nSig ,M);% Beta (corr coef) B = nB+1 for constant column
Ytest = zeros(N, nSig ,M) ./0; %prediction matrix
R2 = zeros(K, nSig ,M) ./0;  % R-squared values from the test cases
SSE = zeros(K, nSig ,M) ./0; % Sum of squared errors
SST = zeros(K, nSig ,M) ./0; % Sum of total squared errors
varsSW = zeros(K,nSig);% Number of variable used in SW Linear Regression
varsPCA = zeros(K,nSig);% Number of components used in PCA
Xvtest = zeros(N, nB, nSig); % transformed Attributes for the PCA analysis
%% ------------------------------------------------------------------
% MODEL LOOP. Testing three regression techniques.
for m = 1:M
    %------------------------------------------------------------------
    % K-FOLD LOOP. This loops through the k-fold groups.
    % Each loop should train the model on K-1 groups, and test on 1 group
    for iK = 1:K % K is the left out K-Fold group for testing
        % SET INDICIES FOR THE TRAINING AND TESTING BASED ON THE K-FOLD
        % Clear the regression index
        iTest = [];
        iTrain = [];    
        for iNotK = 1:K
            if iNotK ~= iK
                % Index for this regression calculation, Training group
                iTrain = [iTrain iKF{iNotK}];
            else
                iTest = [iTest iKF{iK}];
            end
        end 

        % Index of signature to predict with our multi-linear regression model
        for iSig = 1:nSig
            % Dependent variables
            Ytrain = S(iTrain,iSig);
            % Independent variable for training.
            Xtrain = A(iTrain,iA);
            Otrain = ones(size(Xtrain,1),1);
            % Independent variable for testing.
            Xtest = A(iTest,iA);
            Otest = ones(size(Xtest,1),1);
            %--------------------------------------------------------------   
            % TRAINING SEGMENT
            %--------------------------------------------------------------
            % USE ALL DIMENSIONS
            if strcmp(regMethod(m),'All')
                % Train the model on the majority of the data
                Xtrain = [Otrain,Xtrain];
                [b,bint,r,rint,stats] = regress(Ytrain, Xtrain);
                % Store the correlation coefficients, for lumped results
                B(iK,:,iSig, m) = b;
                % Set the regression model coefficients for...
                % the left out K-fold group (iTest).
                Btest = B(iK,:,iSig,m)';
                Xtest = [Otest,Xtest];
            %--------------------------------------------------------------
            % USE STEPWISE LINEAR REGRESSION
            % stepwiselm!!!!!!!!!!! stepwiselm.Formulas.Terms
            elseif strcmp(regMethod(m),'Stepwise')
                SW = stepwiselm(Xtrain,Ytrain,...
                    'PEnter',alpha,'Upper','linear');
                iC=find(sum(SW.Formula.Terms)==1);% Find the selected vars 
                disp(iC);
                % insert the number of variables. Without intercept
                varsSW(iK,iSig) = length(iC); 
                % Calculate the coefficients with only the good variables.
                Btest = SW.Coefficients.Estimate; 
                %Btest = regress(Ytrain,[Otrain,Xtrain(:,iC)]);
                % add in the constant column, add one to rest of
                iC = [1, iC+1];
                B(iK,iC,iSig,m) = Btest';
                % spread the coefficients across all the whole attribute
                Btest = B(iK,:,iSig,m)';
                disp('K-fold group'); disp(iK)
                disp('Signature'); disp(iSig)
                disp(SW);
                % Add constants to testing group
                Xtest = [Otest,Xtest];
            %--------------------------------------------------------------
            % USE PRINCIPAL COMPONENT ANALYSIS TO REDUCE THE DIMENSIONS
            elseif strcmp(regMethod(m),'PCA')
                % Set minimum number of principal components
                % because i don't want to use just one component
                % might as well use as many components as SWLR...
                MinPC = varsSW(iK,iSig);
                [v,~,lamda,tsq,expl,mu] = pca(Xtrain,'Economy',false);
                fcum = cumsum(expl);
                iC = max(MinPC,find(fcum > (1-alpha*100), 1, 'first'));
                varsPCA(iK,iSig) = iC; % store the number of components.
                % Correlation coefficients for active components
                Xv = Xtrain * v(:,1:iC);   %Xv(:,1:iC);
                % Assign coefficients
                B(iK,1:iC+1,iSig, m) = regress(Ytrain, [Otrain,Xv]);
                % Btest used in K-fold stats
                Btest = B(iK,:,iSig, m)';
                % Project the testing data on the principal component space
                Xvtest(iTest,1:iC,iSig) = Xtest * v(:,1:iC);
                Xtest = [Otest,Xvtest(iTest,:,iSig)];
            end
            %---------------------------------------------------------------  
            % TESTING SEGMENT        
            %--------------------------------------------------------------
            % Test the regression model on the left out K-fold group(iTest)
            Ytest(iTest,iSig,m) = Xtest * Btest;
            % Calculate the SSE between actual and MLR model   
            Ydata = S(iTest,iSig);
            % Calculate and store the R^2 values
            Ydiff = Ydata - Ytest(iTest,iSig,m);
            SSE(iK, iSig, m) = sum((Ydiff).^2);
            SST(iK, iSig, m) = sum((Ydata - mean(Ydata)).^2);
            R2(iK,iSig,m) = 1 - ( SSE(iK, iSig, m) / SST(iK, iSig, m) );
        end % Signature loop
    end % K-fold loop
end % model loop
disp(R2);
%% --- PLOTTING RESULTS ---------------------------------------------------
% This seems really bad to me, but the box plot wants a vector of all R2s
R2m1 = [];
R2m2 = [];
R2m3 = [];
for iSig = 1:13
    R2m1 = [R2m1;R2(:,iSig,1)];
    R2m2 = [R2m2;R2(:,iSig,2)];
    R2m3 = [R2m3;R2(:,iSig,3)];
end
R2plot = [R2m1; R2m2; R2m3];
% This big dumb loop is assigning the signature to the R2 vector
c13 = 1;
c21 = 1;
sigs = zeros(3*13*21,1);
for c = 1:3*13*21
    sigs(c) = c13;
    if c21 == 21
        if c13 == 13
            c13 = 0;
        end
        c13 = c13+1;
        c21 = 0;
    end
    c21 = c21+1;
end
re = 13*21;
models = [repmat({'All'},1,re), repmat({'SW'},1,re) repmat({'PCA'},1,re)];
fignum = fignum+1;
figure(fignum); close(fignum); figure(fignum);  
width=15;
height=6.6;
set(gcf,'units','inches','position',[10,10,width,height])
boxplot(R2plot,{sigs, models},'colors', ['k','r','b']);
grid on;
ylim([-1,1]);
%% --- TIME TO COMPLETION -------------------------------------------------
t = toc;
fprintf('. finished - time = %f [seconds] \n',t);
%% --- END ---------------------------------------------------------
