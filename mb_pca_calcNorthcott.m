%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code to perform PCA and LDA and generate Figures 2B-E for 
% ***Fitzgerald et al. (2023)***
% If you make use of this code pleasse cite the above paper
% Thanks to Cristiano Gutta for sharing the original PCA and LDA code. 
% 1. Import individual MB cellline protein expressions
% 2. Normalise expression values
% 3. Perform PCA on measured protein values and plot (incl. coefficient
%    contributions)
% 4. Perform LDA on measured protein values and plot
% 5. Import patient data, normalise and plot on previous PCA and LDA plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mb_pca_calcNorthcott
%% 1. Import individual protein values
mb_proteins_table = readtable('MB_PCA_input_individual_6.xlsx','ReadRowNames',true);
% Proteins:
% 1=APAF1; 2=BAK (BAK1); 3=BAX; 4=BCL2; 5=BCLXL (BLC2L1); 6=BID; 
% 7=PC3 (CASP3); 8=PC9 (CASP9); 9=MCL1; 
% 10=SMAC; 11=XIAP; 12=BIM (BCL2L11); 13=PUMA (BBC3); 14=NOXA (PMAIP1)
% ***NB. Proteins must be in correct order as specified above***

% classify cell lines according to sensitivity: 1=sensitive; 2=resistant
response_classification = [1,2,2,1,1,1];

% Convert input table to array
mb_proteins = table2array(mb_proteins_table);
names_proteins = mb_proteins_table.Properties.VariableNames;
names_cellLines = mb_proteins_table.Properties.RowNames;

% Calculate functional groups
% ***NB. Proteins must be in correct order as specified above***
[mb_funcGroups, names_funcGroups] = calc_funcGroups(mb_proteins);

% Can continue with functional groups or individual protein levels
mb_pca = mb_funcGroups;
names_pca = names_funcGroups;
%mb_pca = mb_proteins;
%names_pca = names_proteins;

%% 2. Normalisation
% Calculate mean, std dev so you can apply these to
% test datasets later on
stdvalue_cellExp = std(mb_pca);
meanvalue_cellExp = mean(mb_pca);

% Zscore normalise data around columns (normalise protein values; data - mean/stddev)
mb_pca_norm = zscore(mb_pca);

%% 3. PCA
% "coeff" are the principal component vectors (also known as loadings or 
% coefficients). These are the eigenvectors of the covariance matrix.
% They identify the contribution of each feature (protein) to each PC
% See https://www.mathworks.com/matlabcentral/answers/270329-how-to-select-the-components-that-show-the-most-variance-in-pca
[coeff_pca,scores,pcavars,tsquared,explained,mu] = pca(mb_pca_norm);

% Verify that original input data * coeff_pca = scores
mb_pca_norm * coeff_pca
scores

%% 3a. PCA plot (Figure 2B)

% Assign scores of 1st 3 PCs to the axes of the plot
x_pca = scores(:,1);
y_pca = scores(:,2);
z_pca = scores(:,3);

% Name figure so can redraw on it later
pca_plot = figure; 

% 3d scatter plot (1st 3 PCs)
for i=1:length(response_classification)
        if response_classification(i)==1
            scatter3(x_pca(i),y_pca(i),z_pca(i),'k','filled'); hold on;
        elseif response_classification(i)==2
            scatter3(x_pca(i),y_pca(i),z_pca(i),'r','filled'); hold on;
        elseif response_classification(i)==3
            scatter3(x_pca(i),y_pca(i),z_pca(i),'b','filled'); hold on;
        end
end

%%% Label axes ('explained' variable is the % that each PC contributes 
%%% to total variance of original data)
xlabel(['PC1 (' num2str(round(explained(1),1)) '%)'])
ylabel(['PC2 (' num2str(round(explained(2),1)) '%)'])
zlabel(['PC3 (' num2str(round(explained(3),1)) '%)'])

%%% Label cell lines
dx = 0.1; % displacement so text does not overlay data points
text(x_pca+dx,y_pca+dx,z_pca+dx,names_cellLines)

%%% Set view
view()
%view(90,0) % view PC2 & PC3

%% 3b. Kaiser criterion and feature contributions to each PC
% Calculate how many PCs need to be included in further analyses according
% to Kaiser Criterion (those >1). Reference in Vetna et al 2020
KC = sum(pcavars>1);

% Pareto chart of contribution of each PC to total variance
% Only the first 95% of the cumulative distribution is displayed.
figure; pareto(explained)

%%% Graph weight coefficients (loadings)
% coeff identifies the contribution of each feature (protein) to each PC
% Each feature influences each PC in a different way, so you can only say 
% "feature x, y, z have the highest influence on PC1", "feature y, t have the 
% highest influence on PC2"…

% biplot (Figure 2E)
% From Matlab: "A biplot allows you to visualize the magnitude and sign of 
% each variable's contribution to the first two or three principal 
% components, and how each observation is represented in terms of those 
% components."
figure
biplot(coeff_pca(:,1:3), 'Scores', scores(:,1:3), 'Varlabels', names_pca,...
    'Marker','.', 'MarkerSize', 20)
xlabel(['PC1 (' num2str(round(explained(1),1)) '%)'])
ylabel(['PC2 (' num2str(round(explained(2),1)) '%)'])
zlabel(['PC3 (' num2str(round(explained(3),1)) '%)'])
view(90,0) % view PC2 & PC3

% barchart of PC contributions
figure
% Set axes limits
limgrapghpos = round(max(max(coeff_pca)),1) +0.1;
limgraphneg = round(min(min(coeff_pca)),1) -0.1;
for pcn = 1:KC
    subplot(1,KC,pcn)
    barh(coeff_pca(:,pcn))
    xlim([limgraphneg limgrapghpos])
    title(['PC' num2str(pcn)])
    hold on;
    if pcn==1
        set(gca,'YTickLabel',names_pca);
    end
end
hold off

%% 4. Linear discriminant analysis (LDA)
% class = classify(Test_Set,Training_Set,train_label,'linear')
% "In order for the covariance matrix of TRAINING to be positive definite, 
% you must at the very least have more observations than variables in Test_Set."

% Use the first KC PCs to determine the LDA space (specify 'linear')
[class,error,posterior,logp,coeff_lda] = ...
    classify(scores(:,1:KC), scores(:,1:KC), response_classification, 'linear');
comparison1 = class;

% Count how many cell lines are correctly classified
countok=0;
for i=1:length(response_classification)
    if response_classification(i)==comparison1(i)
        countok=countok+1;
    end
end

countok

%% 4a. Plot LDA in 2d space (Figure 2C)
% Define a 2D grid of dots and apply LDA to those dots based on 
% response_classifications 
[X,Y] = meshgrid(linspace(min(scores(:,1)),max(scores(:,1)),100),...
    linspace(min(scores(:,2)),max(scores(:,2)),100));
X = X(:); Y=Y(:);
[class,error,posterior,logp,coeff_lda] = ...
    classify([X Y], scores(:,1:2), response_classification, 'linear');

lda_plot = figure('Name','PC1&2');
% first plot the 2d grid
lda1 = gscatter(X,Y,class,'krb','.',1,'off');
hold on
% Now plot our cell lines in this mesh grid
lda2 = gscatter(scores(:,1),scores(:,2), response_classification,'krb', ...
    '.',20,'off');
set(lda1,'LineWidth',2)
legend(lda2,'high sensitivity','low sensitivity','Location','NW')
xlabel(['PC1 (' num2str(round(explained(1),1)) '%)'])
ylabel(['PC2 (' num2str(round(explained(2),1)) '%)'])
% Label the cell lines on the figure
dx = 0.1; % displacement so the text does not overlay the data points
text(scores(:,1)+dx,scores(:,2)+dx,names_cellLines)

% Do the same on the 2nd and 3rd PC axes, as this is where we see
% the separation  
[X,Y] = meshgrid(linspace(min(scores(:,2))-0.25,max(scores(:,2))+0.25,100),...
    linspace(min(scores(:,3))-0.25,max(scores(:,3))+0.25,100));
X = X(:); Y=Y(:);
[class,error,posterior,logp,coeff_lda] = ...
    classify([X Y], scores(:,2:3), response_classification, 'linear');
lda_plot2 = figure('Name','PC2&3');
lda1 = gscatter(X,Y,class,'krb','.',0.1,'off');
hold on
lda2 = gscatter(scores(:,2),scores(:,3), response_classification,'krb', ...
    '.',20,'off');
set(lda1,'LineWidth',2)
legend(lda2,'sensitive','resistant','Location','NW')
xlabel(['PC2 (' num2str(round(explained(2),1)) '%)'])
ylabel(['PC3 (' num2str(round(explained(3),1)) '%)'])
dx = 0.1; 
text(scores(:,2)+dx,scores(:,3)+dx,names_cellLines)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5: Import patient sample data and plot on existing PCA and LDA space (Figure 2D)
% Read in new samples
% NB. Proteins must be in same order as listed above!!!
northcott_protein_table = readtable('Pt data_no CHLA.xlsx','ReadRowNames',true);

% Convert input table to array
northcott_proteins = table2array(northcott_protein_table);
names_proteins = northcott_protein_table.Properties.VariableNames;
names_northcott = northcott_protein_table.Properties.RowNames;

% Calculate functional groups
[testCL_fG, names_funcGroups] = calc_funcGroups(northcott_proteins);

% Z-score normalise using mean, stdev from original dataset
testCL_fGnorm = (testCL_fG-meanvalue_cellExp)./stdvalue_cellExp;

% Multiply by coeff_pca to obtain "scores" = dimensions for PCA space
testCL_scores = testCL_fGnorm/coeff_pca';

% Perform LDA on testCL to classify it ('test_class')
[test_class,error,posterior,logp,coeff_lda] = ...
        classify(testCL_scores(:,1:KC), scores(:,1:KC), response_classification, 'linear');

% For each sample in imported data
for i = 1:size(northcott_proteins,1)
    
    if test_class(i) ==1, dotColour = 'k';
    elseif test_class(i) ==2, dotColour = 'r';
    end
    
    % Plot and label new point on 2D LDA space (PC2vPC3)
    figure(lda_plot2)
    scatter(testCL_scores(i,2),testCL_scores(i,3), dotColour,'x');
    legend off
    
    % Plot and label this point in the PC space
    figure(pca_plot)
    hold on
    scatter3(testCL_scores(i,1),testCL_scores(i,2),testCL_scores(i,3),45,dotColour,'x')
    
end
end