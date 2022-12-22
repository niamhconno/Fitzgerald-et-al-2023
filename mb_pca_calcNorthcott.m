%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in individual protein expressions
% Perform PCA and LDA on measured protein values [plots both]
% Test whether knockout or over-expression of proteins changes their
% classification (resistant <-> sensitive) [output to screen]
% Manually verify this prediction [plots test cell line on LDA&PCA plots]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mb_pca_calcNorthcott
%% 1. Import individual protein values
mb_proteins_table = readtable('MB_PCA_input_individual_6.xlsx','ReadRowNames',true);

% Convert input table to array
mb_proteins = table2array(mb_proteins_table);
names_proteins = mb_proteins_table.Properties.VariableNames;
names_cellLines = mb_proteins_table.Properties.RowNames;
% Proteins:
% 1=APAF1; 2=BAK (BAK1); 3=BAX; 4=BCL2; 5=BCLXL (BLC2L1); 6=BID; 
% 7=PC3 (CASP3); 8=PC9 (CASP9); 9=MCL1; 
% 10=SMAC; 11=XIAP; 12=BIM (BCL2L11); 13=PUMA (BBC3); 14=NOXA (PMAIP1)

% Calculate functional groups
% ***NB. Proteins must be in correct order as specified above***
[mb_funcGroups, names_funcGroups] = calc_funcGroups(mb_proteins);

% Can continue with individual protein levels, or functional groups
mb_pca = mb_funcGroups;

%% 2. Normalisation
% Calculate mean, std dev so you can apply these to
% test datasets later on
stdvalue_cellExp = std(mb_pca);
meanvalue_cellExp = mean(mb_pca);

% Zscore normalise data around columns (normalise protein values; data - mean/stddev)
mb_pca_norm = zscore(mb_pca);

%% 3. PCA
% "coeff" are the principal component vectors.
% These are the eigenvectors of the covariance matrix.
% They identify the contribution of each feature (protein) to each PC
% See https://www.mathworks.com/matlabcentral/answers/270329-how-to-select-the-components-that-show-the-most-variance-in-pca
[coeff_pca,scores,pcavars,tsquared,explained,mu] = pca(mb_pca_norm);

% Verify that original input data * coeff_pca = scores
mb_pca_norm * coeff_pca
scores

%% 3a. PCA plot

% Assign scores of 1st 3 PCs to the axes of the plot
x_pca = scores(:,1);
y_pca = scores(:,2);
z_pca = scores(:,3);

% Colour datapoints according to sensitivity: 1=sensitve; 2=resistant
% Originally taken from responseTMZ variable from Cristiano's code
response_classification = [1,2,2,1,1,1]; 

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
%%% to total variance
xlabel(['PC1 (' num2str(round(explained(1),1)) '%)'])
ylabel(['PC2 (' num2str(round(explained(2),1)) '%)'])
zlabel(['PC3 (' num2str(round(explained(3),1)) '%)'])

%%% Label cell lines
dx = 0.1; % displacement so text does not overlay data points
text(x_pca+dx,y_pca+dx,z_pca+dx,names_cellLines)

%% 3b. Kaiser criterion
% Calculate how many PCs need to be included in further analyses according
% to Kaiser Criterion (those >1). Reference in Vetna et al 2020
KC = sum(pcavars>1);

% pareto chart of contribution of each PC to total variance
% Only the first 95% of the cumulative distribution is displayed.
figure; pareto(explained)

%%% Graph weight coefficients (loadings)
% coeff identifies the contribution of each feature (protein) to each PC
% Each feature influences each PC in a different way, so you can only say 
% feature x, y, z have the highest influence on PC1, feature y, t have the 
% highest influence on PC2…
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
        set(gca,'YTickLabel',names_funcGroups);
    end
end
hold off

%% 4. Linear discriminant analysis (LDA)
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


%% 4a. Plot LDA in 2d space
% Define a 2D grid of dots and apply LDA to those dots based on 
% response_classifications 
[X,Y] = meshgrid(linspace(-3,4,50),linspace(-2.5,2.5,50));
X = X(:); Y=Y(:);
[class,error,posterior,logp,coeff_lda] = ...
    classify([X Y], scores(:,1:2), response_classification, 'linear');

lda_plot = figure;
% first plot the 2d grid
lda1 = gscatter(X,Y,class,'krb','.',1,'off');
hold on
% Now plot our cell lines in this mesh grid
lda2 = gscatter(scores(:,1),scores(:,2), response_classification,'krb', ...
    'v^',[],'off');
set(lda1,'LineWidth',2)
legend(lda2,'sensitive','resistant','Location','NW')
% Label the cell lines on the figure
dx = 0.1; % displacement so the text does not overlay the data points
text(scores(:,1)+dx,scores(:,2)+dx,names_cellLines)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% XX: Plot new samples on existing PCA and LDA space
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
    
    % Plot and label new point on 2D LDA space 
    figure(lda_plot)
    scatter(testCL_scores(i,1),testCL_scores(i,2), dotColour,'x');
    legend off
    %text(testCL_scores(:,1)+dx,testCL_scores(:,2)+dx, names_northcott_test(i))
    
    % Plot and label this point in the PC space
    figure(pca_plot)
    hold on
    scatter3(testCL_scores(i,1),testCL_scores(i,2),testCL_scores(i,3),45,dotColour,'x')
    %text(testCL_scores(:,1)+dx,testCL_scores(:,2)+dx,testCL_scores(:,3)+dx,names_northcott_test(i))
    
end
%%
% Create table with sample names and classifications
sample_class = [names_northcott, num2cell(test_class)];

% Plot protein/funcGroup expression levels in sensitive vs resistant groups
figure
boxplot(northcott_proteins(:,strcmp(names_proteins,'SMAC')),sample_class(:,2))
boxplot(testCL_fG(:,strcmp(names_funcGroups,'BAK_BAX')),sample_class(:,2))

% Plot cell line protein/funcGroup expressions
figure
boxplot(mb_proteins(:,strcmp(names_proteins,'SMAC')), response_classification)
boxplot(mb_funcGroups(:,strcmp(names_funcGroups,'BAK_BAX')), response_classification)

end
