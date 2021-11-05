%% Two-member species simulations with SteadyCom
% The scripts used in setup has been compiled into specific functions to
% reduce clutter and emphasize the analysis part of the process instead.

% clear workspace
clear;

initCobraToolbox(0);
changeCobraSolver('ibm_cplex', 'LP');

%% load model
modelEco = readCbModel(['models' filesep 'iML1515.xml']);
ecModelEco = readCbModel(['models' filesep 'eciML1515_batch.xml']);

% Using two variants of each exchange rxn is a bit tedious and very
% cumbersome especially in the context of communities of multiple GEMs so
% we will be collapsing the forward and reverse exchange rxns of each
% metabolite back into one.
ecModelEco = mergeExcRxns(ecModelEco);

% We will be using the exchange boundaries defined in the tutorial.
modelEco = setBounds(modelEco);
ecModelEco = setBounds(ecModelEco);

% You can verify the ecModel has its bounds modified correctly using:
% sol = solveLP(ecModelEco); printFluxes(ecModelEco, sol.x) %RAVEN function

%% Defining the gene-associated reactions for knockout
% argA, cysE, glyA, hisB, ilvA, leuB, lysA, metA, pheA, proA, serA, thrC, trpC, tyrA
% R, C, G, H, I, L, K, M, F, P, S, T, W, Y
% Arg, Cys, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr
% b2818, b3607, b2551, b2022, b3772, b0073, b2838, b4013, b2599, b0243, b2913, b0004, b1262, b2600
% {'ACGS'}, {'SERAT'}, {'GHMT2r', 'THFAT'}, {'IGPDH', 'HISTP'}, {'DAPDC'}, {'HSST'}, {'PPNDH'}, {'G5SD'}, {'AHGDx', 'PGCD', 'ARHGDx'}, {'THRS'}, {'IGPS', 'PRAIi'}, {'PPND'}


% We will be referring to the 2014 study done previously by Mee et al. that
% explored the syntrophies between obligate auxotrophic pairs of E. coli
% strains. In this setup, there are a total of 14 different amino acid
% auxotrophs.

genes = {'b2818', 'b3607', 'b2551', 'b2022', 'b3772', 'b0073', 'b2838', 'b4013', 'b2599', 'b0243', 'b2913', 'b0004', 'b1262', 'b2600'};
geneNames =  {'argA', 'cysE', 'glyA', 'hisB', 'ilvA', 'leuB', 'lysA', 'metA', 'pheA', 'proA', 'serA', 'thrC', 'trpC', 'tyrA'};
mets = {'arg__L', 'cys__L', 'gly', 'his__L', 'ile__L', 'leu__L', 'lys__L', 'met__L', 'phe__L', 'pro__L', 'ser__L', 'thr__L', 'trp__L', 'tyr__L'};
rxns = {{'ACGS'}, {'SERAT'}, {'ALATA_D2'; 'ALATA_L2'; 'GHMT2r'; 'THRA2'; 'THRA'; 'THFAT'}, {'IGPDH'; 'HISTP'}, {'THRD_L'}, {'IPMD'}, {'DAPDC'}, {'HSST'}, {'PPNDH'}, {'G5SD'}, {'AHGDx'; 'PGCD'; 'ARHGDx'}, {'THRS'}, {'IGPS'; 'PRAIi'}, {'PPND'}};

% First, let us establish these auxotrophic strains and test their
% predicted phenotype in silico.

% For the iML1515 GEM this is done by setting the bounds of rxns without
% isozymes to 0. However, for some gene KOs (glyA and serA) growth persists.
disp('Using iML1515 by fixing rxns with no isozymes to zero flux, we see:');

% Some rxns were actually not sufficient to perturb growth if we had
% knocked out only the rxns without isozymes (for ilvA and leuB) so we
% knocked out the associated rxns since they are the only rxn associated to
% each gene.
for i = 1:length(genes)
    model = changeRxnBounds(modelEco, rxns{i}, 0, 'b');
    sol = optimizeCbModel(model); 
    if sol.f > 0
        disp(['The ' geneNames{i} '-deletion strain can still grow.']);
    end
end

% For the ecModel this is done by setting the bounds of the enzyme usage 
% of that enzyme to 0. Likewise, for some gene KOs growth persists. As is
% the case we saw previously, we have to additionally modify the rxn bounds
% for ilvA and leuB-associated rxns to perturb growth.
disp('Using eciML1515 by fixing enzyme usage of the KOed gene to zero, we see:');

for i = 1:length(genes)
    protein = modelEco.proteinisuniprotID(ismember(modelEco.genes, genes(i)));
    model = changeRxnBounds(ecModelEco, strcat('draw_prot_', protein), 0, 'b');
    sol = optimizeCbModel(model); 
    if sol.f > 0
        model = changeRxnBounds(model, model.rxns(contains(model.rxns, rxns{i})), 0, 'b');
        sol = optimizeCbModel(model); 
        if sol.f > 0
            disp(['The ' geneNames{i} '-deletion strain can still grow.']);
        end
    end
end

%% Generate the auxotrophic strains and define the community for every pair
% In this section we will create a community model structure for each amino
% acid auxotroph pair. Currently, the script `getDeletionStrains` will
% create only unique pairs (for a 14 x 14 matrix include indices above the
% diagonal). As this will take a long time to run, it is not recommended
% but if you wish to proceed you may uncomment the following lines.
% Otherwise, you can continue to the next section where we will be
% inspecting a subset of this.

% [auxModels, commModel]      = getDeletionStrains(modelEco, geneNames, mets, rxns, 'rxns');
% [auxEcModels, commEcModel]  = getDeletionStrains(ecModelEco, geneNames, mets, genes, 'genes');

% If you'd like to inspect all the strains but for a specific pair such as
% with the methionine auxotrophic strain as species 1 in each pair, you can
% download the following files below to run the analysis without having to
% spend time generating the community models:

%websave('models/methionine_commModels.zip','https://chalmersuniversity.box.com/shared/static/oxg2hdgtczaxmozfxzhpdi1mqfhsivc8.zip');
%cd models
%unzip('methionine_commModels.zip');
%cd ..
%load('models/methionine_commModels');

% Otherwise you can also do pair simulations using the strains that were
% used by Chan et al. in the SteadyCom tutorial as defined again below:

% Reactions to be knocked out for amino acid auxotrophy:
metA = {'HSST'};  % essential for methionine biosynthesis
argH = {'ARGSL'};  % essential for arginine biosynthesis
lysA = {'DAPDC'};  % essential for lysine biosynthesis
ilvE = {'PPNDH'};  % essential for phenylalanine biosynthesis

rxns = [metA, argH, lysA, ilvE];
geneNames = {'metA', 'argH', 'lysA', 'ilvE1'};
mets = {'met__L', 'arg__L', 'lys__L', 'phe__L'};
genes = {'b4013', 'b3960', 'b2838', 'b3770'};

% generate community models for SteadyCom tutorial strains
[auxModels, commModel]      = getDeletionStrains(modelEco, geneNames, mets, rxns, 'rxns');
[auxEcModels, commEcModel]  = getDeletionStrains(ecModelEco, geneNames, mets, genes, 'genes');

%% Effect of abundance on community growth rate in SteadyCom
% Here we can assess the influence of relative compositions on the growth
% rate by fixing the abundance of the first member (if you are using the
% SteadyCom tutorial strains, that is the methionine auxotroph) and by
% association the abundance of the second strain. This is simple in the
% case of two members but can become awry with three or more.

% define the abundance values
pcX = 0.0:0.05:1;

% initialize the growth rate values we would like to collect
GR = zeros(length(commModel),length(pcX));
GR_Ec = zeros(length(commModel),length(pcX));

options             = struct();
options.GRguess     = 0.5;  % initial guess for max. HR in the bisection
options.GRtol       = 1e-6;  % tolerance for final GR
options.algorithm   = 1;  % use the default algorithm (simple guessing for bounds, then matlab fzero)

for i = 1:length(commModel) % we are just running for one row
	for j = 1:length(pcX)
	    options.BMcon = [1 0; 0 1];
	    options.BMrhs = [pcX(j); 1-pcX(j)]; %fixing the relative abundances
	    options.BMcsense = ['E', 'E']; %set as equality
	    [sol, result] =     SteadyCom(commModel{1,i},options);
	    GR(i,j) = result.GRmax;
	end
end

% For the ecModel's case the growth rate is lower due to the protein
% constraint, so we will set the initial GR guess value in the optimization
% at 0.25 closer to the predicted growth rate with the ecModel
options.GRguess     = 0.25;

for i = 1:length(commEcModel)
	for j = 1:length(pcX)
	    options.BMcon = [1 0; 0 1];
	    options.BMrhs = [pcX(j); 1-pcX(j)];
	    options.BMcsense = ['E', 'E'];
	    [sol, result] =     SteadyCom(commEcModel{1,i},options);
	    GR_Ec(i,j) = result.GRmax;
	end
end