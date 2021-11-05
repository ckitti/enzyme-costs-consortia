clear;

initCobraToolbox(0)
changeCobraSolver('ibm_cplex', 'LP');

%% Model Construction


% load the ec-iML1515 in the COBRA toolbox
modelEco = readCbModel(['models' filesep 'eciML1515_batch.xml']);

% Before proceeding, make sure all bounds are same between this model and
% the iML1515 from which it is derived:

%modelEco.ub(6084) = 0.2;

modelEco = changeRxnBounds(modelEco, modelEco.rxns(modelEco.ub==Inf), 1000, 'u');
% close all exchange reactions
modelEco = changeRxnBounds(modelEco, modelEco.rxns(contains(modelEco.rxns, 'EX_')), 0, 'u');

% Now we will get the bounds from iML1515 as used in the tutorial
fid         = fopen(['models/exchangeBoundaries_tutorial.txt']);
loadedData  = textscan(fid,'%s %f %f','delimiter','\t', 'HeaderLines',1); fclose(fid);
exRxns      = loadedData{1};
exRxnsLb    = loadedData{2};
exRxnsUb    = loadedData{3};

% In the ecModel, all originally reversible reactions are divided into two
% separate unidirectional reactions: the Forward and Reverse reactions. The
% latter is denoted with '_REV' in the reaction ID.

% Fix forward rxns with upper bounds of their respective original reactions
idx = find(exRxnsUb);
modelEco = changeRxnBounds(modelEco,exRxns(idx),exRxnsUb(idx), 'u');

% Fix reverse rxns with lower bounds of their respective original reactions
idx = find(exRxnsLb);
exRxnsRev = strcat(exRxns, '__');
exRxnsRev = strrep(exRxnsRev, '_e__', '_e_REV');
modelEco = changeRxnBounds(modelEco, exRxnsRev(idx), -exRxnsLb(idx), 'u');

% set the uptake bounds to closely replicate the original study as possible

% These boundaries are based on the original study
modelEco = changeRxnBounds(modelEco, 'EX_glc__D_e_REV', 6, 'u');
modelEco = changeRxnBounds(modelEco, 'EX_o2_e_REV', 18.5, 'u');
modelEco = changeRxnBounds(modelEco, 'EX_cbl1_e_REV', 0.01, 'u');

% These bounds were adapted to the updated E. coli model used in the tutorial
modelEco = changeRxnBounds(modelEco, 'EX_slnt_e_REV', 0, 'u');
modelEco = changeRxnBounds(modelEco, 'EX_sel_e_REV', 0, 'u');
modelEco = changeRxnBounds(modelEco, 'EX_ni2_e_REV', 1000, 'u');

%% Check model structure for consistency

% remove the redundant compartment format from e.g., '_c' to '[c]'
modelEco.mets = regexprep(modelEco.mets, '_e[', '[');
modelEco.mets = regexprep(modelEco.mets, '_c[', '[');
modelEco.mets = regexprep(modelEco.mets, '_p[', '[');
% convert empty cells in cell arrays, if any, to empty strings

fieldToBeCellStr = {'metFormulas'; 'genes'; 'rules'; 'metNames'; 'rxnNames'; 'subSystems'};

for j = 1:numel(fieldToBeCellStr)

modelEco.(fieldToBeCellStr{j})(cellfun(@isempty, modelEco.(fieldToBeCellStr{j}))) = {''};

end

%% Defining the gene-associated reactions for knockout

% modelEco = addReaction(modelEco,{'METt3pp',''},'met__L[c] + h[c] => met__L[p] + h[p]');

% Reactions to be knocked out for amino acid auxotrophy:

argH = {'ARGSLNo1'};  % essential for arginine biosynthesis
lysA = {'DAPDCNo1'};  % essential for lysine biosynthesis
metA = {'HSSTNo1'};  % essential for methionine biosynthesis
ilvE1 = {'PPNDHNo1'};  % essential for phenylalanine biosynthesis (gene b3770)

% find(contains(modelEco.genes, 'b3770'));
% 18
rxnIdx = find(contains(modelEco.rules, 'x(18)'));
ilvE2 = modelEco.rxns(rxnIdx);  % essential for phenylalanine biosynthesis (gene b3770)

% find(contains(modelEco.genes, 'b2599'));
% 87
rxnIdx = find(contains(modelEco.rules, 'x(87)'));
pheA = modelEco.rxns(rxnIdx);  % essential for phenylalanine biosynthesis (gene b2599)

%%
% methionine auxotroph
%Ec1 = changeRxnBounds(modelEco, [argH], 0, 'b');
Ec1 = changeRxnBounds(modelEco, [metA; lysA], 0, 'b');

% phenylalaline auxotroph
Ec2 = changeRxnBounds(modelEco, [metA; ilvE1], 0, 'b');

% lysine auxotroph
Ec3 = changeRxnBounds(modelEco, [lysA; ilvE1], 0, 'b');

%Ec2 = changeRxnBounds(modelEco, ilvE1, 0, 'b');
%Ec2 = changeRxnBounds(modelEco, ilvE2, 0, 'b');
%Ec2 = changeRxnBounds(modelEco, pheA, 0, 'b');


%% Defining the community model boundaries

% Now none of the four organisms can grow alone and they must cross feed 
% each other to survive. We will now have to identify the extracellular 
% metabolites, their corresponding exchange reactions and uptake rates for 
% the E. coli model, which are used later to constrain the community model:

% identify the extracellular metabolites (met[e])
metEx = strcmp(getCompartment(modelEco.mets),'e');

% find the corresponding exchange reactions
rxnExAll = find(sum(modelEco.S ~= 0, 1) == 1);
[rxnEx, ~] = find(modelEco.S(metEx, rxnExAll)');  % need to be in the same order as metEx
rxnEx = rxnExAll(rxnEx);

% find the extracellular metabolites, which were added in the iML1515 update
% that are not accounted for
modelEco = buildRxnEquations(modelEco);
% modelEco.rxnEquations(rxnEx);
extMetsExc = strrep(modelEco.rxnEquations(rxnEx), ' -->', '');
extMetsExc = strrep(extMetsExc, ' <==>', '');

% filter out aforementioned metabolites and add exchange reactions for them
extMets = modelEco.mets(metEx);
extMets = setdiff(extMets,extMetsExc);
modelEco = addExchangeRxn(modelEco, extMets, zeros(1,10), zeros(1,10));

% obtain the exchange boundaries used to constrain the community model
idx = findRxnIDs(modelEco, strcat('EX_',extMets));

% exchange rate
% this method to filter out the reverse reactions only works given how the
% model structure is arranged. In the array of exchange reaction IDs, the
% reverse reactions are evenly spaced (every other term)
rxnExREV = rxnEx(2:2:end)'; 
% get indices for all reverse reactions and the newly added exchange reactions
rxnEx = [rxnExREV; idx];
% match the respective uptake rates with the exchange reactions, again all
% bounds that were negative in the previous model is a positive upper bound
% in the ecModel
lbEx = -modelEco.ub(rxnEx);


%% Creating the community model structure

% Create a community model with the four E. coli tagged as 'Ec1', 'Ec2', 
% 'Ec3', 'Ec4' respectively by calling |createMultipleSpeciesModel|.

nameTagsModel = {'Ec1'; 'Ec2'; 'Ec3'};

EcCom = createMultipleSpeciesModel({Ec1; Ec2; Ec3}, nameTagsModel);
EcCom.csense = char('E' * ones(1,numel(EcCom.mets)))';
 
% Retrieve the names and ids for organism/community exchange reactions/metabolites 
% which are necessary for computation:

[EcCom.infoCom, EcCom.indCom] = getMultiSpeciesModelId(EcCom, nameTagsModel);
disp(EcCom.infoCom);

%% Identify the biomass reaction

% Incorporate also the names and indices for the biomass reactions which 
% are necessary for computing growth:

% specify biomass reaction names from each individual model
rxnBiomass = strcat(nameTagsModel, 'BIOMASS_Ec_iML1515_core_75p37M');  
rxnBiomassId = findRxnIDs(EcCom, rxnBiomass);  % corresponding ids in EcCom

% .spBm for organism biomass reactions
EcCom.infoCom.spBm  = rxnBiomass;  %identify rxn name
EcCom.indCom.spBm   = rxnBiomassId; %identify rxn index

%% Setting the boundaries of the community model structure 

% Before we find the max growth rate, we will first set community and 
% organism-specific uptake rates to be the same as in the original model:

% map the metabolite name
[yn, id] = ismember(strrep(modelEco.mets(metEx), '[e]', '[u]'), EcCom.infoCom.Mcom);  % map the metabolite name
assert(all(yn));  % must be a 1-to-1 mapping

% rename the reverse exchange reactions as normal reactions for the
% community model structure
rxnExID = strrep(modelEco.rxns(rxnEx), '_e_REV', '[e]');
rxnExID = strrep(rxnExID, '[e]', '[u]');

% map the exchange reaction name
[yn, id] = ismember(EcCom.infoCom.EXcom, rxnExID);  
assert(all(yn));  % must be a 1-to-1 mapping

% assign community uptake bounds
EcCom.lb(EcCom.indCom.EXcom(:,1))   = lbEx(id);  
EcCom.ub(EcCom.indCom.EXcom(:,1))   = 1e5;

% assign organism-specific uptake bounds
EcCom.lb(EcCom.indCom.EXsp)         = repmat(lbEx(id), 1, 3);

%% Setting the boundaries of the individual members 

% allow production of anything for each member
exRate = 1;
rxns = {'IEX_met__L[u]tr', 'IEX_lys__L[u]tr', 'IEX_arg__L[u]tr', 'IEX_phe__L[u]tr'};
rxnIdx = find(contains(EcCom.rxns, rxns));

%EcCom = changeRxnBounds(EcCom, EcCom.rxns(rxnIdx), -exRate, 'l');

% Ec1
EcCom = changeRxnBounds(EcCom, {'Ec1IEX_met__L[u]tr'}, -exRate, 'l');
EcCom = changeRxnBounds(EcCom, {'Ec1IEX_lys__L[u]tr'}, -exRate, 'l');
%EcCom = changeRxnBounds(EcCom, {'Ec1IEX_arg__L[u]tr'}, -exRate, 'l');

% Ec2
EcCom = changeRxnBounds(EcCom, {'Ec2IEX_phe__L[u]tr'}, -exRate, 'l');
%EcCom = changeRxnBounds(EcCom, {'Ec2IEX_lys__L[u]tr'}, -exRate, 'l');
%EcCom = changeRxnBounds(EcCom, {'Ec2IEX_arg__L[u]tr'}, -exRate, 'l');
EcCom = changeRxnBounds(EcCom, {'Ec2IEX_met__L[u]tr'}, -exRate, 'l');
%EcCom = changeRxnBounds(EcCom, {'Ec2IEX_lys__L[u]tr'}, -exRate, 'l');

% Ec3
%EcCom = changeRxnBounds(EcCom, {'Ec3IEX_arg__L[u]tr'}, -exRate, 'l');
EcCom = changeRxnBounds(EcCom, {'Ec3IEX_lys__L[u]tr'}, -exRate, 'l');
EcCom = changeRxnBounds(EcCom, {'Ec3IEX_phe__L[u]tr'}, -exRate, 'l');

% allow production of anything for each member
EcCom.ub(EcCom.indCom.EXsp(:)) = 1000;

% print the community uptake bounds to check
printUptakeBoundCom(EcCom, 1);

%% set abundance of first pair member for SteadyCom

pcX = 0.0:0.05:1
GR = [];

options             = struct();
options.GRguess     = 0.25;  % initial guess for max. HR in the bisection
options.GRtol       = 1e-6;  % tolerance for final GR
options.algorithm   = 1;  
% use the default algorithm (simple guessing for bounds, then matlab fzero)

for i = 1:length(pcX)
    options.BMcon = [1 0 0 ;0 1 1];
    options.BMrhs = [pcX(i); 1-pcX(i)];
    options.BMcsense = ['E', 'E'];
    [sol, result] =     SteadyCom(EcCom,options);
    GR = [GR result.GRmax];
end

%% Finding the max growth rate using SteadyCom


options             = struct();
options.GRguess     = 0.5;  % initial guess for max. HR in the bisection
options.GRtol       = 1e-6;  % tolerance for final GR
options.algorithm   = 1;  
% use the default algorithm (simple guessing for bounds, then matlab fzero)

[sol, result] =     SteadyCom(EcCom,options);

for jSp = 1:3
    fprintf('X_%s:  %.6f\n', EcCom.infoCom.spAbbr{jSp}, result.BM(jSp));
end
disp(result);

iter = [0, result.iter0, NaN; result.iter];
for j = 0 : size(iter, 1)
    if j == 0
        fprintf('#iter\tgrowth rate (mu)\tmax. biomass (sum(X))\tmu * sum(X)\tmax. infeasibility\tguess method\n');
    else
        fprintf('%5d\t%16.6f\t%21.6f\t%11.6f\t%18.6e\t%d\n', iter(j,:))
    end
end

%% Analyze flux variability with SteadyComFVA

% To analyze the variability of the organism abundance at various growth 
% rates, we will use the function |SteadyComFVA|.

% percentage of maximum total biomass of the community required. 
% 100 for sum(biomass) = 1 (1 is the default total biomass)
options.optBMpercent = 100;  
n = size(EcCom.S, 2);  % number of reactions in the model

% options.rxnNameList is the list of reactions subject to FVA. Can be 
% reaction names or indices. Use n + j for the biomass variable of the 
% j-th organism. Alternatively, use {'X_j'} for biomass variable of the 
% j-th organism or {'X_Ec1'} for Ec1 (the abbreviation in EcCom.infoCom.spAbbr)

options.rxnNameList = {'X_Ec1'; 'X_Ec2'; 'X_Ec3'};

% define the growth rates we would like to perform FVA for
options.optGRpercent = [89:0.2:99, 99.1:0.1:100];

% perform FVA at various fractions of the maximum growth rate
[fvaComMin,fvaComMax] = SteadyComFVA(EcCom, options);

%% If you want to just plot SteadyComFVA results

% vector of growth rates tested
grComV = options.optGRpercent / 100; % * result.GRmax;
lgLabel = {'{\itEc1 }';'{\itEc2 }';'{\itEc3 }'};
col = [ 95 135 255; 255 0 0; 0 235 0 ]/255;  % color
f = figure;
x = [grComV(:); flipud(grComV(:))];
hold on
for j = 1:3
    y = [fvaComMin(j, :), fliplr(fvaComMax(j, :))];
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
    p(j, 1).Color = col(j, :);
end
tl(1) = title('\underline{SteadyCom}', 'Interpreter', 'latex');
tl(1).Position = [0.7 2 0];
ax(1) = gca;
ax(1).XTick = grComV(1):0.01:grComV(end);
ax(1).YTick = -1:0.1:1;
xlim([grComV(1) grComV(end)])
ylim([0 1])

lg = legend(lgLabel);
lg.Box = 'off';
yl(1) = ylabel('Relative abundance');
xl(1) = xlabel('Community growth rate (h^{-1})');

%% Analyze Pairwise Relationship Using SteadyComPOA

if exist('POAtmp')
    rmdir('POAtmp','s');
end

options.optGRpercent = [99.9 99 90 75 50 0];  % analyze at these percentages of max. growth rate

a = 0.001*(1000.^((0:14)/14));
options.Nstep = sort([a (1-a)]);

[POAtable, fluxRange, Stat, GRvector] = SteadyComPOA(EcCom, options);

% Plot the results (see also Figure 3 in ref. [1]):

nSp = 3;
spLab = {'{\it Ec1 }';'{\it Ec2 }'; '{\it Ec3}'};
mark = {'A', 'B', 'D', 'C', 'E', 'F'};
nPlot = 0;
for j = 1:nSp
    for k = 1:nSp
        if k > j
            nPlot = nPlot + 1;
            ax(j, k) = subplot(nSp-1, nSp-1, (k - 2) * (nSp - 1) + j);
            hold on
            for p = 1:size(POAtable{1, 1}, 3)
                x = [POAtable{j, j}(:, :, p);POAtable{j, j}(end:-1:1, :, p);...
                    POAtable{j, j}(1, 1, p)];
                y = [POAtable{j, k}(:, 1, p);POAtable{j, k}(end:-1:1, 2, p);...
                        POAtable{j, k}(1, 1, p)];
                plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2)
            end
            xlim([0.001 1])
            ylim([0.001 1])
            ax(j, k).XScale = 'log';
            ax(j, k).YScale = 'log';
            ax(j, k).XTick = [0.001 0.01 0.1 1];
            ax(j, k).YTick = [0.001 0.01 0.1 1];
            ax(j, k).YAxis.MinorTickValues=[];
            ax(j, k).XAxis.MinorTickValues=[];
            ax(j, k).TickLength = [0.03 0.01];
            xlabel(spLab{j});
            ylabel(spLab{k});
            tx(j, k) = text(10^(-5), 10^(0.1), mark{nPlot}, 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end
lg = legend(strcat(strtrim(cellstr(num2str(options.optGRpercent(:)))), '%'));
lg.Position = [0.7246 0.6380 0.1700 0.2015];
lg.Box='off';
subplot(3, 3, 3, 'visible', 'off');
t = text(0.2, 0.8, {'% maximum';'growth rate'});
for j = 1:nSp
    for k = 1:nSp
        if k>j
            ax(j, k).Position = [0.15 + (j - 1) * 0.3, 0.8 - (k - 2) * 0.3, 0.16, 0.17];
            ax(j, k).Color = 'none';
        end
    end
end

