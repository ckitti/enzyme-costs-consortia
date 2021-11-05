function modelBounded = setBounds(model)
% script to set the bounds to match between the two model types


% if is an ecModel
if contains(model.modelName, 'ecModel')
    
    % Before proceeding, make sure all bounds are same between this model and
    % the iML1515 from which it is derived:

    model = changeRxnBounds(model, model.rxns(model.ub==Inf), 1000, 'u');
    model = changeRxnBounds(model, model.rxns(model.lb==-Inf), -1000, 'l');
    
    % close all exchange reactions
    model = changeRxnBounds(model, model.rxns(contains(model.rxns, 'EX_')), 0, 'b');
    
	% Now we will get the bounds from iML1515 as used in the tutorial
	fid         = fopen(['models/exchangeBoundaries_tutorial.txt']);
	loadedData  = textscan(fid,'%s %f %f','delimiter','\t','HeaderLines',1); fclose(fid);
	exRxns      = loadedData{1};
	exRxnsLb    = loadedData{2};
	exRxnsUb    = loadedData{3};
	
	% In the ecModel, all originally reversible reactions are divided into 	two
	% separate unidirectional reactions: the Forward and Reverse 	reactions. The
	% latter is denoted with '_REV' in the reaction ID.
	
	% Fix forward rxns with upper bounds of their respective original 	reactions
	idx = find(exRxnsUb);
	model = changeRxnBounds(model,exRxns(idx),exRxnsUb(idx), 'u');
	
	% Fix reverse rxns with lower bounds of their respective original 	reactions
	idx = find(exRxnsLb);
	model = changeRxnBounds(model, exRxns(idx), exRxnsLb(idx),'l');
end

% These boundaries are based on the original study
model = changeRxnBounds(model, 'EX_glc__D_e', -8, 'l');
model = changeRxnBounds(model, 'EX_o2_e', -8, 'l');
model = changeRxnBounds(model, 'EX_cbl1_e', -0.01, 'l');

% These bounds were added to iML1515; set them to 0 as well
model = changeRxnBounds(model, 'EX_slnt_e', 0, 'l');
model = changeRxnBounds(model, 'EX_sel_e', 0, 'l');

% Nickel uptake is essential in iML1515, so we will leave the bounds as is
modelBounded = changeRxnBounds(model, 'EX_ni2_e', -1000, 'l');

% Check model structure for consistency

% remove the redundant compartment format from e.g., '_c' to '[c]'
model.mets = regexprep(model.mets, '_e[', '[');
model.mets = regexprep(model.mets, '_c[', '[');
model.mets = regexprep(model.mets, '_p[', '[');
% convert empty cells in cell arrays, if any, to empty strings

% convert empty cells in cell arrays, if any, to empty strings

fieldToBeCellStr = {'metFormulas'; 'genes'; 'rules'; 'metNames'; 'rxnNames'; 'subSystems'};

for j = 1:numel(fieldToBeCellStr)

model.(fieldToBeCellStr{j})(cellfun(@isempty, model.(fieldToBeCellStr{j}))) = {''};

end


% Additional step for iML1515 update

% The following steps are to address a number of new extracellular
% metabolites that have been added since the first release of iAF1260 used
% in the original tutorial. If there are extracellular metabolites ([e]) that
% do not have a corresponding exchange reaction, the following generation of 
% the community model with createMultipleSpeciesModel will not work. 

% identify the extracellular metabolites (met[e])
metEx = strcmp(getCompartment(model.mets),'e');

% find the corresponding exchange reactions
rxnExAll    = find(sum(model.S ~= 0, 1) == 1);
[rxnEx, ~]  = find(model.S(metEx, rxnExAll)');  % need to be in the same order as metEx
rxnEx       = rxnExAll(rxnEx);

% Below is a quick way to identify those unaccounted metabolites:
model = buildRxnEquations(model);

% Display list of all exchange reactions by equation:
% model.rxnEquations(rxnEx)

% List of all extracellular metabolites with an exchange reaction
extMetsExc = strrep(model.rxnEquations(rxnEx), ' -->', '');
extMetsExc = strrep(extMetsExc, ' <==>', '');

% List all extracellular metabolites and filter out unaccounted ones 
extMets = model.mets(metEx);
extMets = setdiff(extMets,extMetsExc);

% Uncomment the following to see which reactions the unaccounted 
% extracellular metabolites are involved in:
% extMetsRxns = find(contains(model.rxnEquations, extMets))
% model.rxnEquations(extMetsRxns)

% Now add the exchange reactions but set bounds to zero to not confound the
% following results as our intention is to compare with the original study.
% We may adjust these later if needed.
modelBounded = addExchangeRxn(model, extMets, zeros(1,length(extMets)),  zeros(1,length(extMets)));
end