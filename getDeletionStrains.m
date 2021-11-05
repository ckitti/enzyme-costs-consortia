function [models, communityModel] = getDeletionStrains(model, geneNames, mets, ref, refType)
% script to set up the auxotrophic KO strains and the community models

models = cell(numel(ref), 1);
communityModel = cell(numel(ref), numel(ref));

switch refType
    
    case 'rxns'
        for i = 1:numel(ref)
            models{i} = changeRxnBounds(model, ref{i}, 0, 'b');
            models{i}.modelID = strcat('del_', geneNames{i});
        end
        
    case 'genes'
        idxIso = find(contains(model.rules, ' | '));
        for i = 1:numel(ref)
            geneIdx = find(ismember(model.genes, ref(i)));
            rxnIdx = find(contains(model.rules, ['x(' num2str(geneIdx) ')']));
            models{i} = changeRxnBounds(model, model.rxns(rxnIdx(end)), 0, 'b');
            sol = optimizeCbModel(models{i});
            if sol.f > 0
                %rxnIdx = setdiff(rxnIdx, idxIso);
                models{i} = changeRxnBounds(model, model.rxns(rxnIdx), 0, 'b');
            end
            models{i}.modelID = strcat('del_', geneNames{i});
        end
end
        
% identify the extracellular metabolites (met[e])
metEx = strcmp(getCompartment(model.mets),'e');

% find the corresponding exchange reactions
rxnExAll    = find(sum(model.S ~= 0, 1) == 1);
[rxnEx, ~]  = find(model.S(metEx, rxnExAll)');  % need to be in the same order as metEx
rxnEx       = rxnExAll(rxnEx);

% Update the indices for the corresponding exchange reactions
rxnExAll = find(sum(model.S ~= 0, 1) == 1);
[rxnEx, ~] = find(model.S(metEx, rxnExAll)');  % need to be in the same order as metEx
rxnEx = rxnExAll(rxnEx);

% Obtain the exchange rate for each reaction
lbEx = model.lb(rxnEx);


% Creating the community model structure

for i = 1:length(geneNames)
    for j = i:length(geneNames)    
        nameTagsModel = {'Ec1'; 'Ec2'};

        communityModel{i,j} = createMultipleSpeciesModel({models{i}; models{j}}, nameTagsModel);
        communityModel{i,j}.csense = char('E' * ones(1,numel(communityModel{i,j}.mets)))';
 
        % Retrieve the names and ids for organism/community exchange reactions/metabolites 
        % which are necessary for computation:

        [communityModel{i,j}.infoCom, communityModel{i,j}.indCom] = getMultiSpeciesModelId(communityModel{i,j}, nameTagsModel);
        disp(communityModel{i,j}.infoCom);
        
        % specify biomass reaction names from each individual model
        rxnBiomass = strcat(nameTagsModel, 'BIOMASS_Ec_iML1515_core_75p37M');  
        rxnBiomassId = findRxnIDs(communityModel{i,j}, rxnBiomass);  % corresponding ids in communityModel{i,j}

        % .spBm for organism biomass reactions
        communityModel{i,j}.infoCom.spBm  = rxnBiomass;  %identify rxn name
        communityModel{i,j}.indCom.spBm   = rxnBiomassId; %identify rxn index
        
        % Setting the boundaries of the community model structure 
        
        % Before we find the max growth rate, we will first set community and 
        % organism-specific uptake rates to be the same as in the original model:

        % map the metabolite name
        [yn, id] = ismember(strrep(model.mets(metEx), '[e]', '[u]'), communityModel{i,j}.infoCom.Mcom);  
        assert(all(yn));  % make sure it is a 1-to-1 mapping

        % assign community uptake bounds
        communityModel{i,j}.lb(communityModel{i,j}.indCom.EXcom(:,1))   = lbEx(id);  
        communityModel{i,j}.ub(communityModel{i,j}.indCom.EXcom(:,1))   = 1e5;

        % assign organism-specific uptake bounds
        communityModel{i,j}.lb(communityModel{i,j}.indCom.EXsp)         = repmat(lbEx(id), 1, 2);

        % allow production of anything for each member
        exRate = 1;
        
        % allow consumption of lacking amino acid
        communityModel{i,j} = changeRxnBounds(communityModel{i,j}, ['Ec1IEX_',mets{i},'[u]tr'], -exRate, 'l');
        communityModel{i,j} = changeRxnBounds(communityModel{i,j}, ['Ec2IEX_',mets{j},'[u]tr'], -exRate, 'l');
        
        % allow production of anything for each member
        communityModel{i,j}.ub(communityModel{i,j}.indCom.EXsp(:)) = 1000;
    end
end

end