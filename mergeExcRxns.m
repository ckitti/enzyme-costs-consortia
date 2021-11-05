function [newModel, excRxns] = mergeExcRxns(model)
% script to collapse the forward and reverse exchange rxns back into one

excRxns = getExchangeRxns(model);
excRxns(end) = '';
excRxnsREV = excRxns((contains(excRxns, '_REV')));
excRxnsFWD = setdiff(excRxns, excRxnsREV);

for i = 1:length(excRxnsFWD)
    idx_REV = find(ismember(model.rxns, strcat(excRxnsFWD(i), '_REV')));
    idx_FWD = find(ismember(model.rxns, excRxnsFWD(i)));
    
    % get the bounds
    model.lb(idx_FWD) = -1 * model.ub(idx_REV);
    

end

% collapse the irreversible rxns
newModel = removeRxns(model, excRxnsREV);

end