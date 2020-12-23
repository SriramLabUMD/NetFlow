function [model] = defineEcoliMediaM9(model)
mineralNames = {
    % Minerals in M9
    'nh4' % Ammonium
    'ca2' % Calcium
    'mg2' % Magnesium
    'k'   % Potassium
    % Sodium (dead-end rxn)
    'cl'  % Chloride
    'pi'  % Phosphate
    'so4' % Sulfate
    
    % Minerals required for model growth
    'cu2' % Copper II
    'mn2' % Manganese
    'mobd' % Molybdenum
    'cobalt2' % Cobalt
    'fe2' % Iron II
    'fe3' % Iron III
    'zn2' % Zinc 
    }; 

mineralIdxes = false(length(model.rxns),1);
for i = 1:length(mineralNames)
    tmp = sprintf('EX_%s_e',mineralNames{i});
    mineralIdxes = (mineralIdxes | strcmp(model.rxns,tmp));
end

otherNames1 = {
    'h2o' % H2O
    'h'   % H+
    };

otherIdxes1 = false(length(model.rxns),1);
for i = 1:length(otherNames1)
    tmp = sprintf('EX_%s(e)',otherNames1{i});
    otherIdxes1 = (otherIdxes1 | strcmp(model.rxns,tmp));
end

otherNames2 = {
    'cbl1'
    'ni2'
    };
otherIdxes2 = false(length(model.rxns),1);
for i = 1:length(otherNames2)
    tmp = sprintf('EX_%s_e',otherNames2{i});
    otherIdxes2 = (otherIdxes2 | strcmp(model.rxns,tmp));
end

co2Idx = strcmp(model.rxns,'EX_co2_e'); % CO2
glcIdx = strcmp(model.rxns,'EX_glc_e'); % Glucose
o2Idx = strcmp(model.rxns,'EX_o2_e'); % Oxygen

% Set model bounds
exRxns = strncmp('EX_',model.rxns,3);
model.lb(exRxns) = 0;
model.lb(mineralIdxes) = -1000;
model.lb(otherIdxes1) = -100;
model.lb(otherIdxes2) = -0.01;
model.lb(co2Idx) = -1000;
model.lb(glcIdx) = -10;
model.lb(o2Idx) = -18.5;
end