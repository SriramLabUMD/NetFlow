function [modelSub,vSub] = extractSubnetwork(model,vRef,mets,options)
%Extacts a subnetwork containing only the given metabolites from a full
%model
%Created by Sean Mack 3/30/20
%Modified by Sean Mack 8/5/20
%-Added simplified naming for lumped sink/source reactions
lumpIdenticalRxns = 1;
keepRxns = [];

if exist('options','var')
    if (isfield(options,'filename'))
        filename = options.filename;
    end
    if (isfield(options,'lumpIdenticalRxns'))
        lumpIdenticalRxns = options.lumpIdenticalRxns;
    end
    if isfield(options,'keepRxns')
        keepRxns = options.keepRxns;
    end
end

if islogical(mets)
    sub = mets;
elseif isnumeric(mets)
    sub = ismember(model.mets,model.mets(mets));
else
    sub = ismember(model.mets,mets);
end
modelSub = removeMetabolites(model,model.mets(~sub));
subRxns = any(model.S(sub,:),1);
vSub = vRef(subRxns,:);
neg = all(vSub<=0,2);
if any(neg)
    vSub(neg,:) = -1*vSub(neg,:);
    modelSub.S(:,neg) = -1*modelSub.S(:,neg);
    modelSub.lb(neg) = -1*modelSub.lb(neg);
    modelSub.ub(neg) = -1*modelSub.ub(neg);
%     idx = find(neg);
%     for i = 1:length(idx)
%         modelSub.rxns{idx(i)} = [modelSub.rxns{idx(i)},'_r'];
%     end
end
if lumpIdenticalRxns == 1
    keepRxns = find(ismember(find(subRxns),keepRxns)); 
    
    [modelSub.S,vSub,rMap] = combineparallel(modelSub.S,vSub,0,keepRxns);
    rxns = cell(length(vSub),1);
    [c,lb,ub] = deal(zeros(length(vSub),1));
    for i = 1:length(vSub)
        rIdx = find(rMap(i,:));
        mIdx = find(modelSub.S(:,i));
        if length(mIdx) == 1 && length(rIdx) > 1 %lumped exchange of single metabolite
            if modelSub.S(mIdx,i) < 0
                rxns{i} = [modelSub.mets{mIdx},'_Sink'];
            else
                rxns{i} = [modelSub.mets{mIdx},'_Source'];
            end
        elseif length(rIdx)>3
            reactIdx = find(modelSub.S(:,i)<0);
            prodIdx = find(modelSub.S(:,i)>0);
            rxns{i} = modelSub.mets{reactIdx(1)};
            for j = 2:length(reactIdx)
               rxns{i} = [rxns{i}, '_and_', modelSub.mets{reactIdx(j)}];
            end
            rxns{i} = [rxns{i}, '_to_', modelSub.mets{prodIdx(1)}];
            for j = 2:length(prodIdx)
                rxns{i} = [rxns{i}, '&', modelSub.mets{prodIdx(j)}];
            end 
        else
            rxns{i} = modelSub.rxns{rIdx(1)};
            for j = 2:length(rIdx)
                rxns{i} = [rxns{i},'--',modelSub.rxns{rIdx(j)}];
            end
        end
        c(i) = max(modelSub.c(rIdx));
        lb(i) = max(sum(modelSub.lb(rIdx)),-1000);
        ub(i) = max(sum(modelSub.ub(rIdx)),1000);
    end
    modelSub.rxns = rxns;
    modelSub.rxnNames = rxns;
    modelSub.c = c;
    modelSub.lb = lb;
    modelSub.ub = ub;
end
modelSub.rev = any(vSub<0,2);

if exist('filename','var')
    exportModeltoExcel(modelSub.S,modelSub.mets,vSub,modelSub.rxns,modelSub.rev,filename);
end
end
function [Sp,vp,rMap,rMat] = combineparallel(S,v,antiflag,keepRxns)
% Combines identical and reverse reactions in the stoichiometric matrix
% Inputs are S (the stoichiometric matrix) and v (the vector of fluxes)
% Outputs are Sc and vc, the stoichiometric matrix and flux vector with
% identical and reverse reactions combined.

% Modified by Sean Mack 11/26/19 - removed duplicate parallel/antiparallel
% function calls
% Modified by Sean Mack 12/20/19 - added parallel pathways as output for
% tracking model reduction
% Modified by Sean Mack 1/17/20 - added capability to handle multiple flux
% vectors
% Modified by Sean Mack 2/14/20 - adjusted code to run through all
% parallel pairs at once and accurately report original index pairs
% Modified by Sean Mack 2/24/20 - revamped code to speed up and report full
% sets of parallel and antiparallel reactions
if nargin < 3
    antiflag = 1;
end
if nargin < 4
    keepRxns = [];
end

vp = v;

[parPair(:,1),parPair(:,2)] = parallel(S);

keepPPairs = any(ismember(parPair,keepRxns),2);
parPair = parPair(~keepPPairs,:);

rMap = speye(size(S,2));
for i = 1:size(parPair,1)
    p = parPair(i,:);
    if any(rMap(p(1),:))
        rMap(p(2),:) = rMap(p(2),:) + rMap(p(1),:);
        rMap(p(1),:) = 0;
        f = norm(S(:,p(2)))/norm(S(:,p(1)));
        for j = 1:size(vp,2)
            vp(p(2),j) = vp(p(2),j) + vp(p(1),j)/f;
        end
        vp(p(1),:) = 0;
    end
end
if antiflag
    [revPair(:,1),revPair(:,2)] = antiparallel(S);
    keepRPairs = any(ismember(revPair,keepRxns),2);
    revPair = revPair(~keepRPairs,:);
    for i = 1:size(revPair,1)
        a = revPair(i,:);
        if all(any(rMap(a,:),2))
            rMap(a(2),:) = rMap(a(2),:) - rMap(a(1),:);
            rMap(a(1),:) = 0;
            f = norm(S(:,a(2)))/norm(S(:,a(1)));
            for j = 1:size(vp,2)
                vp(a(2),j) = vp(a(2),j) - vp(a(1),j)/f;
            end
            vp(a(1),:) = 0;
        end
    end
end
rMat = speye(size(rMap));
tf = ~any(rMap,2);
rMat(tf,:) = [];
rMap(tf,:) = [];
Sp = S*rMat';
vp = rMat*vp;
end
function [r,c] = parallel(A)
%PARALLEL Determines the location of parallel vectors in the column space
%of a matrix
An = zeros(size(A));
for i = 1:size(A,2)
    if sum(A(:,i).^2) ~= 0
        An(:,i) = A(:,i)/norm(A(:,i));
    end
end
R = An'*An;
R = triu(R);
R = R - eye(size(R,1),size(R,2));
[r,c] = find(R >= 0.99999999);

end
function [r,c] = antiparallel(A)
%ANTIPARALLEL Determines the location of antiparallel vectors in the column space
%of a matrix
An = zeros(size(A));
for i = 1:size(A,2)
    if sum(A(:,i).^2) ~= 0
        An(:,i) = A(:,i)./sqrt(sum(A(:,i).^2));
    end
end
R = An'*An;
R = triu(R);
R = R - eye(size(R,1),size(R,2));
[r,c] = find(R <= -0.99999999);

end
function [] = exportModeltoExcel(S,mets,flux,rxns,rev,filename)
if isempty(rxns)
    rxns = (1:size(S,2))';
end
if isempty(rev)
    rev = zeros(size(S,2),1);
end

s = {' '};
arrowI = '->';
arrowR = '<->';
plus = '+';
S = full(S);

reaction = cell(size(S,2),1);
for j = 1:size(S,2)
    
    reactants = [];
    products = [];
    [u,~] = find(S(:,j) < 0);
    [a,~] = find(S(:,j) > 0);
    
    for i = 1:length(u)
        reactants = strcat(reactants,num2str(abs(S(u(i),j))),s,mets(u(i)),s);
        if i ~= length(u)
            reactants = strcat(reactants,s,plus,s);
        end
    end
    
    for k = 1:length(a)
        products = strcat(products,num2str(S(a(k),j)),s,mets(a(k)),s);
        if k ~= length(a)
            products = strcat(products,s,plus,s);
        end
    end
    if rev(j) == 1
        reaction{j,1} = strcat(reactants,s,arrowR,s,products);
    else
        reaction{j,1} = strcat(reactants,s,arrowI,s,products);
    end
end

reactionChar = cell(size(S,2),1);
for i = 1:length(reaction)
    reactionChar{i,1} = char(reaction{i});
end

warning('off','MATLAB:xlswrite:AddSheet')
xlswrite(filename,rxns,'Reactions','A1');
xlswrite(filename,reactionChar,'Reactions','B1');
xlswrite(filename,flux,'Reactions','C1');
xlswrite(filename,mets,'Metabolites','A1');
end