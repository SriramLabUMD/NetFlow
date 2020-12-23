%% Script for generating the data from Mack and Sriram, Met. Eng. Comm. (2020)
clear
clc
%% Run toy model with NetFlow
load('toyMFA.mat')

vToy = [10 10 5 5 3 2 8 10 13 -3]'; %toy v1

[flow] = NetFlow(toy,vToy);

%Identify metabolite sets and create submodels
toE = flow.metLinks(:,5);
[toy2E,vToE] = extractSubnetwork(toy,vToy,toE,'Toy Subnetwork - toE.xlsx');

B2D = flow.metLinks(2,:) & flow.metLinks(:,4)';
[toyB2D,vB2D] = extractSubnetwork(toy,vToy,B2D,'Toy Subnetwork - B2D.xlsx');

save('toySubs.mat','toy2E','toyB2D')

%% Case Study 1: CCM Model
load('iDM2014_CCM.mat')

% Define media conditions
importRxns = find(strcmp(model.subSystems,'Uptake'));
importMets = cell2mat(model.prodIdx(importRxns));
exportRxns = find(strcmp(model.subSystems,'Export'));
exportMets = cell2mat(model.reactIdx(exportRxns));
model.ub(importRxns) = [0;10]; %[CO2;Glc]

% Identify reversible reaction pairs
rev = [0,0];
nRxns = length(model.rxns);
for i = 1:nRxns
    for j = 1:nRxns
        tmp = ~any(model.S(:,i) + model.S(:,j));
        if tmp
            if all(rev(:,2) ~= i)
                rev = [rev;i,j];
            end
        end
    end
end
rev = rev(2:end,:);

% Define flux ratios
model.ub(42) = 10; %Glc uptake

edRxns = strcmp(model.subSystems,'ED');
anapRxns = strcmp(model.subSystems,'Anap');
model.ub(edRxns) = 0; %Block ED activity
model.ub(anapRxns) = 0; %Block anaplerotic activity
model.ub([63,83,71,69]) = 0; %Force irrevesible flux by blocking the reverse reactions GLCt2pp, PGI_rev, LDH_rev, & ICDH_rev 

empVppp = [82,66]; %PGI (EMP) & GND (PPP) 
tcaVferm = [79,81,72]; %PDH (TCA), PFL (TCA/For), & LDH (Lac)
glxVtca = [68,70]; %ICDH (TCA) & ICL (Glx)

empRatio = [8,2]; %EMP:PPP ratio defined as 4:1
tcaRatio = [3,1,1]; %TCA to fermentation ratio defined as 3:1
glxRatio = [1,1]; %Glx:TCA ratio defined as 1:1

%EMP and PPP
model.lb(empVppp) = empRatio;

%TCA and Ferm
model.c(:) = 0;
model.c(tcaVferm(1)) = 1;
sol = optimizeCbModel(model);
tmpV = sol.f;
model.lb(tcaVferm) = tmpV.*(tcaRatio./5);

%Glx and TCA
model.c(:) = 0;
model.c(glxVtca(2)) = 1;
sol = optimizeCbModel(model);
tmpV = sol.f;
model.lb(glxVtca) = tmpV.*(glxRatio./2);

%Base Efflux
model.c(:) = 0;
model.c(45) = 1; %Succinate efflux
sol = optimizeCbModel(model);
vCCM(:,1) = sol.x;

%LDH KO
tmpUB = model.ub; tmpLB = model.lb;
model.ub(72) = 0; model.lb(72) = 0;
sol = optimizeCbModel(model);
vCCM(:,2) = sol.x;
model.ub = tmpUB; model.lb = tmpLB;

%ICDH KO
model.ub(68) = 0; model.lb(68) = 0;
sol = optimizeCbModel(model);
vCCM(:,3) = sol.x;

%LDH & ICDH KO
model.ub(72) = 0; model.lb(72) = 0;
sol = optimizeCbModel(model);
vCCM(:,4) = sol.x;
model.ub = tmpUB; model.lb = tmpLB;

% Find net fluxes for reversible
for i = 1:length(rev)
    v = (vCCM(rev(i,1),:)-vCCM(rev(i,2),:));
    neg = v < 0;
    vCCM(rev(i,2),neg) = abs(v(neg));
    vCCM(rev(i,2),~neg) = 0;
    vCCM(rev(i,1),~neg) = v(~neg);
    vCCM(rev(i,1),neg) = 0;
end

% Run NetFlow
[base] = NetFlow(model,vCCM(:,1));
[ldhKO] = NetFlow(model,vCCM(:,2));
[icdhKO] = NetFlow(model,vCCM(:,3));
[doubleKO] = NetFlow(model,vCCM(:,4));

% Define reactants in KO reactions for pathway isolation
ldhReact = 47; %pyr
icdhReact = 40; %icit
bothReact = [47,40]; %pyr & icit

zFlux = all(vCCM == 0,2);
modelR = removeRxns(model,model.rxns(zFlux), 'metFlag', false);
vAct = vCCM(~zFlux,:);

% Identify metabolites contributing carbon to succinate
succ.links(:,1) = base.metLinks(:,54);
succ.links(:,2) = ldhKO.metLinks(:,54);
succ.links(:,3) = icdhKO.metLinks(:,54);
succ.links(:,4) = doubleKO.metLinks(:,54);

% Identify metabolite yield on succinate
succ.yield(:,1) = base.metYield(:,54);
succ.yield(:,2) = ldhKO.metYield(:,54);
succ.yield(:,3) = icdhKO.metYield(:,54);
succ.yield(:,4) = doubleKO.metYield(:,54);

% Identify metabolites connecting KO reactants and succinate
KO.links(:,1) = (base.metLinks(ldhReact,:) | ldhKO.metLinks(ldhReact,:))' & (succ.links(:,2) | succ.links(:,1));
KO.links(:,2) = (base.metLinks(icdhReact,:) | icdhKO.metLinks(icdhReact,:))' & (succ.links(:,3) | succ.links(:,1));
KO.links(:,3) = any((base.metLinks(bothReact,:) | doubleKO.metLinks(bothReact,:))',2) & (succ.links(:,4) | succ.links(:,1));

% Extract subnetworks
options.lumpIdenticalRxns = 0;
[ccmLinks,vLinks] = extractSubnetwork(modelR,vAct,succ.links(:,1),options);

options.keepRxns = find(ismember(modelR.rxns,model.rxns(72)));
[ldhLinks,vLDH] = extractSubnetwork(modelR,vAct,KO.links(:,1),options);

options.keepRxns = find(ismember(modelR.rxns,model.rxns(68)));
[icdhLinks,vICDH] = extractSubnetwork(modelR,vAct,KO.links(:,2),options);

options.keepRxns = find(ismember(modelR.rxns,model.rxns([68,72])));
[doubleLinks,vDouble] = extractSubnetwork(modelR,vAct,KO.links(:,3),options);

save('ccmSubs.mat','ccmLinks','ldhLinks','icdhLinks','doubleLinks')

%% Case Study 2: Apler, et al Knockout Analysis
load('alperKO.mat')
load('iDM2014_lycop.mat')

% Define media conditions
[model] = defineEcoliMediaAlper(model);
BM = find(model.c);
lycopEx = 550; %Lycopene exchange

% Define rxns to KO
for i = 1:length(alper.rxns)
    rIdx{i,1} = find(ismember(model.rxns,alper.rxns{i}));
end
rIdx{6} = [rIdx{2},rIdx{3}]; %gdhA & aceE
rIdx{7} = [rIdx{2},rIdx{5}]; %gdhA & talB
rIdx{8} = [rIdx{2},rIdx{1}]; %gdhA & fdhF
rIdx{9} = [rIdx{2},rIdx{3},rIdx{5}]; %gdhA, aceE, & talB
rIdx{10} = [rIdx{2},rIdx{3},rIdx{1}]; %gdhA, aceE, & fdhF

% Simulate all KOs with MOMA
nKO = length(rIdx);
vKO = zeros(length(model.rxns),nKO);
[lycopene,growth] = deal(zeros(nKO,1));
for i = 1:nKO
    modelKO(i,1) = model;
    modelKO(i).ub(rIdx{i}) = 0;
    modelKO(i).lb(rIdx{i}) = 0;
    [solDel, solWT, totalFluxDiff, solStat] = MOMA(model, modelKO(i),'max',0,1);
    vKO(:,i) = solDel.x;
    lycopene(i,1) = solDel.x(lycopEx);
    growth(i,1) = solDel.x(BM);
    
    %Eliminate infeasible loop of ADK4 and NTP10
    if vKO(58,i) == 0 %ATPHs
        vKO([24,412],i) = 0; %ADK4 & NTP10 
    end
end
vWT = solWT.x;
growth(:,2) = growth(:,1)/solWT.x(BM);

% Run fluxes through NetFlow
%Run WT flux
[flowWT] = NetFlow(model,vWT);

% Analyze KO fluxes
KOtoRun = [2,6,10]; %gdhA; gdhA & aceE; gdhA, aceE, & fdhF
KOreact = [214,369,194]; %Glu, Pyr, & For
KOname = {'Single','Double','Triple'};
for i = 1:length(KOtoRun)
    %Run NetFlow on the corresponsing KO flux
    idx = KOtoRun(i);
    [flowKO(i)] = NetFlow(modelKO(idx),vKO(:,idx));
    
    %Identify metabolite chain
    refMet = 431; %Lycopene
    
    %Identify carbon-contributing metabolites and carbon yield on lycopene
    lycop(i).links = flowKO(i).metLinks(:,refMet);
    lycop(i).yield = flowKO(i).metYield(:,refMet);
    
    %Identify all carbon-containing metabolites participating in the KO
    % reaction(s)
    KO(i).mets = find(any(model.S(:,rIdx{idx}),2) & model.nCarbon>0);
    KO(i).react = KOreact(1:i); %reactants only
    
    %Identify metabolites recieving carbon from metabolites involved in
    %KO reaction(s)
    KO(i).links = any(flowKO(i).metLinks(KO(i).mets,:),1);
    KO(i).rLinks = any(flowKO(i).metLinks(KO(i).react,:),1); %reactants only
    
    %Extract subnetworks
    vRef = [vWT,vKO(:,idx)];
    zFlux = all(vRef == 0,2);
    modelR = removeRxns(model,model.rxns(zFlux), 'metFlag', false);
    vRef = vRef(~zFlux,:);
    
    %Set options to lump identical reactions in subnetworks but ensure KO
    %reaction(s) are kept seperate
    options.lumpIdenticalRxns = 1;
    options.keepRxns = find(ismember(modelR.rxns,model.rxns(rIdx{idx})));
    
    %Lycopene links
%     options.filename = sprintf('Alper %s KO - Lycopene Links.xlsx',KOname{i});
    [modelLycL(i),vLycL{i}] = extractSubnetwork(modelR,vRef,lycop(i).links,options);
    
    %Lycopene yield
%     options.filename = sprintf('Alper %s KO - Lycopene Yield.xlsx',KOname{i});
    yieldMets(:,i) = (lycop(i).yield>0.001);
    [modelLycY(i),vLycY{i}] = extractSubnetwork(modelR,vRef,yieldMets(:,i),options);
    
    %Lycopene yield and KO reactant(s)
%     options.filename = sprintf('Alper %s KO - Lycopene Yield & KO Mets.xlsx',KOname{i});
    yieldMets(KO(i).react,i) = 1;
    [modelLYK(i),vLYK{i}] = extractSubnetwork(modelR,vRef,yieldMets(:,i),options);
    
    %Lycopene and KO reactant links
%     options.filename = sprintf('Alper %s KO - Lycopene Yield & KO React Link.xlsx',KOname{i});
    bothMets(:,i) = yieldMets(:,i) & KO(i).rLinks';
    [modelBoth(i),vBoth{i}] = extractSubnetwork(modelR,vRef,bothMets(:,i),options);
end
save('alperSubs.mat','modelLycL','modelLycY','modelLYK','modelBoth')
