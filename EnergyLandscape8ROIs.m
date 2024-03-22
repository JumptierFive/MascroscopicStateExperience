%% estimation energy landscape for each subject
clear all;
clc;
threshold =0.0; %for binarization, above (below) which ROI activity is defined to be +1 (-1).
rootDir=['D:\dataN\ReplayFMRI\wholeTimeSeriesS01\'];
allFiles=filename_list(rootDir,'sub_*.mat');
load(allFiles{1})
roiN=length(subTS(1,:));
values=[-1,1];
allStates=generate_states(roiN);
stateNet=getNetStates(allStates);

EnergyLandscape={};
EnergyLandscape.stateNet=stateNet;
EnergyLandscape.allStates=allStates;
allBin=[];
% rootDir=['D:\dataN\ReplayFMRI\wholeTimeSeriesS01\'];
% allFiles=filename_list(rootDir,'sub_*DMN.mat');
% for subid=1:1:length(allFiles)
%     load(allFiles{subid});
%     subTS=subTS';
%     binarizedData = pfunc_01_Binarizer(subTS,threshold);
%     allBin=[allBin,binarizedData];
% end

rootDir=['D:\dataN\ReplayFMRI\wholeTimeSeriesS02\'];
allFiles=filename_list(rootDir,'sub_*DMN.mat');
for subid=1:1:length(allFiles)
    load(allFiles{subid});
    subTS=subTS';
    binarizedData = pfunc_01_Binarizer(subTS,threshold);
    allBin=[allBin,binarizedData];
end

EnergyLandscape.allBins=allBin;
disp('reading time series done');

[h,J] = pfunc_02_Inferrer_ML(allBin);
[probN, prob1, prob2, rD, r] = pfunc_03_Accuracy(h, J, allBin);
EnergyLandscape.h=h;
EnergyLandscape.J=J;
EnergyLandscape.probAll=[probN,prob1,prob2];
EnergyLandscape.rD=rD;
EnergyLandscape.r=r;

disp('pairwise MEM model done');


tempEnergy=[];
allstatesN=length(allStates(:,1));
for i=1:1:allstatesN
    energy = energyCalculate(h,J,allStates(i,:));
    tempEnergy=[tempEnergy;energy];
end
EnergyLandscape.EnergyNode=tempEnergy;
disp('enegy for each node calculation done');

localMin=[];%local mins in energy landcape as a typical states in neural dynamics
for i=1:1:allstatesN
    xx=find(stateNet(i,:)==1);
    if sum(tempEnergy(i)>tempEnergy(xx))==0
        localMin=[localMin,i];
    end
end
EnergyLandscape.localMins=localMin;

allSets=[]; %tree organziation for the energy landscape
ctEnergy=tempEnergy;
ctNet=stateNet;
for i=1:1:1000000
    if sum(sum(ctNet))==0
        break;
    end
    aa=max(ctEnergy);
    xx=find(ctEnergy==aa);
    ctEnergy(xx)=-1000;
    ctNet(xx,:)=0;ctNet(:,xx)=0;
    G = graph(ctNet);
    bins=conncomp(G);
    allSets=[allSets;bins(localMin)];
    
end
EnergyLandscape.UnconnectedGraph=allSets;

lN=numel(localMin);
barrier=zeros(lN,lN);% barriers between any two local mins 
ppath={};
G=graph(stateNet);
for ki=1:1:lN
    for ji=1:1:lN
        if ki==ji
            continue
        end
        allPath = pathof(G,localMin(ki),localMin(ji));
        energyL=[];
        for pi=1:1:length(allPath)
            energyL=[energyL;max(tempEnergy(allPath{pi}(2:end)))];
        end
        barrier(ki,ji)=min(energyL);
    end
end
EnergyLandscape.barrier=barrier;

BasinLabels=[];% basins for enegylanscape
for i=1:1:allstatesN
    [basinN] = nodeBasin(stateNet,i,localMin,tempEnergy);
    BasinLabels(i)=basinN;
end
EnergyLandscape.basinLabel=BasinLabels;
disp('energy landscape done')
save('EnergyLandscapeS02.mat','EnergyLandscape');

%% belows are functions refferring to energy landscapes model calculations


function paths = findPaths(adjM, startN, endN, n)
%% the function locates paths between startN node and endN node within n steps
paths={};poolPathC={startN};
for zi=1:1:n
    newPoolPath={};
    for ki=1:1:length(poolPathC)
        % find current node's neighbors and judge if they had been reached
        xx=find(adjM(poolPathC{ki}(end),:)==1);
        xx=setdiff(xx,poolPathC{ki});
        if isempty(xx)
            continue;
        end
        
        for ji=1:1:length(xx)
            newPoolPath{end+1}=[poolPathC{ki}(:)',xx(ji)];
            if xx(ji)==endN
                paths{end+1}=[poolPathC{ki}(:)',endN];
            end
        end
    end
    poolPathC=newPoolPath;
    length(poolPathC)
end
end





function [basinN] = nodeBasin(adjM, startN, reachSet, Energy)
%% the function locates paths between startN node and endN node within n steps
ctNode=startN;
for ki=1:1:10000000000
    ctx=find(reachSet==ctNode);
    if ~isempty(ctx)
        basinN=reachSet(ctx);
        break;
    end
    xx=find(adjM(ctNode,:)==1);
    ctE=Energy(xx);
    ex=find(ctE==min(ctE));
    ctNode=xx(ex);
end
end





function all_states=generate_states(roiN)
all_states=zeros(2^roiN,roiN);
for i=1:1:2^roiN
    binary_representation = de2bi(i-1,roiN);
    all_states(i,:)=2*binary_representation-1;
end
end



function energy = energyCalculate(h,J,state)
% Ising model canculate energy
fieldEnergy=-1*sum(state.*h');
interEnergy=-1*sum(sum(state'*state.*J))/2;
energy=fieldEnergy+interEnergy;
end

function [statesNet] = getNetStates(States)
netN=length(States(:,1));
statesNet=zeros(netN,netN);
for i=1:1:netN
    diffState=sum(abs(States-States(i,:)),2);
    xx=find(diffState==2);
    statesNet(i,xx)=1;
end
end

function pth=pathof(graph,startn,endn)
stop=0;
n=0;
while stop~=1
    n=n+1;
    Temp=shortestpath(graph,startn,endn);
    eidx=findedge(graph,Temp(1:end-1),Temp(2:end));
    if n~=1
        if length(Temp)==length(pth{n-1,1})
            if Temp==pth{n-1,1}
                stop=1;
            else
                pth{n,1}=Temp;
                graph.Edges.Weight(eidx)=100;
            end
        else
            pth{n,1}=Temp;
            graph.Edges.Weight(eidx)=100;
        end
    else
        pth{n,1}=Temp;
        graph.Edges.Weight(eidx)=100;
    end
    clear Temp eidx;
end
end