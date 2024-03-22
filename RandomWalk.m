%% random walk simulation on the energy landscape
clear all;
clc;
load('EnergyLandscape.mat');


aaaaaa

% calculate for transfer probility
transferProb=EnergyLandscape.stateNet;
[xx,~]=size(EnergyLandscape.stateNet);
for i=1:1:xx
    for j=1:1:xx
        if EnergyLandscape.stateNet(i,j)==1
            transferProb(i,j)=min(1,exp(EnergyLandscape.EnergyNode(i)-EnergyLandscape.EnergyNode(j)));
        end
    end
end
ctS=sum(transferProb,2);
transferProbM=[];
for ki=1:1:xx
    transferProbM=[transferProbM;transferProb(ki,:)./ctS(ki)];
end
KM=sum(transferProbM,2);
transferS={};
for i=1:1:xx
    prob=transferProbM(i,:);xn=find(prob>0);prob(prob==0)=[];
    transferX{i}.prob=cumsum(prob);
    transferX{i}.nei=xn;
end
% random walk on the energy landscape
allWalksSim={};
parfor i=1:10
    i
    aa=randperm(length(transferX));ctN=aa(1);ctWalk=[];ctWM=[];
    for ki=1:1:10000
        ctWalk=[ctWalk;ctN];
        ctWM=[ctWM;EnergyLandscape.basinLabel(ctN)];
        ctM=rand(1,1);
        for mi=length(transferX{ctN}.prob):-1:1
            if ctM>transferX{ctN}.prob(mi)
                ctN=transferX{ctN}.nei(mi);
                break;
            end
        end
    end
    binT=ctWM;binT(binT==34)=1;binT(binT==136)=1;binT(binT==121)=2;binT(binT==223)=2;binT(binT==82)=3;binT(binT==175)=4;
    allWalksSim{i}.allStates=ctWalk;allWalksSim{i}.Mins=ctWM;allWalksSim{i}.Basin=binT;
end
save('randomWalksSimulation.mat','allWalksSim');
