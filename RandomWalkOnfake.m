%% random walk simulation on the energy landscape
clear all;
clc;
load('EnergyLandscape.mat');


allWalksSim={};
for i=1:1000
    i
    
    J=EnergyLandscape.J;
    for mi=1:1:length(J)
        for hi=1:1:length(J)
            if mi<=hi
                continue;
            end
            permT=normrnd(0,0.05,1,1);
            J(mi,hi)=J(mi,hi)+permT;J(hi,mi)=J(hi,mi)+permT;
        end
    end
    allWalksSim{i}.J=J;
    
    tempEnergy=[];
    allstatesN=length(EnergyLandscape.allStates(:,1));
    for zi=1:1:allstatesN
        energy = energyCalculate(EnergyLandscape.h,J,EnergyLandscape.allStates(zi,:));
        tempEnergy=[tempEnergy;energy];
    end
    EnergyLandscape.EnergyNode=tempEnergy;
    
    % calculate for transfer probility
    transferProb=EnergyLandscape.stateNet;
    [xx,~]=size(EnergyLandscape.stateNet);
    for ii=1:1:xx
        for jj=1:1:xx
            if EnergyLandscape.stateNet(ii,jj)==1
                transferProb(ii,jj)=min(1,exp(EnergyLandscape.EnergyNode(ii)-EnergyLandscape.EnergyNode(jj)));
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
    for ti=1:1:xx
        prob=transferProbM(ti,:);xn=find(prob>0);prob(prob==0)=[];
        transferX{ti}.prob=cumsum(prob);
        transferX{ti}.nei=xn;
    end
    % random walk on the energy landscape

    
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
save('randomWalksSimulationFakeNew005.mat','allWalksSim');



function energy = energyCalculate(h,J,state)
% Ising model canculate energy
fieldEnergy=-1*sum(state.*h');
interEnergy=-1*sum(sum(state'*state.*J))/2;
energy=fieldEnergy+interEnergy;
end
