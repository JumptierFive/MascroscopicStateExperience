%% time coupling between memory replay and whole brain dynamics
clear all;
clc
load('EnergyLandscape.mat');
load('leftRestSimilarityLof.mat');
load('rightRestSimilarityLof.mat');

rootDir=['D:\dataN\ReplayFMRI\wholeTimeSeriesS01\'];
allFiles=filename_list(rootDir,'sub_*DMN.mat');
threshold =0.0; %for binarization, above (below) which ROI activity is defined to be +1 (-1).
allBinC={};
allBinM={};
allBinT={};
AN=121;BN=136;
for subid=1:1:length(allFiles)
    load(allFiles{subid});
    subTS=subTS';
    binarizedData = pfunc_01_Binarizer(subTS,threshold);
    binC=[];
    binM=[];
    for i=1:1:length(binarizedData(1,:))
        index=find(ismember(EnergyLandscape.allStates,binarizedData(:,i)','rows'));
        bas=EnergyLandscape.basinLabel(index);
        binC=[binC;index];binM=[binM;bas];
    end
    binT=binM;binT(binM==121)=1;binT(binM==136)=1;binT(binM==121)=2;binT(binM==136)=2;binT(binM==82)=3;binT(binM==175)=4;
    allBinC{subid,1}=binC;allBinM{subid,1}=binM;allBinT{subid,1}=binT;
end

rootDir=['D:\dataN\ReplayFMRI\wholeTimeSeriesS02\'];
allFiles=filename_list(rootDir,'sub_*DMN.mat');
for subid=1:1:length(allFiles)
    load(allFiles{subid});
    subTS=subTS';
    binarizedData = pfunc_01_Binarizer(subTS,threshold);
    binC=[];
    binM=[];
    for i=1:1:length(binarizedData(1,:))
        index=find(ismember(EnergyLandscape.allStates,binarizedData(:,i)','rows'));
        bas=EnergyLandscape.basinLabel(index);
        binC=[binC;index];binM=[binM;bas];
    end
    binT=binM;binT(binM==121)=1;binT(binM==136)=1;binT(binM==121)=2;binT(binM==136)=2;binT(binM==82)=3;binT(binM==175)=4;
    allBinC{subid,2}=binC;allBinM{subid,2}=binM;allBinT{subid,2}=binT;
end
info='this is the network timeseries for whole resting state brain, allBinC for the states calculated by pairwise MEM method, allBinM for all basins in energy landscape, allBinT for all large basins in tree structure analysis';
realStateTS.allBinC=allBinC;realStateTS.allBinM=allBinM;realStateTS.allBinT=allBinT;realStateTS.info=info;
save('stateSeries.mat','realStateTS');

allPor1=[];
memoPor1=[];
for subid=1:1:24
    for sid=1
        numel(allBinM{subid,sid})
        allBinM{subid,sid}=allBinM{subid,sid}(1:end-1);
        similarity=(mean(leftRestSimilarity{subid,sid},2)+mean(rightRestSimilarity{subid,sid},2))./2;
        similarity=similarity(2:end);
        xx=find(similarity>0.5);
        tpBT=allBinM{subid,sid}(xx);
        x1=find(allBinM{subid,sid}==AN);x2=find(allBinM{subid,sid}==BN);
        x1=numel(x1)./numel(allBinM{subid,sid});x2=numel(x2)./numel(allBinM{subid,sid});
%         x1=numel(x1);x2=numel(x2);
        allPor1=[allPor1;[x1,x2]];
        x1=find(tpBT==AN);x2=find(tpBT==BN);
        x1=numel(x1)./numel(tpBT);x2=numel(x2)./numel(tpBT);
        memoPor1=[memoPor1;[x1,x2]];
    end
end


allPor2=[];
memoPor2=[];
for subid=1:1:24
    for sid=2
        numel(allBinM{subid,sid})
        allBinM{subid,sid}=allBinM{subid,sid}(1:end-2);
        similarity=(mean(leftRestSimilarity{subid,sid},2)+mean(rightRestSimilarity{subid,sid},2))./2;
        similarity=similarity(3:end);
        xx=find(similarity>0.5);
        tpBT=allBinM{subid,sid}(xx);
        x1=find(allBinM{subid,sid}==AN);x2=find(allBinM{subid,sid}==BN);
        x1=numel(x1)./numel(allBinM{subid,sid});x2=numel(x2)./numel(allBinM{subid,sid});
%         x1=numel(x1);x2=numel(x2);
        allPor2=[allPor2;[x1,x2]];
        x1=find(tpBT==AN);x2=find(tpBT==BN);
        x1=numel(x1)./numel(tpBT);x2=numel(x2)./numel(tpBT);
        memoPor2=[memoPor2;[x1,x2]];
    end
end

allPor=(allPor1+allPor2)./2;
memoPor=(memoPor1+memoPor2)./2;