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

%% 121 223  136  34
tar1=34;tar2=121;
allTM=0;allTMP=0;
for i=1:1:24
    for si=1:1:2
        tpAN=allBinM{i,si};
        tpN=tpAN(1);
        for ki=1:1:length(tpAN)-1
            if tpAN(ki)~=tpAN(ki+1)
                tpN=[tpN,tpAN(ki+1)];
            end
        end
        tm=0;tmp=0;
        for mi=1:1:length(tpN)-1
            if tpN(mi)==tar1
                tm=tm+1;
                if tpN(mi+1)==tar2
                    tmp=tmp+1;
                end
            end
        end
        allTM=allTM+tm;allTMP=allTMP+tmp;
    end
end
allTMP./allTM