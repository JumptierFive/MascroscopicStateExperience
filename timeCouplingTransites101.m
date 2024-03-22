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
allBinTr={};
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
    
    binTr=zeros(length(binM),1);
    for mi=2:1:numel(binM)
        if binM(mi-1)~=binM(mi)
            binTr(mi)=1;
        end      
    end
    
    
    binT=binM;binT(binM==34)=1;binT(binM==136)=1;binT(binM==121)=2;binT(binM==223)=2;binT(binM==82)=3;binT(binM==175)=4;
    allBinC{subid,1}=binC;allBinM{subid,1}=binM;allBinT{subid,1}=binT;allBinTr{subid,1}=binTr;
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
    binTr=zeros(length(binM),1);
    for mi=2:1:numel(binM)
        if binM(mi-1)~=binM(mi)
            binTr(mi)=1;
        end
    end
    
    binT=binM;binT(binM==34)=1;binT(binM==136)=1;binT(binM==121)=2;binT(binM==223)=2;binT(binM==82)=3;binT(binM==175)=4;
    allBinC{subid,2}=binC;allBinM{subid,2}=binM;allBinT{subid,2}=binT;allBinTr{subid,2}=binTr;
end
info='this is the network timeseries for whole resting state brain, allBinC for the states calculated by pairwise MEM method, allBinM for all basins in energy landscape, allBinT for all large basins in tree structure analysis';
realStateTS.allBinC=allBinC;realStateTS.allBinM=allBinM;realStateTS.allBinT=allBinT;realStateTS.info=info;
save('stateSeries.mat','realStateTS');

allPor1=[];
memoPor1=[];
for subid=1:1:24
    for sid=1
        similarity=(mean(leftRestSimilarity{subid,sid},2)+mean(rightRestSimilarity{subid,sid},2))./2;
        xx=find(similarity>0.5);
        tpBT=allBinTr{subid,sid}(xx);
        x1=find(allBinTr{subid,sid}==1);x2=find(allBinTr{subid,sid}==2);x3=find(allBinTr{subid,sid}==3);x4=find(allBinTr{subid,sid}==4);x5=find(allBinTr{subid,sid}==5);
        x6=find(allBinTr{subid,sid}==6);x7=find(allBinTr{subid,sid}==7);x8=find(allBinTr{subid,sid}==8);x9=find(allBinTr{subid,sid}==9);x10=find(allBinTr{subid,sid}==10);
        ctX=[numel(x1),numel(x2),numel(x3),numel(x4),numel(x5),numel(x6),numel(x7),numel(x8),numel(x9),numel(x10)]./numel(allBinTr{subid,sid});
        allPor1=[allPor1;ctX];
        x1=find(tpBT==1);x2=find(tpBT==2);x3=find(tpBT==3);x4=find(tpBT==4);x5=find(tpBT==5);x6=find(tpBT==6);x7=find(tpBT==7);x8=find(tpBT==8);x9=find(tpBT==9);x10=find(tpBT==10);
        ctX=[numel(x1),numel(x2),numel(x3),numel(x4),numel(x5),numel(x6),numel(x7),numel(x8),numel(x9),numel(x10)]./numel(tpBT);
        memoPor1=[memoPor1;ctX];
    end
end


allPor2=[];
memoPor2=[];
for subid=1:1:24
    for sid=2
        similarity=(mean(leftRestSimilarity{subid,sid},2)+mean(rightRestSimilarity{subid,sid},2))./2;
        xx=find(similarity>0.5);
        tpBT=allBinTr{subid,sid}(xx);
        x1=find(allBinTr{subid,sid}==1);x2=find(allBinTr{subid,sid}==2);x3=find(allBinTr{subid,sid}==3);x4=find(allBinTr{subid,sid}==4);x5=find(allBinTr{subid,sid}==5);
        x6=find(allBinTr{subid,sid}==6);x7=find(allBinTr{subid,sid}==7);x8=find(allBinTr{subid,sid}==8);x9=find(allBinTr{subid,sid}==9);x10=find(allBinTr{subid,sid}==10);
        ctX=[numel(x1),numel(x2),numel(x3),numel(x4),numel(x5),numel(x6),numel(x7),numel(x8),numel(x9),numel(x10)]./numel(allBinTr{subid,sid});
        allPor2=[allPor2;ctX];
        x1=find(tpBT==1);x2=find(tpBT==2);x3=find(tpBT==3);x4=find(tpBT==4);x5=find(tpBT==5);x6=find(tpBT==6);x7=find(tpBT==7);x8=find(tpBT==8);x9=find(tpBT==9);x10=find(tpBT==10);
        ctX=[numel(x1),numel(x2),numel(x3),numel(x4),numel(x5),numel(x6),numel(x7),numel(x8),numel(x9),numel(x10)]./numel(tpBT);
        memoPor2=[memoPor2;ctX];
    end
end

allPor=allPor1+allPor2;
memoPor=memoPor1+memoPor2;