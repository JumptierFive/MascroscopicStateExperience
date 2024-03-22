%% Permutation test for resting-state scans similar to each figure item activtion map
clear all;
clc;
rootDir='D:\dataN\';
template=spm_read_vols(spm_vol('LeftTemplate.nii'));template=reshape(template,[1,53*63*52]);
xx=find(template>0);
load('LeftHippVols.mat');
load('leftLofCleared.mat');
load('leftPermHippVols.mat');
leftRestSimilarity={};

for subid=1:24
    subid
    for sid=1:1:2
        restFile=[rootDir,'\Session',dec2base(sid,10,2),'\sub',dec2base(subid,10,2),'Rest\bcWGSdswranrest.nii'];
        restScans=spm_read_vols(spm_vol(restFile));
        restSimi=[];
        for mid=1:1:15
            restS=[];
            for resti=1:1:length(restScans)
                restScan=restScans(:,:,:,resti);
                restScan=reshape(restScan,[1,53*63*52]);
                restScan=restScan(xx);
                SimiAll=[];
                for tid=HippsCoor{subid,sid}{1,mid}'
                    R=corrcoef(allHippVols{subid,sid}{1,mid}(tid,:)',restScan');
                    
                    P=corrcoef([restScan',leftPermHipp{subid,sid}{1,tid}']);
                    P=P(1,2:end);P(isnan(P))=[];
                    nn=find(P<R(2));
                    SimiAll=[SimiAll;numel(nn)/numel(P)];
                end
                restS=[restS;mean(SimiAll)];
            end
            restSimi=[restSimi,restS];
        end
        leftRestSimilarity{subid,sid}=restSimi;
    end
end
save('leftRestSimilarityLof.mat','leftRestSimilarity');




%% Permutation test for resting-state scans similar to each figure item activtion map
clear all;
clc;
rootDir='D:\dataN\';
template=spm_read_vols(spm_vol('RightTemplate.nii'));template=reshape(template,[1,53*63*52]);
xx=find(template>0);
load('RightHippVols.mat');
load('rightLofCleared.mat');
load('rightPermHippVols.mat');


rightRestSimilarity={};

for subid=1:24
    subid
    for sid=1:1:2
        restFile=[rootDir,'\Session',dec2base(sid,10,2),'\sub',dec2base(subid,10,2),'Rest\bcWGSdswranrest.nii'];
        restScans=spm_read_vols(spm_vol(restFile));
        restSimi=[];
        for mid=1:1:15
            restS=[];
            for resti=1:1:length(restScans)
                restScan=restScans(:,:,:,resti);
                restScan=reshape(restScan,[1,53*63*52]);
                restScan=restScan(xx);
                SimiAll=[];
                for tid=HippsCoor{subid,sid}{1,mid}'
                    R=corrcoef(allHippVols{subid,sid}{1,mid}(tid,:)',restScan');
                    
                    P=corrcoef([restScan',rightPermHipp{subid,sid}{1,tid}']);
                    P=P(1,2:end);P(isnan(P))=[];
                    nn=find(P<R(2));
                    SimiAll=[SimiAll;numel(nn)/numel(P)];
                end
                restS=[restS;mean(SimiAll)];
            end
            restSimi=[restSimi,restS];
        end
        rightRestSimilarity{subid,sid}=restSimi;
    end
end
save('rightRestSimilarityLof.mat','rightRestSimilarity');