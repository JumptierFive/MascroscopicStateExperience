%% Permutation test for resting-state scans similar to each figure item activtion map
clear all;
clc;
rootDir='D:\dataN\';
template=spm_read_vols(spm_vol('LeftTemplate.nii'));template=reshape(template,[1,53*63*52]);
xx=find(template>0);
load('leftLofCleared.mat');


template=spm_read_vols(spm_vol('rRandom625_WithAALBoundary_3mm.nii'));
template=reshape(template,[1,53*63*52]);
lbCoor={};
for li=1:1:625
    xx=find(template==li);
    lbCoor{li}=xx;
end



leftRestSimilarity={};

for subid=1:24
    subid
    for sid=1:1:2
        load(['ROI625/subid',dec2base(subid,10,2),'Session',dec2base(sid,10,2),'.mat']);
        
        restFile=[rootDir,'\Session',dec2base(sid,10,2),'\sub',dec2base(subid,10,2),'Rest\bcWGSdswranrest.nii'];
        restScans=spm_read_vols(spm_vol(restFile));
        restSimi=[];
        
        allPerm={};
        for tid=1:1:8
            perDir=filename_list([rootDir,'ReplayFMRI\FirstLevelPerm\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2)],'*nii');
            ctPerm=[];
            for zid=1:1:100
                tpPerm=spm_read_vols(spm_vol(perDir{zid}));
                tpPerm=reshape(tpPerm,[1,53*63*52]);
                ctPerm=[ctPerm;tpPerm];
            end
            allPerm{tid}=ctPerm;
        end
        
        
        
        midAll={};
        parfor mid=1:15
            mid
            allAct=[];

            for tid=1:1:8
                fileName=[rootDir,'ReplayFMRI\FirstLevel\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2),'\taskid_',dec2base(tid,10,3),'.nii'];
                tempT=spm_read_vols(spm_vol(fileName));
                tempT=reshape(tempT,[1,53*63*52]);
                allAct=[allAct;tempT];
            end
            
            
            roiAll={};
            orderN=HippsCoor{subid,sid}{1,mid}';
            for roi=1:625
                roi
                ctACT=allAct(:,lbCoor{roi});
                restS=[];
                coorCt=lbCoor{roi};
                for resti=1:length(restScans)
                    restScan=restScans(:,:,:,resti);
                    restScan=reshape(restScan,[1,53*63*52]);
                    restScan=restScan(coorCt);
                    SimiAll=[];


                    for tid=orderN
                        P=corrcoef([restScan',ctACT(tid,:)',allPerm{tid}(:,coorCt)']);
                        R=P(1,2);
                        if isnan(R)
                            continue;
                        end
                        P=P(1,3:end);P(isnan(P))=[];
                        nn=find(P<R);
                        SimiAll=[SimiAll;numel(nn)/numel(P)];
                    end
                    restS=[restS;mean(SimiAll)];
                end
                roiAll{roi}=restS;
            end
            midAll{mid}=roiAll;
        end
        save(['ROI625Similarity\subid',dec2base(subid,10,2),'_Sessio',dec2base(sid,10,2),'.mat'],'midAll');
    end
    
    
end