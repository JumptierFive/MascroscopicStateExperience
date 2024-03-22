%% Permutation test for resting-state scans similar to each figure item activtion map
clear all;
clc;
load('allSL.mat');
headFF=spm_vol('template.nii');

allSLN=allSL;
indN=[];
for i=1:1:length(allSLN)
    if numel(allSLN{i})>=3
        indN=[indN,i];
    end
end

% aa=allSLN{93199};
% template=zeros(53,63,52);
% template(aa)=1;
% % for i=1:1:length(aa)
% %     template(aa(i,1),aa(i,2),aa(i,3))=1;
% % end
% an=spm_vol('template.nii');
% an.fname='templateTT.nii';
% spm_write_vol(an,template);
% aaaa
rootDir='D:\dataN\';

for subid=1:24
    for sid=1:1:2
        restFile=[rootDir,'\Session',dec2base(sid,10,2),'\sub',dec2base(subid,10,2),'Rest\bcWGSdswranrest.nii'];
        restFF=spm_vol(restFile);
        restScans=spm_read_vols(restFF);
%         for ski=1:1:length(restFile)
%             mn=['restScan',dec2base(ski,10,4),'=restScans(:,:,:,ski)'];
%             eval(mn);
%         end
        outputDir=[rootDir,'\ReplayFMRI\RestingStateSimilarityPermTest\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2)];mkdir(outputDir);
        
        for resti=1:1:length(restFF)
            restScan=restScans(:,:,:,resti);
            allMEM=zeros(15,numel(indN));
            parfor tid=1:8
                
                perFile=[rootDir,'\ReplayFMRI\FirstLevelPerm\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2)];
                permFiles=filename_list(perFile,'*.nii');
                permTemp={};
                for pi=1:1:100
                    permTemp{pi}=spm_read_vols(spm_vol(permFiles{pi}));
                end
                
                meTemps={};
                for mid=1:1:15
                    meTemps{mid}=spm_read_vols(spm_vol([rootDir,'\ReplayFMRI\FirstLevel\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2),'\taskid_',dec2base(mid,10,3),'.nii']));
                end
                
                allMemNSI=[];
                for si=1:length(indN)
                    allPermN=[];
                    coorCT=allSLN{indN(si)};
                    for pi=1:1:100
                   
                        R=corrcoef(restScan(coorCT),permTemp{pi}(coorCT));
                        allPermN=[allPermN,R(2)];
                    end
                    allMemN=[];
                    for mid=1:1:15
                        R=corrcoef(restScan(coorCT),meTemps{mid}(coorCT));
                        xxm=find(allPermN<R(2));
                        xxn=numel(xxm);
                        allMemN=[allMemN;xxn];
                    end
                    allMemNSI=[allMemNSI,allMemN];
                end
                allMEM=allMEM+allMemNSI;
            end

            headFF=[outputDir,'\rest',dec2base(resti,10,4),'.mat'];
            save(headFF,'allMEM');
        end
        

        
    end
end
