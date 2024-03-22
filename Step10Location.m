%% location analysis for memory replay
clear all;
clc;
% rootDir='D:\dataN\ReplayFMRI\RestingStateSimilarityPermTest';
% subMEM={};
% 
% %% rearrange the similarity matrix
% for subid=1:1:24
%     subid
%     for sid=1:1:2
%         allFiles=filename_list([rootDir,'\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2)],'rest*.mat');
% 
%         for mid=1:1:15
%             midN=[];
%             for rid=1:1:length(allFiles)
%                 load([allFiles{rid}]);
%                 midN=[midN;allMEM(mid,:)./8];
%             end
%             subMEM{subid,sid}{mid}=midN;
%         end
%     end
% end
% save('subTRMEM.mat','subMEM');
load('subTRMEM.mat');
%% similarity between whole brain
allSimilarity=cell(24,2);

atlas=spm_read_vols(spm_vol('rDosenbach160_3mm.nii'));
atlas=reshape(atlas,[1,53*63*52]);
atlasN={};
for i=1:1:160
    atlasN{i}=find(atlas==i);
end
rootDir='D:\dataN';
for subid=1:1:24
    for sid=1:1:2
        %% calcualting ROI means during resting state 
        restFile=[rootDir,'\Session',dec2base(sid,10,2),'\sub',dec2base(subid,10,2),'Rest\bcWGSdswranrest.nii'];
        restFF=spm_vol(restFile);
        restScans=spm_read_vols(restFF);
        allROIs=[];
        for rid=1:1:length(restFF)
            ROI=[];
            tempAt=reshape(restScans(:,:,:,rid),[53,63,52]);
            for roi=1:1:160
                ROI=[ROI,sum(tempAt(atlasN{roi}))/numel(atlasN{roi})];
            end
            allROIs=[allROIs;ROI];
        end
        
        %% calculating rois in each activaiton map for memory items
        allTROIs={};
        for tid=1:1:8
            itROI=[];
            for iti=1:1:15
                
                fileName=[rootDir,'\ReplayFMRI\FirstLevel\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2),'\taskid_',dec2base(tid,10,3),'.nii'];
                templ=reshape(spm_read_vols(spm_vol(fileName)),[1,53*63*52]);
                ROI=[];
                for roi=1:1:160
                    ROI=[ROI,sum(templ(atlasN{roi}))/numel(atlasN{roi})];
                end
                itROI=[itROI;ROI];
            end
            allTRIs{tid}=itROI;
        end
        
        
        
        
        
    end
end