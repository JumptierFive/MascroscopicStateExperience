%% compare activation similarity between real and dummy sitations
clear all;
clc;
% rootDir='D:\dataN\ReplayFMRI\FirstLevel';
% similarityHipp={};
% mask=spm_read_vols(spm_vol('RightTemplate.nii'));
% mask=reshape(mask,[1,53*63*52]);
% xyz=find(mask>0);
% 
% allHippVols={};
% for subid=1:1:24
%     subid
%     for sid=1:1:2
%         for mid=1:1:15
%             hippVols=[];
%             
%             for tid=1:1:8
%                 tempT=spm_read_vols(spm_vol([rootDir,'\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2),'\taskid_',dec2base(mid,10,3),'.nii']));
%                 tempT=reshape(tempT,[1,53*63*52]);
%                 tempT=tempT(xyz);
%                 hippVols=[hippVols;tempT];
%             end
%             allHippVols{subid,sid}{1,mid}=hippVols;
%         end
%     end
% end
% save('RightHippVols.mat','allHippVols');

%% test if the similar stimuli aroused more similarity activation map in left/right hippcampus than that in between different stimulis





