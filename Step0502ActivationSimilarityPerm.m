%% compare activation similarity between real and dummy sitations
clear all;
clc;
rootDir='D:\dataN\ReplayFMRI\FirstLevelPerm';
similarityHipp={};
mask=spm_read_vols(spm_vol('LeftTemplate.nii'));
mask=reshape(mask,[1,53*63*52]);
xyz=find(mask>0);

allHippVols={};
for subid=1:1:24
    subid
    for sid=1:1:2
        allVols={};
        for tid=1:1:8
            hippVols=[];
            for kid=1:1:100
                tempT=spm_read_vols(spm_vol([rootDir,'\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2),'\taskid_',dec2base(kid,10,3),'.nii']));
                tempT=reshape(tempT,[1,53*63*52]);
                tempT=tempT(xyz);
                hippVols=[hippVols;tempT];
            end
            allVols{tid}=hippVols;
        end
        
        similarityMat=cell(8,8);
        for xi=1:1:8
            xi
            for yi=1:1:8
                if xi==yi
                    continue;
                end
                allCor=[];
                for pid=1:1:100
                    for ppid=1:1:100
                        
                        R=corrcoef(allVols{xi}(pid,:)',allVols{yi}(ppid,:)');
                        allCor=[allCor;R(2)];
                    end
                end
                similarityMat{xi,yi}=allCor;similarityMat{yi,xi}=allCor;
            end
        end
        allHippVols{subid,sid}=similarityMat;
    end
end
save('LeftHippVolsPerm.mat','allHippVols');


clear all;
clc;
rootDir='D:\dataN\ReplayFMRI\FirstLevelPerm';
similarityHipp={};
mask=spm_read_vols(spm_vol('RightTemplate.nii'));
mask=reshape(mask,[1,53*63*52]);
xyz=find(mask>0);

allHippVols={};
for subid=1:1:24
    subid
    for sid=1:1:2
        allVols={};
        for tid=1:1:8
            hippVols=[];
            for kid=1:1:100
                tempT=spm_read_vols(spm_vol([rootDir,'\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2),'\taskid_',dec2base(kid,10,3),'.nii']));
                tempT=reshape(tempT,[1,53*63*52]);
                tempT=tempT(xyz);
                hippVols=[hippVols;tempT];
            end
            allVols{tid}=hippVols;
        end
        
        similarityMat=cell(8,8);
        for xi=1:1:8
            xi
            for yi=1:1:8
                if xi==yi
                    continue;
                end
                allCor=[];
                for pid=1:1:100
                    for ppid=1:1:100
                        
                        R=corrcoef(allVols{xi}(pid,:)',allVols{yi}(ppid,:)');
                        allCor=[allCor;R(2)];
                    end
                end
                similarityMat{xi,yi}=allCor;similarityMat{yi,xi}=allCor;
            end
        end
        allHippVols{subid,sid}=similarityMat;
    end
end
save('RightHippVolsPerm.mat','allHippVols');

