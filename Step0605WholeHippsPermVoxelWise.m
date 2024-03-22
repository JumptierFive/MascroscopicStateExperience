%% Permutation test for resting-state scans similar to each figure item activtion map
clear all;
clc;
rootDir='D:\dataN\';
template=spm_read_vols(spm_vol('LeftTemplate.nii'));template=reshape(template,[1,53*63*52]);
x1=find(template>0);

template=spm_read_vols(spm_vol('RightTemplate.nii'));template=reshape(template,[1,53*63*52]);
x2=find(template>0);

leftPermHipp={};
rightPermHipp={};
for subid=1:24
    subid
    for sid=1:1:2
        for tid=1:1:8
            perFile=[rootDir,'\ReplayFMRI\FirstLevelPerm\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2),'\task_',dec2base(tid,10,2)];
            permFiles=filename_list(perFile,'*.nii');
            permTemp=[];
            permTemp1=[];
            for pi=1:1:100
                tempTN=reshape(spm_read_vols(spm_vol(permFiles{pi})),[1,53*63*52]);
                permTemp=[permTemp;tempTN(x1)];
                permTemp1=[permTemp1;tempTN(x2)];
            end
            leftPermHipp{subid,sid}{1,tid}=permTemp;
            rightPermHipp{subid,sid}{1,tid}=permTemp1;
        end
    end
end

save('leftPermHippVols.mat','leftPermHipp');
save('rightPermHippVols.mat','rightPermHipp');




