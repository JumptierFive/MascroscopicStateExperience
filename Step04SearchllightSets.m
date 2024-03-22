clear all;
clc;
rootDir='D:\Work\dataN\ReplayFMRI\FirstLevel';
% allFiles=filename_list(rootDir,'sub*');
% template=ones(53,63,52);
% mask=spm_read_vols(spm_vol('HippMask.nii'));
% for i=1:1:length(allFiles)
%     i
%     for si=1:1:2
%         for taskid=1:1:8
%             for mid=1:1:15
%                 templateH=spm_vol([allFiles{i},'\Session',dec2base(si,10,2),'\task_',dec2base(taskid,10,2),'\taskid_',dec2base(mid,10,3),'.nii']);
%                 templateN=spm_read_vols(templateH);
%                 templateN(templateN~=0)=1;
%                 template=template.*templateN;
%             end
%         end
%     end
% end
% template=template.*mask;
% templateH.fname='template.nii';
% spm_write_vol(templateH,template);

template=spm_read_vols(spm_vol('template.nii'));

allSL={};
for li=1:1:53*63*52
    if template(li)==0
        allSL{li}=[];
        continue;
    end
    tempRadius=[];
    
    [ki,kki,kkki] = ind2sub([53,63,52],li);
    for mi=1:1:53
        for mmi=1:1:63
            for mmmi=1:1:52
                if template(mi,mmi,mmmi)==0
                    continue;
                end
                if (ki-mi)^2+(kki-mmi)^2+(kkki-mmmi)^2<=3^2
                    tempRadius=[tempRadius;[mi,mmi,mmmi]];
                end
            end
        end
    end
    allSL{li}=sub2ind([53,63,52],tempRadius(:,1),tempRadius(:,2),tempRadius(:,3));
end

save('allSL.mat','allSL');