%% ROIs calculation 
% clear all;
% clc;
% rootDir='D:\dataN\';
% template=spm_read_vols(spm_vol('rDosenbach160_3mm.nii'));
% lbCoor={};
% for li=1:1:160
%     xx=find(template==li);
%     lbCoor{li}=xx;
% end
% 
% ROI160={};
% for subid=1:1:24
%     subid
%     for sid=1:1:2
%         restFile=[rootDir,'\Session',dec2base(sid,10,2),'\sub',dec2base(subid,10,2),'Rest\bcWGSdswranrest.nii'];
%         restScans=spm_read_vols(spm_vol(restFile));
%         tempROI=[];
%         for resti=1:1:length(restScans)
%             roiN=[];
%             funcA=restScans(:,:,:,resti);
%             for roi=1:1:160
%                 roiN=[roiN,mean(funcA(lbCoor{roi}))];
%             end
%             tempROI=[tempROI;roiN];
%         end
%         ROI160{subid,sid}=tempROI;
%     end
% end
% save('ROI160.mat','ROI160');



clear all;
clc;
rootDir='D:\dataN\';
template=spm_read_vols(spm_vol('rRandom625_WithAALBoundary_3mm.nii'));
template=reshape(template,[1,53*63*52]);
lbCoor={};
for li=1:1:625
    xx=find(template==li);
    lbCoor{li}=xx;
end


for subid=1:1:24
    subid
    for sid=1:1:2
        restFile=[rootDir,'\Session',dec2base(sid,10,2),'\sub',dec2base(subid,10,2),'Rest\bcWGSdswranrest.nii'];
        restScans=spm_read_vols(spm_vol(restFile));
        allScans=[];
        for scn=1:1:length(restScans)
            ctScans=reshape(restScans(:,:,:,scn),[1,53*63*52]);
            allScans=[allScans;ctScans];
        end
        
        allROI={};
        for tid=1:1:625
            allROI{tid}=allScans(:,lbCoor{tid});
        end
        save(['ROI625/subid',dec2base(subid,10,2),'Session',dec2base(sid,10,2),'.mat'],'allROI');
    end
end