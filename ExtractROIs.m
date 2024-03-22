%% extract 12 ROIs
% Anterior medial PFC	1	55	26
% Ventro-medial PFC	−3	40	0
% Left SFG	−14	36	59
% Right SFG	17	35	58
% Left ITG	−62	−33	−20
% Right ITG	66	−17	−19
% Left parahippocampal gyrus	−22	−26	−21
% Right parahippocampal gyrus	25	−26	−18
% PCC	−2	−29	39
% Left lateral parietal	−47	−71	35
% Right lateral parietal	54	−61	36
% Posterior cingulate	3	−53	6

clear all;
clc;
rootDir='D:\dataN\';
% make seven ROIs template
niHead=spm_vol([rootDir,'Session01\sub01Rest\bcWGSdswranrest.nii']);
image=spm_read_vols(niHead(1,:));
image(:,:,:)=0;

allInfo=readtable('264ROIs.xls');
coorROI=[allInfo.X,allInfo.Y,allInfo.Z];
coorROI={};
for i=1:1:8
    xx=find(allInfo.SubRegion==i);
    coorROI{i}=[allInfo.X(xx),allInfo.Y(xx),allInfo.Z(xx)];
end


%coorROI=[1 55 26;-3 40 0;-14 36 59;17 35 58; -62 -33 -20;66 -17 -19; -22 -26 -21; 25 -26 -18; -2 -29 39;-47 -71 35;64 -61 36;3 -53 6];
allXYZcoors={};
tempImage=image;
for i=1:1:length(coorROI)
    i
    % coorXYZ=mni2xyz(niHead(1),coorROI(i,:));
    for mi=1:1:length(coorROI{i}(:,1))
        mi
        coorXYZ=mni2cor(coorROI{i}(mi,:),niHead(1).mat);
        coorXYZroi=ROIcoors(niHead(1),coorXYZ,6/3);
        for ki=1:1:length(coorXYZroi(:,1))
            tempImage(coorXYZroi(ki,1),coorXYZroi(ki,2),coorXYZroi(ki,3))=i;
        end
    end
end

header=niHead(1);
header.fname=['wholeBrain_ROI.nii'];
spm_write_vol(header,tempImage);

allXYZcoors={};
for i=1:1:8
    tempIm=tempImage;tempIm(tempIm~=i)=0;tempIm(tempIm==i)=1;
    allXYZcoors{i}=tempIm;
end



outputDir=['D:\dataN\ReplayFMRI\wholeTimeSeriesS01\'];mkdir(outputDir);
allFiles=filename_list([rootDir,'Session01'],'sub*rest');
for i=1:1:length(allFiles)
    disp(['subid_',dec2base(i,10,2),'  start']);
    subTS=[];
    strN=spm_vol([allFiles{i},'\bcWGSdswranrest.nii']);
    for sti=1:1:length(strN)
        vols=spm_read_vols(strN(sti));
        tempN=[];
        for bgi=1:1:length(allXYZcoors)
            sumM=sum(sum(sum(vols.*allXYZcoors{bgi})))/sum(sum(sum(allXYZcoors{bgi})));
            tempN=[tempN,sumM];
        end
        subTS=[subTS;tempN];
    end
    save([outputDir,'sub_',dec2base(i,10,2),'DMN.mat'],'subTS');
end

outputDir=['D:\dataN\ReplayFMRI\wholeTimeSeriesS02\'];mkdir(outputDir);
allFiles=filename_list([rootDir,'Session02'],'sub*rest');
for i=1:1:length(allFiles)
    disp(['subid_',dec2base(i,10,2),'  start']);
    subTS=[];
    strN=spm_vol([allFiles{i},'\bcWGSdswranrest.nii']);
    for sti=1:1:length(strN)
        vols=spm_read_vols(strN(sti));
        tempN=[];
        for bgi=1:1:length(allXYZcoors)
            sumM=sum(sum(sum(vols.*allXYZcoors{bgi})))/sum(sum(sum(allXYZcoors{bgi})));
            tempN=[tempN,sumM];
        end
        subTS=[subTS;tempN];
    end
    save([outputDir,'sub_',dec2base(i,10,2),'DMN.mat'],'subTS');
end


function [XYZcoOrds] = mni2xyz(vol,MNIcoOrds)

[m n] = size(MNIcoOrds);
if n > m
    MNIcoOrds = MNIcoOrds';
end

XYZcoOrds = round(inv(vol.mat)*[MNIcoOrds; 1]);
end

function [coordinates] = ROIcoors(vol,coorXYZ,distance)
coordinates=[];
for i=1:1:vol.dim(1)
    for k=1:1:vol.dim(2)
        for j=1:1:vol.dim(3)
            if ((i-coorXYZ(1))^2+(k-coorXYZ(2))^2+(j-coorXYZ(3))^2)<distance^2
                coordinates=[coordinates;[i,k,j]];
            end
        end
    end
end
end


function coordinate = mni2cor(mni, T)

if isempty(mni)
    coordinate = [];
    return;
end

if nargin == 1
	T = ...
        [-4     0     0    84;...
         0     4     0  -116;...
         0     0     4   -56;...
         0     0     0     1];
end

coordinate = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]*(inv(T))';
coordinate(:,4) = [];
coordinate = round(coordinate);
end