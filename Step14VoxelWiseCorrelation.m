%% network contructions for 160 ROIs and two hippocampus activations memory similarities
clear all;
clc;
rootDir='D:\dataN\';
load('leftRestSimilarityLof.mat');
load('rightRestSimilarityLof.mat');
load('ROI160.mat');


Network160={};
allRL=[];

allRR=[];
nm=0;
parfor subid = 1:1:24
    subid
    ctRL=zeros(1,53*63*52);
    ctRR=zeros(1,53*63*52);
    for sid=1:1:2
        restFile=[rootDir,'\Session',dec2base(sid,10,2),'\sub',dec2base(subid,10,2),'Rest\bcWGSdswranrest.nii'];
        restScans=spm_read_vols(spm_vol(restFile));
        restScans=reshape(restScans,[53*63*52,length(restScans)]);
       
        
        leftHipp=mean(leftRestSimilarity{subid,sid}')';
        rightHipp=mean(rightRestSimilarity{subid,sid}')';
        allTs=[leftHipp,rightHipp,restScans'];
        ctRL=ctRL+corr(allTs(:,1),allTs(:,3:end));
        ctRR=ctRR+corr(allTs(:,2),allTs(:,3:end));
    end
    allRL=[allRL;ctRL];
    allRR=[allRR;ctRR];
    nm=nm+1;
end
[~,~,~,DD1]=ttest(allRL);
[~,~,~,DD2]=ttest(allRR);
mask=spm_read_vols(spm_vol('rAAL116_1mm.nii'));

mask(mask>0)=1;
str=spm_vol('leftTemplate.nii');
str.fname='leftCorres2MR.nii';
spm_write_vol(str,reshape(DD1.tstat,[53,63,52]).*mask);
str.fname='rightCorres2MR.nii';
spm_write_vol(str,reshape(DD2.tstat,[53,63,52]).*mask);


