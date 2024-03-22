%% %% network contructions for 160 ROIs and two hippocampus activations memory similarities
clear all;
clc;
rootDir='D:\dataN\';
load('leftRestSimilarityLof.mat');
load('rightRestSimilarityLof.mat');


allRL=[];

allRR=[];
nm=0;
for subid = 1:1:24
    subid
    ctRL=zeros(1,625);
    for sid=1:1:2
        load(['ROI625Similarity\subid',dec2base(subid,10,2),'_Sessio',dec2base(sid,10,2),'.mat']);
        wbSimi=zeros(length(midAll{1,1}{1,1}),625);
        for miid=1:1:15
            wbSimi=wbSimi+cell2mat(midAll{1,miid});
        end
        
        leftHipp=mean(leftRestSimilarity{subid,sid}')';
        rightHipp=mean(rightRestSimilarity{subid,sid}')';
        Hipp=rightHipp+leftHipp;%Hipp=Hipp(1:end,:);
        ctRL=ctRL+corr(Hipp,wbSimi(1:end,:));
    end
    allRL=[allRL;ctRL];

    nm=nm+1;
end
[~,~,~,DD1]=ttest(allRL);
mask=spm_read_vols(spm_vol('rRandom625_WithAALBoundary_3mm.nii'));

str=spm_vol('leftTemplate.nii');
template=spm_read_vols(str);
template(:,:,:)=0;
template=reshape(template,[1,53*63*52]);
mask=reshape(mask,[1,53*63*52]);
for roi=1:1:625
    template(mask==roi)=DD1.tstat(roi);
end


str.fname='newROI625SimilarityStrikeWithHipp.nii';
spm_write_vol(str,reshape(template,[53,63,52]));