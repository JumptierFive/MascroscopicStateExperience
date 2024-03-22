%% find timepoints where most Replay happens
clear all;
clc;
rootDir='D:\dataN\ReplayFMRI\RestingStateSimilarityPermTest';
ReplayTR={};
for subid=1:1:24
    for sid=1:1:2
        memTR=[];
        allFiles=filename_list([rootDir,'\subid',dec2base(subid,10,2),'\Session',dec2base(sid,10,2)],'rest*.mat');
        for rid=1:1:length(allFiles)
            load([allFiles{rid}]);
            xxn=max(max(max(allMEM)));
            memTR=[memTR;xxn/8];
        end
        ReplayTR{subid}{sid}=memTR;
    end
end
