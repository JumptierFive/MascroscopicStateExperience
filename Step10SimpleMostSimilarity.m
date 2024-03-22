%% location analysis for memory replay
clear all;
clc;
load('subTRMEM.mat');
subAll=[];
for subid=1:1:24
    meanS=zeros(1,552);
    for taskid=1:1:15
        meanS=meanS+max(subMEM{subid,1}{1,taskid});
    end
    meanS=meanS./15;
    subAll=[subAll;meanS];
end