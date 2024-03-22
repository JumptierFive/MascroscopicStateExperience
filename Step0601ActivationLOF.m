%% lof extreme value detection for hippcampus actiation
clear all;
clc
load('leftHippVols.mat');
k=2;

HippsCoor={};
for subid=1:1:24
    for sid=1:1:2
        for mid=1:1:15
            data=allHippVols{subid,sid}{1,mid};
            lof_scores=calculate_lof(data,k);
            standr=mean(lof_scores)+1.5*std(lof_scores);
            xx=find(lof_scores<standr);
            HippsCoor{subid,sid}{1,mid}=xx;
        end
    end
end
save('leftLofCleared.mat','HippsCoor');



load('RightHippVols.mat');
k=3;

HippsCoor={};
for subid=1:1:24
    for sid=1:1:2
        for mid=1:1:15
            data=allHippVols{subid,sid}{1,mid};
            lof_scores=calculate_lof(data,k);
            standr=mean(lof_scores)+1.5*std(lof_scores);
            xx=find(lof_scores<standr);
            HippsCoor{subid,sid}{1,mid}=xx;
        end
    end
end
save('rightLofCleared.mat','HippsCoor');



function lof_scores = calculate_lof(data, k)
    % 输入参数：
    % - data: 输入数据矩阵，每行是一个样本，每列是一个特征
    % - k: 邻近样本的数量

    [num_samples, num_features] = size(data);
    
    % 初始化 LOF 得分
    lof_scores = zeros(num_samples, 1);
    R=corrcoef(data');
    for i = 1:num_samples
        % 计算样本 i 与其他样本的距离
        distances = R(:,i).*-1;
        % 找到 k 个最近邻样本的索引
        [~, idx] = mink(distances, k + 1); 
        % 计算局部密度
        local_density_i = 1 / mean(distances(idx(2:end)));
        
        % 计算相对密度和 LOF
        relative_density_i = zeros(k, 1);
        for j = 1:k
            idx_j = setdiff(idx, idx(j + 1));  % 排除第 j 个邻近样本
            local_density_j = 1 / mean(distances(idx_j));
            relative_density_i(j) = local_density_i / local_density_j;
        end
        
        lof_scores(i) = mean(relative_density_i);
    end
end