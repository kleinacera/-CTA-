function [ATR] = AverageTrueRange(High,Low,Close,N)
% 计算真实波动 

TR = [];
for i = 1:length(High)
    if i == 1
        TR(i) = max([High(i)-Low(i)]);
    else
        TR(i) = max([High(i)-Low(i),High(i)-Close(i-1),Close(i-1)-Low(i)]);
    end
end

[~,ATR] = movavg(TR,5,N,'e');
end

