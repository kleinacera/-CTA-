function [answer] = cap(x,a,b)
% CAP 取最接近的值，其中a<b
% 输入
%     x      被处理变量
%     a      下界
%     b      上界
% 输出
%     answer 被处理后的变量 
answer = min(max(x,a),b);
end

