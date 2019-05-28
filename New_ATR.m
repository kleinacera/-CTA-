function [dailyVor] = New_ATR(O,H,L,C)
% 求真实波动
% 输入
%     O   每日开盘价
%     H   每日最高价
%     L   每日最低价
%     C   每日收盘价
%     N   真实波动率回测长度
u = [];d =[];c = [];TR = [];
for i = 1:length(O)
    u(i) = log(H(i)/O(i));
    d(i) = log(L(i)/O(i));
    c(i) = log(C(i)/O(i));
    dailyVor(i) = 0.511*(u(i)-d(i))^2-0.019*(c(i)*(u(i)+d(i))-2*u(i)*d(i))-0.383*c(i)^2; % 日内真实波动率
  %  TR(i) = dailyVor(i)*O(i)*100; % 日内真实波动幅度
end
% 
% ATR = [];
% for i = 1:length(TR)
%     if i == 1
%         ATR(i) = TR(i);
%     else
%         ATR(i) = ((N-1)*ATR(i-1)+TR(i))/N;
%     end
% end

