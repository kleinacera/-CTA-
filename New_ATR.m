function [dailyVor] = New_ATR(O,H,L,C)
% ����ʵ����
% ����
%     O   ÿ�տ��̼�
%     H   ÿ����߼�
%     L   ÿ����ͼ�
%     C   ÿ�����̼�
%     N   ��ʵ�����ʻزⳤ��
u = [];d =[];c = [];TR = [];
for i = 1:length(O)
    u(i) = log(H(i)/O(i));
    d(i) = log(L(i)/O(i));
    c(i) = log(C(i)/O(i));
    dailyVor(i) = 0.511*(u(i)-d(i))^2-0.019*(c(i)*(u(i)+d(i))-2*u(i)*d(i))-0.383*c(i)^2; % ������ʵ������
  %  TR(i) = dailyVor(i)*O(i)*100; % ������ʵ��������
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

