function [answer] = cap(x,a,b)
% CAP ȡ��ӽ���ֵ������a<b
% ����
%     x      ���������
%     a      �½�
%     b      �Ͻ�
% ���
%     answer �������ı��� 
answer = min(max(x,a),b);
end

