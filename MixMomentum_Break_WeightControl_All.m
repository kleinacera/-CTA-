function MixMomentum_Break_WeightControl_All(Begin,End,Account)
% ���д�����Ʒ��Ȩ�ط���
% ����
%     Begin    �ز⿪ʼʱ��,str��ʽ,��'2016-01-01'
%     End      �ز����ʱ��,str��ʽ,��'2017-12-31'
%     Account  �˻���ʼ���

%% �����ʽ�
PortfolioNavList = [];
SingleAccount = Account/4;
for i = 1:4
    [PortfolioNav,Date] = MixMomentum_Break_WeightControl(Begin,End,SingleAccount,i);
    PortfolioNavList(:,i) = PortfolioNav;
end
PortfolioFinalNav = sum(PortfolioNavList,2);

%% ����ָ��
% ������
CumulativeReturn = PortfolioFinalNav(end)/PortfolioFinalNav(1)-1;
% �껯����
AnnualReturn = (1+CumulativeReturn)^(245/length(PortfolioFinalNav))-1;
% ������
dret = price2ret(PortfolioFinalNav);
StandardDeviation = std(dret)*sqrt(245);
% ���ձ���
SharpRatio = AnnualReturn/StandardDeviation;
% �س�
DrawDown = [];
for i = 1:length(PortfolioFinalNav)
    DrawDown(i) = PortfolioFinalNav(i)/max(max(PortfolioFinalNav(1:i)))-1;
end
% ���س�
MaxDrawDown = min(min(DrawDown));

%% ��ͼ
subplot(2,1,1)
plot(Date,PortfolioFinalNav)
dateaxis('x',2)
subplot(2,1,2)
plot(Date,DrawDown)
dateaxis('x',2)

%% ������
disp(['������ = ',num2str(CumulativeReturn)])
disp(['�껯���� = ',num2str(AnnualReturn)])
disp(['������ = ',num2str(StandardDeviation)])
disp(['���ձ��� = ',num2str(SharpRatio)])
disp(['���س� = ',num2str(MaxDrawDown)])
end
