function [PortfolioFinalNav,Date] = MixMomentum_Break_WeightControl(Begin,End,Account,i);
% ������Ʒ�ڲ���Ȩ�ط���
% ����
%     Begin    �ز⿪ʼʱ��,str��ʽ,��'2016-01-01'
%     End      �ز����ʱ��,str��ʽ,��'2017-12-31'
%     Account  �˻���ʼ���

%% �����ʽ�
load IndustrialComdty
load EnergyComdty
load MetalComdty
load BondComdty

if i == 1
    Comdty = IndustrialComdty;
elseif i == 2
    Comdty = EnergyComdty;
elseif i == 3
    Comdty = MetalComdty;
else
    Comdty = BondComdty;
end
    
for j = 1:length(Comdty)
    [PortfolioNav,Date] = Momentum_Break_WeightControl(Comdty(j),Begin,End,Account/length(Comdty));
    PortfolioNavList(:,j) = PortfolioNav;
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

