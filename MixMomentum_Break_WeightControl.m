function [PortfolioFinalNav,Date] = MixMomentum_Break_WeightControl(Begin,End,Account,i);
% 大类商品内部等权重分配
% 输入
%     Begin    回测开始时间,str格式,如'2016-01-01'
%     End      回测结束时间,str格式,如'2017-12-31'
%     Account  账户起始金额

%% 分配资金
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

%% 计算指标
% 总收益
CumulativeReturn = PortfolioFinalNav(end)/PortfolioFinalNav(1)-1;
% 年化收益
AnnualReturn = (1+CumulativeReturn)^(245/length(PortfolioFinalNav))-1;
% 波动率
dret = price2ret(PortfolioFinalNav);
StandardDeviation = std(dret)*sqrt(245);
% 夏普比率
SharpRatio = AnnualReturn/StandardDeviation;
% 回撤
DrawDown = [];
for i = 1:length(PortfolioFinalNav)
    DrawDown(i) = PortfolioFinalNav(i)/max(max(PortfolioFinalNav(1:i)))-1;
end
% 最大回撤
MaxDrawDown = min(min(DrawDown));

%% 画图
subplot(2,1,1)
plot(Date,PortfolioFinalNav)
dateaxis('x',2)
subplot(2,1,2)
plot(Date,DrawDown)
dateaxis('x',2)

%% 输出结果
disp(['总收益 = ',num2str(CumulativeReturn)])
disp(['年化收益 = ',num2str(AnnualReturn)])
disp(['波动率 = ',num2str(StandardDeviation)])
disp(['夏普比率 = ',num2str(SharpRatio)])
disp(['最大回撤 = ',num2str(MaxDrawDown)])
end

