function [PortfolioNav,Date] = Momentum_Break_WeightControl(Comdty,Begin,End,Account)
% Momentum 基于趋势的CTA策略
% 输入
%     Comdty   期货代码,str格式,如'RB.SHF'
%     Begin    回测开始时间,str格式,如'2016-01-01'
%     End      回测结束时间,str格式,如'2017-12-31'
%     coef     该期货合约每一手的数量，如螺纹钢1手10吨，coef=10
%     Account  账户起始金额

%% 从wind和本地提取数据
load IndustrialComdty
load EnergyComdty
load MetalComdty
load BondComdty
load FutureBaseData
load ChinaGC
load ChinaFutureSettlePrice
load ChinaFutureOpenPrice
load ChinaFutureHighPrice
load ChinaFutureLowPrice
load ChinaFutureClosePrice

w = windmatlab;
[~,~,~,Date] = w.tdays(Begin,End,'TradingCalendar=DCE');
TradingBegin = find(ismember(ChinaFutureOpenPrice.dates,Date(1)));
TradingEnd = find(ismember(ChinaFutureOpenPrice.dates,Date(end)));
AnalystBegin = TradingBegin-20;
ComdtyPosition = find(ismember(ChinaFutureOpenPrice.header,Comdty));
Settle = ChinaFutureOpenPrice.data(TradingBegin:TradingEnd,ComdtyPosition);
High = ChinaFutureHighPrice.data(AnalystBegin:TradingEnd,ComdtyPosition);
High_New = ChinaFutureHighPrice.data(TradingBegin:TradingEnd,ComdtyPosition);
Low = ChinaFutureLowPrice.data(AnalystBegin:TradingEnd,ComdtyPosition);
Low_New = ChinaFutureLowPrice.data(TradingBegin:TradingEnd,ComdtyPosition);
Close = ChinaFutureClosePrice.data(AnalystBegin:TradingEnd,ComdtyPosition);
coef = FutureBaseData.mul(find(ismember(FutureBaseData.code,Comdty)));

%% 计算均线
N = 5;M = 20; % N 端期均线、真实波动率和最大突破值的回测长度; M 长期均线的时间长度
[S,L]=movavg(Close,N,M,'e');
Short = S(21:end);
Long = L(21:end);
[TR] = AverageTrueRange(High,Low,Close,N); % 求出真实日内波动
ATR = TR(21:end);

%% 交易信号识别
% % 买入为1,2,3....，卖空为-1,-2,-3,...，空仓为0 其中的数字为仓位数量
Weight = zeros(length(Settle),1); %手数
for i = M:length(Settle)-1
    if Weight(i-1) == 0
        if Short(i)>Long(i) && High_New(i)==max(High_New(i-M+1:i))
            Weight(i) = 1;  
        elseif Short(i)<Long(i) && Low_New(i)==min(Low_New(i-M+1:i))
            Weight(i) = -1;
        end
    else
        if Weight(i-1) >= 1
            if Weight(i-1)>Weight(i-2)
                BuyStand = Settle(i); %购入价格
            end
            if Short(i)>=Long(i)
                if High_New(i)>=BuyStand+0.5*ATR(i) %持续涨
                    Mul = floor((High_New(i)-BuyStand)/(0.5*ATR(i)));
                    Weight(i) = Weight(i-1)+Mul;
                elseif Low_New(i)<=BuyStand-2*ATR(i)||Low_New(i)==min(Low_New(i-M+1:i))
                     Weight(i) = 0;
                else
                    Weight(i) = Weight(i-1);
                end
            else
                Weight(i) = 0;
            end
        else
            if Weight(i-1)<Weight(i-2)
                SellStand = Settle(i); %购入价格
            end
            if Short(i)<=Long(i)
                if Low_New(i)<=SellStand-0.5*ATR(i)
                    Mul = floor((SellStand-Low_New(i))/(0.5*ATR(i)));
                    Weight(i) = Weight(i-1)-Mul;
                elseif High_New(i)>=SellStand+2*ATR(i-1)||High_New(i)==max(High_New(i-M+1:i))
                    Weight(i) = 0;
                else
                    Weight(i) = Weight(i-1);
                end
            else
                Weight(i) = 0;
            end
        end
    end
end

%% 根据权重计算净值
PortfolioNav(1:2,1) = Account; %净值
PortfolioWeight = zeros(length(Weight),1); %手数
for i = 3:length(Weight)
    if Weight(i-1)==0
        if Weight(i-2)==0 %空仓
            PortfolioNav(i,1) = PortfolioNav(i-1,1)*(1+(ChinaGC.data(find(ismember(ChinaGC.dates,Date(i))),1)/365));
        else %平仓
            Amount = PortfolioWeight(i-2)*coef*Settle(i);
            [FinalCost] = Cost(Comdty,PortfolioWeight(i-2),Amount);
            PortfolioNav(i,1) = PortfolioNav(i-1,1)+PortfolioWeight(i-2)*coef*(Settle(i)-Settle(i-1))+(PortfolioNav(i-1,1)-Amount/5)*(ChinaGC.data(find(ismember(ChinaGC.dates,Date(i))),1)/365)-FinalCost;
        end
    else %开仓
        Unit = round((0.01*PortfolioNav(i-1))/(ATR(i-1)*coef)); % 求出单位加仓量
        PortfolioWeight(i-1) = cap(PortfolioWeight(i-2)+(Weight(i-1)-Weight(i-2))*Unit,-1*round(PortfolioNav(i-1)/(coef*Settle(i))),1*round(PortfolioNav(i-1)/(coef*Settle(i)))); % 限制最大仓位
        Amount = PortfolioWeight(i-2)*coef*Settle(i);
        [FinalCost] = Cost(Comdty,PortfolioWeight(i-2),Amount);
        PortfolioNav(i,1) = PortfolioNav(i-1,1)+PortfolioWeight(i-2)*coef*(Settle(i)-Settle(i-1))+(PortfolioNav(i-1,1)-Amount/5)*(ChinaGC.data(find(ismember(ChinaGC.dates,Date(i))),1)/365)-FinalCost;
    end
end

%% 计算指标
% 总收益
CumulativeReturn = PortfolioNav(end)/PortfolioNav(1)-1;
% 年化收益
AnnualReturn = (1+CumulativeReturn)^(250/length(PortfolioNav))-1;
% 波动率
dret = price2ret(PortfolioNav);
StandardDeviation = std(dret)*sqrt(250);
% 夏普比率
SharpRatio = AnnualReturn/StandardDeviation;
% 回撤
DrawDown = [];
for i = 1:length(PortfolioNav)
    DrawDown(i) = PortfolioNav(i)/max(max(PortfolioNav(1:i)))-1;
end
% 最大回撤
MaxDrawDown = min(min(DrawDown));

%% 画图
subplot(3,1,1)
plot(Date,PortfolioNav)
dateaxis('x',2)
subplot(3,1,2)
plot(Date,DrawDown)
dateaxis('x',2)
subplot(3,1,3)
bar(Date,PortfolioWeight)
dateaxis('x',2)

%% 输出结果
disp(['总收益 = ',num2str(CumulativeReturn)])
disp(['年化收益 = ',num2str(AnnualReturn)])
disp(['波动率 = ',num2str(StandardDeviation)])
disp(['夏普比率 = ',num2str(SharpRatio)])
disp(['最大回撤 = ',num2str(MaxDrawDown)])
end