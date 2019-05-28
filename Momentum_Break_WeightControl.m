function [PortfolioNav,Date] = Momentum_Break_WeightControl(Comdty,Begin,End,Account)
% Momentum �������Ƶ�CTA����
% ����
%     Comdty   �ڻ�����,str��ʽ,��'RB.SHF'
%     Begin    �ز⿪ʼʱ��,str��ʽ,��'2016-01-01'
%     End      �ز����ʱ��,str��ʽ,��'2017-12-31'
%     coef     ���ڻ���Լÿһ�ֵ������������Ƹ�1��10�֣�coef=10
%     Account  �˻���ʼ���

%% ��wind�ͱ�����ȡ����
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

%% �������
N = 5;M = 20; % N ���ھ��ߡ���ʵ�����ʺ����ͻ��ֵ�Ļزⳤ��; M ���ھ��ߵ�ʱ�䳤��
[S,L]=movavg(Close,N,M,'e');
Short = S(21:end);
Long = L(21:end);
[TR] = AverageTrueRange(High,Low,Close,N); % �����ʵ���ڲ���
ATR = TR(21:end);

%% �����ź�ʶ��
% % ����Ϊ1,2,3....������Ϊ-1,-2,-3,...���ղ�Ϊ0 ���е�����Ϊ��λ����
Weight = zeros(length(Settle),1); %����
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
                BuyStand = Settle(i); %����۸�
            end
            if Short(i)>=Long(i)
                if High_New(i)>=BuyStand+0.5*ATR(i) %������
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
                SellStand = Settle(i); %����۸�
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

%% ����Ȩ�ؼ��㾻ֵ
PortfolioNav(1:2,1) = Account; %��ֵ
PortfolioWeight = zeros(length(Weight),1); %����
for i = 3:length(Weight)
    if Weight(i-1)==0
        if Weight(i-2)==0 %�ղ�
            PortfolioNav(i,1) = PortfolioNav(i-1,1)*(1+(ChinaGC.data(find(ismember(ChinaGC.dates,Date(i))),1)/365));
        else %ƽ��
            Amount = PortfolioWeight(i-2)*coef*Settle(i);
            [FinalCost] = Cost(Comdty,PortfolioWeight(i-2),Amount);
            PortfolioNav(i,1) = PortfolioNav(i-1,1)+PortfolioWeight(i-2)*coef*(Settle(i)-Settle(i-1))+(PortfolioNav(i-1,1)-Amount/5)*(ChinaGC.data(find(ismember(ChinaGC.dates,Date(i))),1)/365)-FinalCost;
        end
    else %����
        Unit = round((0.01*PortfolioNav(i-1))/(ATR(i-1)*coef)); % �����λ�Ӳ���
        PortfolioWeight(i-1) = cap(PortfolioWeight(i-2)+(Weight(i-1)-Weight(i-2))*Unit,-1*round(PortfolioNav(i-1)/(coef*Settle(i))),1*round(PortfolioNav(i-1)/(coef*Settle(i)))); % ��������λ
        Amount = PortfolioWeight(i-2)*coef*Settle(i);
        [FinalCost] = Cost(Comdty,PortfolioWeight(i-2),Amount);
        PortfolioNav(i,1) = PortfolioNav(i-1,1)+PortfolioWeight(i-2)*coef*(Settle(i)-Settle(i-1))+(PortfolioNav(i-1,1)-Amount/5)*(ChinaGC.data(find(ismember(ChinaGC.dates,Date(i))),1)/365)-FinalCost;
    end
end

%% ����ָ��
% ������
CumulativeReturn = PortfolioNav(end)/PortfolioNav(1)-1;
% �껯����
AnnualReturn = (1+CumulativeReturn)^(250/length(PortfolioNav))-1;
% ������
dret = price2ret(PortfolioNav);
StandardDeviation = std(dret)*sqrt(250);
% ���ձ���
SharpRatio = AnnualReturn/StandardDeviation;
% �س�
DrawDown = [];
for i = 1:length(PortfolioNav)
    DrawDown(i) = PortfolioNav(i)/max(max(PortfolioNav(1:i)))-1;
end
% ���س�
MaxDrawDown = min(min(DrawDown));

%% ��ͼ
subplot(3,1,1)
plot(Date,PortfolioNav)
dateaxis('x',2)
subplot(3,1,2)
plot(Date,DrawDown)
dateaxis('x',2)
subplot(3,1,3)
bar(Date,PortfolioWeight)
dateaxis('x',2)

%% ������
disp(['������ = ',num2str(CumulativeReturn)])
disp(['�껯���� = ',num2str(AnnualReturn)])
disp(['������ = ',num2str(StandardDeviation)])
disp(['���ձ��� = ',num2str(SharpRatio)])
disp(['���س� = ',num2str(MaxDrawDown)])
end