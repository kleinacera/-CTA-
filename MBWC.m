function [PortfolioNav,Date] = MBWC(Comdty,Begin,End,Account)
% Momentum �������Ƶ�CTA����
% ����
%     Comdty   �ڻ�����,str��ʽ,��'RB.SHF'
%     Begin    �ز⿪ʼʱ��,str��ʽ,��'2016-01-01'
%     End      �ز����ʱ��,str��ʽ,��'2017-12-31'
%     coef     ���ڻ���Լÿһ�ֵ������������Ƹ�1��10�֣�coef=10
%     Account  �˻���ʼ���

%% ��wind��ȡ����
load FutureBaseData
w = windmatlab;
[~,~,~,Date] = w.tdays(Begin,End,'TradingCalendar=DCE');
Open = w.wsd('CF.CZC','open',Begin,End);
High = w.wsd('CF.CZC','high',Begin,End);
Low = w.wsd('CF.CZC','low',Begin,End);
Close = w.wsd('CF.CZC','close',Begin,End);
coef = FutureBaseData.mul(find(ismember(FutureBaseData.code,Comdty)));

%% �������
N = 5;M = 20; % N ���ھ��ߡ���ʵ�����ʺ����ͻ��ֵ�Ļزⳤ��; M ���ھ��ߵ�ʱ�䳤��
[Short,Long]=movavg(Close,N,M,'e');
ATR = AverageTrueRange(High,Low,Close,N); % �����ʵ���ڲ���

%% �����ź�ʶ��
% % ����Ϊ1,2,3....������Ϊ-1,-2,-3,...���ղ�Ϊ0 ���е�����Ϊ��λ����
Weight = zeros(length(Open),1); %����
for i = M:length(Open)
    if Weight(i-1) == 0
        if Short(i)>Long(i) && High(i)==max(High(i-M+1:i))
            Weight(i) = 1;
        elseif Short(i)<Long(i) && Low(i)==min(Low(i-M+1:i))
            Weight(i) = -1;
        end
    else
        if Weight(i-1) >= 1
            if Weight(i-1)>Weight(i-2)
                BuyStand = Open(i); %����۸�
            end
            if Short(i)>=Long(i)
                if High(i)>=BuyStand+0.5*ATR(i) %������
                    Mul = floor((High(i)-BuyStand)/(0.5*ATR(i)));
                    Weight(i) = Weight(i-1)+Mul;
                elseif Low(i)<=BuyStand-2*ATR(i)||Low(i)==min(Low(i-M+1:i))
                     Weight(i) = 0;
                else
                    Weight(i) = Weight(i-1);
                end
            else
                Weight(i) = 0;
            end
        else
            if Weight(i-1)<Weight(i-2)
                SellStand = Open(i); %����۸�
            end
            if Short(i)<=Long(i)
                if Low(i)<=SellStand-0.5*ATR(i)
                    Mul = floor((SellStand-Low(i))/(0.5*ATR(i)));
                    Weight(i) = Weight(i-1)-Mul;
                elseif High(i)>=SellStand+2*ATR(i-1)||High(i)==max(High(i-M+1:i))
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
            PortfolioNav(i,1) = PortfolioNav(i-1,1);
        else %ƽ��
            Amount = PortfolioWeight(i-1)*coef*Open(i);
            [FinalCost] = Cost(Comdty,PortfolioWeight(i-2),Amount);
            PortfolioNav(i,1) = PortfolioNav(i-1,1)+PortfolioWeight(i-2)*coef*(Open(i)-Open(i-1))-FinalCost;
        end
    else %����
        Unit = round((0.01*PortfolioNav(i-1))/(ATR(i-1)*coef)); % �����λ�Ӳ���
        PortfolioWeight(i-1) = cap(PortfolioWeight(i-2)+(Weight(i-1)-Weight(i-2))*Unit,-10*round(PortfolioNav(i-1)/(coef*Open(i))),10*round(PortfolioNav(i-1)/(coef*Open(i)))); % ��������λ
        Amount = PortfolioWeight(i-1)*coef*Open(i);
        [FinalCost] = Cost(Comdty,PortfolioWeight(i-2),Amount);
        PortfolioNav(i,1) = PortfolioNav(i-1,1)+PortfolioWeight(i-2)*coef*(Open(i)-Open(i-1))-FinalCost;
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