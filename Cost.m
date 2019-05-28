function [FinalCost] = Cost(Comdty,Weight,Amount)
% 求出交易成本
% 输入
%      Comdty    期货品种
%      Weight    手数
%      Amount    交易金额
% 输出
%      FinalCost 总交易成本

load FutureBaseData
Fee_amount = FutureBaseData.fee_amount(find(ismember(FutureBaseData.code,Comdty),1));
Fee_lot = FutureBaseData.fee_lot(find(ismember(FutureBaseData.code,Comdty),1));
FinalCost = Amount*Fee_amount+Weight*Fee_lot;

end