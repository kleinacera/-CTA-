function [FinalCost] = Cost(Comdty,Weight,Amount)
% ������׳ɱ�
% ����
%      Comdty    �ڻ�Ʒ��
%      Weight    ����
%      Amount    ���׽��
% ���
%      FinalCost �ܽ��׳ɱ�

load FutureBaseData
Fee_amount = FutureBaseData.fee_amount(find(ismember(FutureBaseData.code,Comdty),1));
Fee_lot = FutureBaseData.fee_lot(find(ismember(FutureBaseData.code,Comdty),1));
FinalCost = Amount*Fee_amount+Weight*Fee_lot;

end