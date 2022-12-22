%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate functional groups
function [funcGroups, names_funcGroups] = calc_funcGroups(proteins) 
% Protein order
% 1: APAF-1
% 2: BAK (BAK1)
% 3: BAX
% 4: BCL2
% 5: BCLXL (BLC2L1)
% 6: BID
% 7: PC3 (CASP3)
% 8: PC9 (CASP9)
% 9: MCL1
% 10: SMAC
% 11: XIAP
% 12: BIM
% 13: PUMA (BBC3)
% 14: NOXA

    APAF_C9 = proteins(:,1) .* proteins(:,8);
    BAK_BAX = proteins(:,2) + proteins(:,3);
    MCL1_BCL2_BCLXL = proteins(:,9) + proteins(:,4) + proteins(:,5);
    BID = proteins(:,6);
    XIAP_C3 = proteins(:,11) ./ proteins(:,7)
    SMAC = proteins(:,10);
    BIM = proteins(:,12);
    PUMA = proteins(:,13);
    NOXA = proteins(:,14);

    funcGroups = [APAF_C9 XIAP_C3 SMAC BAK_BAX MCL1_BCL2_BCLXL BIM BID PUMA NOXA];
    names_funcGroups = {'APAF-1*Caspase-9' 'XIAP/C3' 'SMAC' 'BAK+BAX' 'MCL1+BCL2+BCLXL' 'BIM' 'BID' 'PUMA' 'NOXA' };

end