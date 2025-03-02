function [obelezja] = izdvoji_obelezja(bin)
    [N, M] = size(bin);
    sums = sum(bin);
    [max_sum, max_ind] = max(sums);
    [row,col] = find(bin);
    col_min = min(col);
    obelezja1 = max_ind-col_min;
    ivice = edge(bin(:,1:max_ind), 'Canny');
    obelezja2 = sum(sum(ivice));

    
%     figure();
%     subplot(1,2,1);
%     imshow(bin(:,1:max_ind));
%     title('Izdvajanje obeležja');
%     subplot(1,2,2);
%     imshow(ivice);
%     title('Izdvajanje obeležja');

    obelezja = [obelezja1, obelezja2];
end