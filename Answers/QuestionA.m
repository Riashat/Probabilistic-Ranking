load('tennis_data.mat')
[kk, ii] = sort(G, 'descend');

np = 107;

count_won = zeros(np,1);
count_lost = zeros(np,1);
ratio = zeros(np,1);
for i = 1:np
   count_won(i) = sum(kk(:,1)==i);
   count_lost(i) = sum(kk(:,2)==i);    
   ratio(i) = count_won(i)/ ( count_won(i) + count_lost(i) );
end
% 
% players = [1:1:107]';
% ratio_players = [players, ratio];
% 
% A = sortrows(ratio_players, 2);


