% make a bar plot from vector P and annotate with player names from W
load('tennis_data.mat')
[kk, ii] = sort(G, 'descend');

np = 107;
barh(kk(np:-1:1))
set(gca,'YTickLabel',W(ii(np:-1:1)),'YTick',1:np,'FontSize',8)
axis([0 1 0.5 np+0.5])

%%%

count_won = zeros(np,1);
count_lost = zeros(np,1);
ratio = zeros(np,1);
for i = 1:np
   count_won(i) = sum(kk(:,1)==i);
   count_lost(i) = sum(kk(:,2)==i);    
   ratio(i) = count_won(i)/ ( count_won(i) + count_lost(i) );
end

players = [1:1:107]';
ratio_players = [players, ratio];

A = sortrows(ratio_players, 2);
