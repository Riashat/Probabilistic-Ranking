load('tennis_data.mat')
% [kk, ii] = sort(G, 'descend');


np = 107;

count_won = zeros(np,1);
count_lost = zeros(np,1);
ratio = zeros(np,1);
for i = 1:107
   count_won(i) = sum(G(:,1)==i);
   count_lost(i) = sum(G(:,2)==i);    
   ratio(i) = count_won(i)/ ( count_won(i) + count_lost(i) );
end

[kk, ii] = sort(ratio, 'descend');

np = 107;
barh(kk(np:-1:1))
set(gca,'YTickLabel',W(ii(np:-1:1)),'YTick',1:np,'FontSize',10)
axis([0 1 0.5 np+0.5])

title('Average Probability (Empirical Ratio) of each player winning', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Empirical Ratio (Average Probability) of Player to Win', 'FontSize', 15);
ylabel('Player Name', 'FontSize', 15);
