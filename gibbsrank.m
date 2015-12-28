clear all
clc
load tennis_data

randn('seed',27); % set the pseudo-random number generator seed

M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

pv = 0.5*ones(M,1);           % prior skill variance 

w = zeros(M,1);               % set skills to prior mean
independent_samples = zeros(M,100); %because we want 100 samples - each 107x1 vector is a sample of 107 players

iterations = 1000;
P_1 = zeros(length(iterations), 1);
P_7 = zeros(length(iterations), 1);
P_10 = zeros(length(iterations), 1);
P_15 = zeros(length(iterations), 1);
P_100 = zeros(length(iterations), 1);
P_107 =  zeros(length(iterations), 1);

for iters = 1:iterations

  % First, sample performance differences given the skills and outcomes
  
  t = nan(N,1); % contains a t_g variable for each game
  for g = 1:N   % loop over games
    s = w(G(g,1))-w(G(g,2));  % difference in skills
    t(g) = randn()+s;         % performace difference sample
    while t(g) < 0  % rejection sampling: only positive perf diffs accepted
      t(g) = randn()+s; % if rejected, sample again
    end
  end 
 
  
  % Second, jointly sample skills given the performance differences
  
  m = nan(M,1);  % container for the mean of the conditional 
                 % skill distribution given the t_g samples
  for g = 1:M
   m(g) = t'*((g==G(:,1)) - (g==G(:,2)));
  end
  
  iS = zeros(M,M); % container for the sum of precision matrices contributed
                   % by all the games (likelihood terms)
   %build the iS matrix
                
    for i = 1:M
        for j = 1:i
            if i==j  
                iS(i,j) = sum(i==G(:,1)) + sum(i==G(:,2));
            else
                iS(i,j) = -sum((i==G(:,1)).*(j==G(:,2))+(i==G(:,2)).*(j==G(:,1)));
                iS(j,i) = iS(i,j);
            end
        end
    end
           
  iSS = diag(1./pv) + iS; % posterior precision matrix
  % prepare to sample from a multivariate Gaussian
  % Note: inv(M)*z = R\(R'\z) where R = chol(M);
  iR = chol(iSS);  % Cholesky decomposition of the posterior precision matrix
  mu = iR\(iR'\m); % equivalent to inv(iSS)*m but more efficient
    
  % sample from N(mu, inv(iSS))
  w = mu + iR\randn(M,1);
    
  P_1(iters)=w(1);
  P_7(iters)=w(5);
  P_10(iters)=w(10);
  P_15(iters)=w(15);
  P_100(iters)=w(100);
  P_107(iters) = w(107);
  
  %to make sure that samples are roughly 15 units apart to be independent
  if mod(iters+9,10)==0     
    independent_samples(:,(iters+9)/10)=w;
  end
  
end

autocorr(independent_samples(4,:))

figure(1)

subplot(2,2,1)
plot(1:iterations, P_1, 'r.')
title('Sampled Player Skill for Player 1')

subplot(2,2,2)
plot(1:iterations, P_7, 'b.')
title('Sampled Player Skill for Player 7')

subplot(2,2,3)
plot(1:iterations, P_10, 'g.')
title('Sampled Player Skill for Player 10')

subplot(2,2,4)
plot(1:iterations, P_15, 'r.')
title('Sampled Player Skill for Player 15')
suptitle('Pseudo Random Seed = 27')


figure (2)

subplot(2,2,1)
plot(1:iterations, P_100, 'm.')
title('Sampled Player Skill for Player 100')

subplot(2,2,2)
plot(1:iterations, P_107, 'k.')
title('Sampled Player Skill for Player 107')
suptitle('Pseudo Random Seed = 40')

k = mean(independent_samples,2);    %mean along the rows 
prob = zeros(M,size(independent_samples,2));
Sum = 0;

%Calculate average probaility
for P1=1:M
   for P2=1:M
       for l = 1:size(independent_samples,2)
       prob(P1, l) = normcdf(independent_samples(P1,l) - independent_samples(P2,l));
       end
   end
end

mean_prob= mean(prob,2);

%Plot this ranking system 
[kk, ii] = sort(mean_prob, 'descend');

np = 107;
figure(3)
barh(kk(np:-1:1), 'r')
set(gca,'YTickLabel',W(ii(np:-1:1)),'YTick',1:np,'FontSize',6)
axis([0 1 0.5 np+0.5])
title('Probabilistic Ranking - Average Probability of Winning based on Gibbs Sampling', 'FontSize', 13, 'FontWeight', 'bold')
xlabel('Probability of a player winning against other players', 'FontSize', 12);
ylabel('Player Name', 'FontSize', 12);



%TABLE 4x4
P1 = [16,1,5,11];
player_skills = independent_samples(P1, :);
X = zeros(4,4);
for q = 1:4
    for w = 1:4
        X(q,w) = mean(normcdf( player_skills(q,:) - player_skills(w,:) ));
    end
end
    



higher_skills = zeros(4,4);
for q = 1:4
    for w = 1:4
        higher_skills(q,w) = mean(player_skills(q,:)>player_skills(w,:));
    end
end






