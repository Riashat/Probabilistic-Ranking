% Created by: Dima Damen
% Date: December 2008

% This function takes a probability density function, and draws a set of
% samples of size n and plots the histogram of those samples along with the
% original pdf showing convergence.

% examples of some pdf you can experiment with are:
% pdf = @(x) 0.3*normpdf (x,0.3,0.05)+0.7*normpdf(x,0.7,0.2);
% pdf = @(x) normpdf (x, 4, 3);
% pdf = @(x) 0.1*normpdf (x,0.3,0.1)+0.3*normpdf(x,0.6,0.3)+0.6*normpdf(x,0.8,0.2);
% pdf = @(x) exppdf (x, 2);
% pdf = @(x) 0.4*exppdf(x,2) + 0.6*normpdf (x, 4, 2);
% pdf = @(x) betapdf (x, 2, 1);
% pdf = @(x) gampdf (x,1.5,2);

function MCMC (pdf, nsamples)

    %proppdf is the distribution from which you can sample/propose (Q)
    delta = 1;
    proppdf = @(x,y) unifpdf (y-x, -delta, 1);
    proprnd = @(x) x + rand*2*delta - delta;

    % metropolis-hastings algorithm 
    x = mhsample (0.4, nsamples, 'pdf', pdf, 'proppdf', proppdf, 'proprnd', proprnd);

    % d is the domain of the function you want to run
    d = min(x):0.01:max(x);

    % this draws the actual (hidden) distribution
    d_x = pdf (d);

    n = hist (x,50);
    prod = max(n)/max(d_x);
    yVal = d_x * prod;
    figure; hist(x,50);
    hold on;
    plot (d, yVal, 'r-', 'linewidth', 3);
end