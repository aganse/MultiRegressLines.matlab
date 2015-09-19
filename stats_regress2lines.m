% Two-phase straight-line regression: fitting two straight lines to data.
% This script uses Monte Carlo analysis to compute solution statistics for
% the 2-line regression script regress2lines.m
% ----------------
% Andy Ganse, 2006, http://staff.washington.edu/aganse
% Applied Physics Laboratory, University of Washington, Seattle
% ----------------
%
% Usage:
%   With no arguments, run a demonstration using made-up data:
%      stats_regression2lines();
%   Or with arguments, specify x and y column vectors of data points:
%      stats_regression2lines(x,y);
%   All output simply goes to the screen (stdout).
%   If the data have more than 100 x,y pairs then the Monte Carlo (MC)
%   statistical analysis is skipped by default.  Nothing prevents you
%   from going into this script and reenabling the MC analysis for >100
%   points, merely concerns about computation time.
%
function stats_regress2lines(x_in,y_in)
clear all;

% Create the example data or bring in the actual data:
if nargin==0   %  <-- no input args means run the example
    % Create some synthetic, noisy data:
    % (note we construct this data to be sorted in X)
    %x=[1:2:100]';  % evenly spaced x locations
    x=sort(100*rand(50,1));  % random x locations
    % Some choices for "true" model m=[a1,b1,a2,b2,x0] for y=ax+b :
    %mtrue=[3,8,-.5,148,40]';
    mtrue=[5,-10,-2,168.5,25.5]';
    %mtrue=[3,0,-3,300,50]';
    %mtrue=[0,0,0,0,68]';  % note this will mess up horribly since all zeros
    icutoff=max(find(x<mtrue(5)));
    sigmatrue=10;
    ytrue=zeros(length(x),1);  % just dimensioning the variable
    y=zeros(length(x),1);  % just dimensioning the variable
    % calculate the values of the two connecting lines at the given x's:
    ytrue(1:icutoff)=mtrue(1)*x(1:icutoff)+mtrue(2);
    ytrue(icutoff+1:length(x))=mtrue(3)*x(icutoff+1:length(x))+mtrue(4);
    % add some independent Gaussian noise to them:
    noise=sigmatrue*randn(length(x),1);  % gaussian with mean=0 & stdev=sigmatrue
    y=ytrue+noise;
else   %  <-- inputting actual data vectors x and y
    % Check that x_in and y_in are column vectors, if not make them so:
    if size(x_in,2)~=1
        x_in=x_in';
    end
    if size(y_in,2)~=1
        y_in=y_in';
    end
    % Sort the data in X, just in case it's not already:
    data=sortrows([x_in,y_in]);
    x=data(:,1);
    y=data(:,2);
end


% Calc the parameters for the 2 regression lines:
[m, sumsqr, i_intersect, optG] = regress2lines(x,y);


% Linearize the problem at the solution to calculate covariances:
% calc sample stdev of data based on residuals (assuming constant stdev):
s=sqrt( sumsqr/(length(y)-4) );  % ie minus the 4 linear parameters 
% calc model covariance matrix (sqrts of diagonal are the model stdevs)
% for just the four line parameters, independent of the x_intercept param:
Cm=s^2*inv(optG'*optG);


% Just printing out all the results from above...
fprintf(1,'\nStatistics for two-phase straight-line regression:\n');  
fprintf(1,'Model is two lines:  y=a1*x+b1 and y=a2*x+b2, intersecting at x0.\n');
fprintf(1,'-------------------------------------------------------------------------\n');
if exist('mtrue','var')
    fprintf(1,'Playing with synthetic data; true model for m=[a1,b1,a2,b2,x0] is:\n');
    fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g \n',mtrue');
end
fprintf(1,'Note we assume that the data noise are Gaussian with constant unknown stdev.\n');
if exist('sigmatrue','var')
    fprintf(1,'Known population standard deviation for synthetic true data = %9.3g\n',sigmatrue);
end
fprintf(1,'Data sample standard deviation s computed from residuals = %9.3g\n', s);
fprintf(1,'Solution via regress2lines.m for m=[a1,b1,a2,b2,x0]:\n');
fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g  \n',m');
fprintf(1,'Associated stdevs for m=[a1 b1 a2 b2], computed from covariance matrix below:\n');
fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g  \n',sqrt(diag(Cm)));
fprintf(1,'Model covariance matrix for m=[a1 b1 a2 b2], computed from s and optimal G :\n');
fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g  \n',Cm);
fprintf(1,'Model correlation matrix for m=[a1 b1 a2 b2], computed from cov matrix :\n');
sigCm=sqrt(diag(Cm));
for i=1:4
    for j=1:4
        R(i,j)=Cm(i,j)/sigCm(i)/sigCm(j);
    end
end
fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g  \n',R);


% Exploring total, explained, and unexplained sum-sqrs:
ymodel=optG*m(1:4);
ymodelmean=mean(ymodel);
ydatamean=mean(y);
SSresids = sum( (y-ymodel).^2 );  % ie unexplained
SSmodel = sum( (ymodel-ymodelmean).^2 );  % ie explained
SStotal = sum( (y-ydatamean).^2 );
fprintf(1,'\n\nSStotal should equal SSmodel+SSresids, so SStotal-SSmodel-SSresids should = 0:\n');
fprintf(1,['SStotal = ' num2str(SStotal) '\n']);
fprintf(1,['SSmodel = ' num2str(SSmodel) '\n']);
fprintf(1,['SSresids = ' num2str(SSresids) '\n']);
fprintf(1,['SStotal-SSmodel-SSresids = ' num2str(SStotal-SSmodel-SSresids) '  <-- (consider: close enough to zero?)\n']);
fprintf(1,['Coefficient of determination R^2 = SSmodel/SStotal = ' num2str(SSmodel/SStotal) '\n\n\n']);


if length(x)>100  %  <-- too much for Monte Carlo then
% if you want to turn on the Monte Carlo analysis for more points, just
% increase this 100 here to a larger number.  But it might run for a long time!
    
    fprintf(1,'More than 100 data points, so skipping Monte Carlo analysis.\n');
    plot(x,y,'*'); hold on;
    plot([x(1),m(5)]',m(1)*[x(1),m(5)]'+m(2),'r',...
        [m(5),x(end)]',m(3)*[m(5),x(end)]'+m(4),'r','linewidth',2);
    hold off;
    return;
    
else  %  <-- do the Monte Carlo statistical analysis
    
    numMCruns=1e4;  % 1000 runs yields about 3% error on MC estimate,
                    % 10000 runs yields about 1% error on MC estimate.
                    % it would be preferable to do 1e4, 1e5, or even 1e6
                    % runs, but that likely takes way too long.
    
    % Now let's do Monte Carlo to estimate the stats of all five model params
    fprintf(1,'Starting %d Monte Carlo runs to compute Cm (each run calls regress2lines)\n', numMCruns);
    for j=1:numMCruns
        if mod(j,numMCruns/10)==0
            fprintf(1,'%2d%%... ',(100*j/numMCruns));
        end
        ynew(:,j)=ymodel+s*randn(length(ymodel),1);  % use estimated (sample) stats
        %ynew(:,j)=ytrue+10*randn(length(ytrue),1);  % use true (population) stats
        [MC_m(:,j), MC_sumsqr(j), MC_i_intersect(j), MC_G] = regress2lines(x,ynew(:,j));
        drawnow();  % (allows Ctrl-C to interrupt heavy calculations...)
    end
    fprintf(1,'\n\nFull-scale Monte Carlo means for m=[a1,b1,a2,b2,x0] via %g fwd problem runs:\n',numMCruns);
    fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g \n',mean(MC_m'));
    Cm=cov(MC_m');
    fprintf(1,'Associated standard deviations for m=[a1 b1 a2 b2 x0] based on cov matrix below:\n');
    fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g \n',sqrt(diag(Cm)));    
    fprintf(1,'Monte Carlo model cov matrix for m=[a1 b1 a2 b2 x0] via %g fwd problem runs:\n',numMCruns);
    fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g \n',Cm);
    fprintf(1,'Model correlation matrix for m=[a1 b1 a2 b2], computed from cov matrix :\n');
    sigCm=sqrt(diag(Cm));
    for i=1:5
        for j=1:5
            R(i,j)=Cm(i,j)/sigCm(i)/sigCm(j);
        end
    end
    fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g \n',R);


    % Plot the data, solution fit lines, and all the MC fit lines in background
    plot(x,y,'*');
    hold on;
    for i=1:numMCruns
        plot([x(1),MC_m(5,i)]',...
            MC_m(1,i)*[x(1),MC_m(5,i)]'+MC_m(2,i),...
            [MC_m(5,i),x(end)]',...
            MC_m(3,i)*[MC_m(5,i),x(end)]'+MC_m(4,i),'Color',[.75 .75 .75]);
    end
    plot(x,y,'*');
    plot([x(1),m(5)]',m(1)*[x(1),m(5)]'+m(2),'r',...
        [m(5),x(end)]',m(3)*[m(5),x(end)]'+m(4),'r','linewidth',2);
    title('blue points = data, red line = solution, grey lines = Monte Carlo runs');
    hold off;


    % Now let's compare x_intersect row & column in Cm above to that calculated
    % here based only upon stats of m(1-4)={a1,b1,a2,b2} and the constraint:
    % x_intersect = ( m(4)-m(2) ) / ( m(1)-m(3) )
    % It's not a linear relationship between x_intersect and m(1-4), because
    % it's not a linear combination of model parameters (a1,b1,a2,b2).
    % So again we use Monte Carlo, but this one is much faster than above.
    % Let's explore the results of this to see if they may be useful when
    % the above is not available due to much longer calculation time...?
    %
    % Simulate samples of m(1-4) by simulating zero-mean deviates based on Cm
    % and Chol decomp, adding those onto the mean m soln.  As above, we build
    % up a histogram of samples and from that find stats.
    % Use same number of MC runs as the above MC (easier to compare that way)
    fprintf(1,'\n\nWait another minute for "alternate/cheap" Monte Carlo (faster than above)...\n');
    L=chol(Cm(1:4,1:4));
    U=randn(numMCruns,4);
    Z=U*L;
    CZ=cov(Z);
    a1=m(1)+Z(:,1);
    b1=m(2)+Z(:,2);
    a2=m(3)+Z(:,3);
    b2=m(4)+Z(:,4);
    x0=( b2-b1 )./( a1-a2 );
    Cm2=cov([a1,b1,a2,b2,x0]);
    fprintf(1,'\n"Alternate/cheap" Monte Carlo means for m=[a1,b1,a2,b2,x0] via G and Chol decomp:\n');
    fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g \n',mean([a1,b1,a2,b2,x0]));
    fprintf(1,'Associated stdevs for m=[a1,b1,a2,b2,x0] based on cov matrix below:\n');
    fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g \n', sqrt(diag(Cm2)));
    fprintf(1,'Monte Carlo model cov matrix for m=[a1,b1,a2,b2,x0] via G and Chol decomp:\n');
    fprintf(1,'%9.4f  %9.4f  %9.4f  %9.4f  %9.4f \n',Cm2);
    fprintf(1,'Model correlation matrix for m=[a1 b1 a2 b2], computed from cov matrix :\n');
    sigCm=sqrt(diag(Cm2));
    for i=1:5
        for j=1:5
            R(i,j)=Cm2(i,j)/sigCm(i)/sigCm(j);
        end
    end
    fprintf(1,'%9.4g  %9.4g  %9.4g  %9.4g %9.4g \n',R);

end