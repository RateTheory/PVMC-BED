%% Select Optimal Columns with Cone Method
% Selin Bac, Usama Saleem, Shaama M. Sharada University of Southern
% California, 2025
clc;
clear all;
close all;
pyenv("Version",'3.10', 'ExecutionMode','OutOfProcess')


% Change this block
K = 7;                  % Total Number of columns to sample in addition to stationary points
k=  0;                  % number of columns via bond distance

writematrix(k, 'k.txt') % Number of columns to sample with bond distance
method='bond_distance';  % Choose method
%method= 'euclidean' 

pyrunfile('bond_selection.py')

system = 'Ir';      % Reaction to recover (see MatrixMATs)
half_reaction = false;  % Must be true when using symmetric half-reactions (Sn2)
partition = false;      % Enables roughly even sampling of product and reactant regions.
generate_plots = false;

writelines(strcat(system), 'sys.txt')

prefix = 'MatrixMATs/';
load(strcat(prefix,system, '.mat'));



strue=s;
sbd= importdata('bd_s.txt');                                               % // example data
[~, ind] = min(abs(bsxfun(@minus, strue(:), sbd(:).')), [], 1); %'// compute all differences
                                                            % // and find index of minimizer
result = strue(ind); 
disp(sbd')
disp(result)
Ibd=zeros(1, length(sbd));
for n=1:length(sbd)
    element=result(n);
    i=find(s==element);
    Ibd(n)=i;
end

i = cone_method_selin(Xtrue,s, K, half_reaction, partition, Ibd);
j= cone_method_mod(Xtrue,s, K, half_reaction, partition, Ibd);
m= cone_method(Xtrue,s, K-k, half_reaction, partition);

Ibd;

I=union(j, Ibd)


%I=cone_method(Xtrue,s, K, half_reaction, partition);


%% Perform PVMC and HVMC

% Declared Parameters for HVMC
options.gamma0=0;
options.eta=1.01;
options.p=1; %0.5
options.eigtol=1e-7;
options.epsilon=0;
options.d = 2; 
options.niter = 5000;
options.gammain=1e-16;
options.c=1;
options.exit_tol=1e-8;
options.m=size(Xtrue,1);
options.n=size(Xtrue,2);
options.stepsize=0.1;
options.itersvt=100;
options.tau=0.001;
options.maxCol=5;

%Declared Parameters for PVMC
options_QS.d = 2; 
options_QS.p = 1;
options_QS.exit_tol = 1e-8;
options_QS.niter = 5000;
options_QS.polynomial_degree = 4; % Default 4
options_QS.lambda = 0.997; %[0,1] Determines strength of polynomial interpolation
options_QS.gammamin=1e-16;

m=size(Xtrue,1);
n=size(Xtrue,2);
R=eye(size(Xtrue,2),size(Xtrue,2));

% Sample Columns
sampmask_c = false(m,n);
for j=1:length(I)
    sampmask_c(:,I(j))=true(m,1);
end

sigma=0; % Noise parameter
m = size(Xtrue,1);
n = size(Xtrue,2);

samples_c = Xtrue(sampmask_c)+sigma*randn(size(Xtrue(sampmask_c)));
noise=sigma*randn(size(Xtrue(sampmask_c)));
Xinit_c=zeros(m,n);
Xinit_c(sampmask_c) = samples_c;

% Perform PVMC
[Xvmc2,error2,error_qs2,Qi,Qf,S,iter_info] = pvmc_step(Xinit_c,Xtrue,sampmask_c,I,s,options_QS);

% Perform HVMC to Compare
[Xvmc,cost,update,error] = vmc_step(Xinit_c,sampmask_c,samples_c,options,Xtrue);

error_t=norm(Xvmc-Xtrue,'fro')/norm(Xtrue,'fro');
error_t2 =norm(Xvmc2-Xtrue,'fro')/norm(Xtrue,'fro');
error_qs_t=norm(Qf*S-Xtrue,'fro')/norm(Xtrue,'fro');

fprintf('HVMC NRMSE = %1.2e\n',mean(error_t));
fprintf('PVMC X NRMSE = %1.2e\n',mean(error_t2));
fprintf('QS NRMSE = %1.2e\n',mean(error_qs_t));

Xpredictpoly = Qi*S;
Xpredictpoly(:,I) = Xtrue(:,I);
Xpolynorm = norm((Xvmc2-Xpredictpoly)./norm(Xpredictpoly,'fro'),'fro');
QS_original_error = norm((Xtrue-Xpredictpoly)./norm(Xtrue,'fro'),'fro');
fprintf('Polynomial Interpolation NRMSE  = %1.2e\n',QS_original_error);

save(strcat('ResultMATs/results_',num2str(k), system, '.mat'), ...
    'Xinit_c','Xvmc','Xvmc2','Xtrue','s', ...
    'system','K', 'I', 'options', 'options_QS')

save(strcat('results',num2str(k), '.mat'), ...
    'Xinit_c','Xvmc','Xvmc2','Xtrue','s', ...
    'system','K', 'I', 'options', 'options_QS')



%% Plot Results
plot_results; 
%plot_results
