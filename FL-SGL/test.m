clear, clc;

mex('dpotrf_.c', '-lmwblas','-lmwlapack');
mex('dpotrs_.c', '-lmwblas','-lmwlapack');


load('test_data.mat');
Xtrain = Data.Xtrain;
Xtest = Data.Xtest;
Ytrain = Data.Ytrain;
Ytest = Data.Ytest;

Iters = length(Xtrain);
Times = length(Xtrain{1});
Tasks = size(Ytrain{1}{1},2);
Dimen = size(Xtrain{1}{1},2);

Ncross = 5;
Print = 0;
sigma = 10;

% set opts
opts.Print = Print;
opts.MaxIter = 1000;
opts.Tol = 1e-3;
opts.rho = 10;
opts.Dimen = Dimen;
opts.Times = Times;
opts.sigma = sigma;

lambda_1 = 10; 
lambda_2 = 10;
lambda_3 = 10;

admm_mets = {'two', 'multi'};

for j = 1:length(admm_mets),
    fprintf('ADMM %s-block\n', admm_mets{j});
    for i = 1:Iters
        Xtraini = Xtrain{1,i};
        Xtesti = Xtest{1,i};
        Ytraini = Ytrain{1,i};
        Ytesti = Ytest{1,i};
        for t = 1:Times
            for k = 1:Tasks
                Ytraint{1,k}{1,t} = Ytraini{1,t}(:,k);
                Ytestt{1,k}{1,t} = Ytesti{1,t}(:,k);
            end
        end

        for k = 1:Tasks
            fprintf('Iters,%d,Tasks,%d\n',i,k);
            Ytraink = Ytraint{1,k};
            Ytestk = Ytestt{1,k};
            if strcmp(admm_mets{j}, 'two')
                [W] = FLADMM_TB(Xtraini, Ytraink, lambda_1, lambda_2,lambda_3, opts);
            else
                [W] = FLADMM_MB(Xtraini, Ytraink, lambda_1, lambda_2,lambda_3, opts);
            end
        end
    end
end





