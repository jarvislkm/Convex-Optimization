%% Initialize
clear all;
num_iter = 5;
% Nlist = [5, 10, 20, 40]
Nlist = [5, 10, 20, 40, 80]
result_truths = zeros([size(Nlist,2),num_iter]);
result_SDRs = zeros([size(Nlist,2),num_iter]);
%% SDR QCQP 
for idx_N = 1:size(Nlist,2)
    N = Nlist(idx_N)
    for idx_iter = 1: num_iter
        % Adjacency Matrix
        A = randi(2,N,N) - 1;
        A = A - tril(A,-1) + triu(A,1)';
        % Laplacian Matrix
        D = diag(A * ones(N,1));
        L = D - A;
        % find the ground truth
        maxCut = 0;
        for i = 1:2^N
            x = de2bi([i-1], N);
            x = x *2 -1;
            currCut = x*L*x'/4;
            if currCut > maxCut
                maxCut = currCut;
                result = x;
            end
        end
        % CVX
        cvx_begin
        %   cvx_precision high
           variable X(N,N) symmetric;
           maximize( trace(L*X)/4 )
           subject to:
           for i = 1:N
                e = zeros([1,N]);
                e(i) = 1;
                E = e'*e;
                trace(E*X) == 1;
           end
            X == semidefinite(N);
        cvx_end
        % SDR result
        [U,S,V] = svd(X);
        X_1 = U(:,1) * S(1,1) * V(:,1)';

        result_truth = maxCut
        result_SDR = trace(L*X)/4
        result_truths(idx_N, idx_iter) = result_truth;
        result_SDRs(idx_N, idx_iter) = result_SDR;
    end
end

%%
result_truths;
result_SDRs;
result_SDRs./result_truth
bound_ratio = mean(result_SDRs./result_truths, 2)

figure
plot(log(Nlist),bound_ratio,'-o')
title('Max cut : log N vs (upper bound/groud truth','FontSize',20)
xlabel('log N','FontSize',20)
ylabel('ratio','FontSize',20)

