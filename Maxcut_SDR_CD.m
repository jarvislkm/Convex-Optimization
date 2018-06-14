%% Initialize
clear all;
num_iter = 20;
Nlist = [4, 8, 16, 20]
result_truths = zeros([size(Nlist,2),num_iter]);
result_SDRs = zeros([size(Nlist,2),num_iter]);
result_CDs = zeros([size(Nlist,2),num_iter]);
p=0.50
%% SDR QCQP 
for idx_N = 1:size(Nlist,2)
    N = Nlist(idx_N)
    for idx_iter = 1: num_iter
        idx_iter
        % Adjacency Matrix
        A = randi(2,N,N) - 1;
        Mask = rand(N,N)>=1-p;
        A = A.*Mask;
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
        maxCut_truth = maxCut;
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

        % Coordinate descent
        x = randsample([-1,1], N,true);
        change =true;
        maxCut = x*L*x'/4;
        while change == true
            change = false;
            for i = 1:N
                x(i) = -1 * x(i);
                currCut =  x*L*x'/4;
                if currCut > maxCut
                    maxCut = currCut;
                    change =true;
                else
                    x(i) = -1 * x(i);
                end
            end
        end
        maxCut_CD = maxCut;
        
        result_truth = maxCut_truth
        result_SDR = trace(L*X)/4
        result_CD = maxCut_CD
        result_truths(idx_N, idx_iter) = result_truth;
        result_SDRs(idx_N, idx_iter) = result_SDR;
        result_CDs(idx_N, idx_iter) = result_CD;
    end
end

%%
result_truths;
result_SDRs;
result_SDRs./result_truth;
bound_ratio = mean(result_SDRs./result_truths, 2)
CD_ration = mean(result_CDs./result_truths,2)
figure
plot(log(Nlist),bound_ratio,'r-o',log(Nlist),CD_ration,'b-o' )
title('Max cut : log N vs (SDR relax bound/GT, CD /GT)','FontSize',20)
xlabel('log N','FontSize',20)
ylabel('ratio','FontSize',20)
legend({'SDR relaxation bound ratio','CD solver result ration'},'Fontsize',15);

%%

figure
plot(log(Nlist),mean(result_SDRs,2),'r-o'...
    ,log(Nlist),mean(result_truths,2),'g-o',...
    log(Nlist),mean(result_CDs,2),'b-o')
title('Max cut : log N vs Maxcut','FontSize',20)
xlabel('log N','FontSize',15)
ylabel('Maxcut','FontSize',15)
legend({'SDR relaxation bound','Grouth Truth','CD solver'},'Fontsize',15);



