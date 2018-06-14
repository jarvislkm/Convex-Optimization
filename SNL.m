%% basic SNL problem simulation
clear all;
anchor_number = 5;
sensor_number = 45;
global_record = [];
anchor = rand(anchor_number,2)*100;
sensor = rand(sensor_number,2)*100;
average_time = 1;

p_start = 30;
p_end = 30;
p_step = 5;

for threshold = p_start:p_step:p_end
    record = [];
for times = 1:1:average_time
%% 
% figure
% plot(anchor(:,1), anchor(:,2),'bo', 'LineWidth',3);
% hold on;
% plot(sensor(:,1), sensor(:,2),'ro', 'LineWidth',3);
% hold on;
% xlim([0,100]);
% ylim([0,100]);
% title('position of anchor and sensor');
%%
kernal = sensor*sensor';
Z = [eye(2) sensor';
     sensor kernal];
indicator = eye(sensor_number);
%%
M_ss = zeros(sensor_number+2, sensor_number+2, sensor_number, sensor_number);
for i = 1:1:sensor_number
    for j = 1:1:sensor_number
        M_ss(3:end, 3:end, i, j) = (indicator(i,:)-indicator(j,:))'*(indicator(i,:)-indicator(j,:));
    end
end
M_sa = zeros(sensor_number+2, sensor_number+2, anchor_number, sensor_number);
for i = 1:1:anchor_number
    for j = 1:1:sensor_number
        M_sa(:,:,i,j) = [anchor(i,:)'; indicator(j,:)']*[anchor(i,:) indicator(j,:)];
    end
end
%% Ground Truth
Distance_ss = zeros(sensor_number, sensor_number);
for i = 1:1:sensor_number
    for j = 1:1:i-1
        Distance_ss(j,i) = trace(M_ss(:,:,i,j)*Z);
        Distance_ss(i,j) = Distance_ss(j,i);
    end
end


Distance_sa = zeros(anchor_number, sensor_number);
for i = 1:1:anchor_number
    for j = 1:1:sensor_number
        Distance_sa(i,j) = trace(M_sa(:,:,i,j)*Z);
    end
end
%% hist
Distance_sa_vec = reshape(Distance_sa, size(Distance_sa,1)*size(Distance_sa, 2), 1);
Distance_ss_vec = reshape(Distance_ss, size(Distance_ss,1)*size(Distance_ss, 2), 1);

% figure 
% max_sa = max(max(hist(Distance_sa)));
% hist(Distance_sa)
% hold on
sa_threshold = prctile(Distance_sa_vec, threshold);
% plot([sa_threshold sa_threshold], [0, max_sa], 'r.-');
% title('hist of anchor-sensor distance');
% 
% figure
% max_ss = max(max(hist(Distance_ss)));
% hist(Distance_ss);
% hold on
ss_threshold = prctile(Distance_ss_vec, threshold);
% plot([ss_threshold ss_threshold], [0, max_ss], 'r.-');
% title('hist of sensor-sensor distance');
%% CVX
n = sensor_number;
m = anchor_number;
cvx_begin
   cvx_precision high
   variable X(n+2,n+2) symmetric;
   minimize( 0 )
   subject to:
    for i = 1:2
        X(i,i) == 1;
    end
    for i = 1:n
        for j = 1:n
            if Distance_ss(i,j)<= ss_threshold
                abs(trace(M_ss(:,:,i,j)*X) - Distance_ss(i,j)) <= 10;
            end
        end
    end
    for i = 1:m
        for j = 1:n
            if Distance_sa(i,j)<= sa_threshold
                abs(trace(M_sa(:,:,i,j)*X) - Distance_sa(i,j)) <= 10;
            end
        end
    end
    X == semidefinite(n+2);
cvx_end
%% measurement
diff = trace((X-Z)*(X-Z)');
disp('The difference of reconstruction matrix is: ')
disp(diff)

%% reconstruction without rank-2 approx
% figure
% plot(sensor(3:end,1), sensor(3:end,2),'r*', 'LineWidth',5);
% hold on;
% plot(X(3:end,1), X(3:end,2),'Bo', 'LineWidth',2);
% hold on
% xlim([0,100]);
% ylim([0,100]);
% title('position of anchor and sensor');

%% Rank 2 approximation
if isnan(max(max(X)))
    diff_2 = diff;
else
    [U,S,V] = svds(X);
    Z_2 = U(:,1:2) * S(1:2,1:2) * V(:,1:2)';
    pos = Z_2(3:end,1:2);
    diff_2 = trace((pos-sensor)*(pos-sensor)');
end
disp('The difference of reconstruction matrix(rank2 approx) is: ')
disp(diff_2)
record = [record diff_2];
end
figure
plot(sensor(:,1), sensor(:,2),'r*', 'LineWidth',5);
hold on
plot(Z_2(3:end,1), Z_2(3:end,2),'bo', 'LineWidth',2);
hold on
plot(anchor(:,1), anchor(:,1),'c*', 'LineWidth',5);
xlim([0,100]);
ylim([0,100]);
legend('position', 'reconstruction','anchor');
title('Position of anchor and sensor');

global_record = [global_record; record];
end


