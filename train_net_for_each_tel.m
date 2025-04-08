% This script is used to fit a BP neural network using the collected data to replace interactive calculations.

% Use the data collected from information_collection.m to train a neural network.
inputFileName = './00_result/net_fit_data/fittingDataResult-Tel1.mat'; 
train_data = importdata(inputFileName);

[~, name, ~] = fileparts(inputFileName); 
number = regexp(name, 'Tel(\d+)', 'tokens'); 
number = number{1}{1}; 

% Randomly shuffle the data
rowrank1 = randperm(size(train_data, 2));
ptrain = train_data(:, rowrank1);
input_data = ptrain(1:8, :);
output_data = ptrain(9, :);

% Split the training set and validation set
cv = cvpartition(size(input_data, 2), 'HoldOut', 0.2); % 80% for train，20% for test
idx = cv.test;

input_train = input_data(:, ~idx);
output_train = output_data(:, ~idx);
input_val = input_data(:, idx);
output_val = output_data(:, idx);

disp('The data has been processed.')

% train
net = newff(input_train, output_train, [15,30,15,5]);
net.trainFcn = 'trainbr';                % trainlm/trainbr/trainscg 。
net.trainParam.epochs = 350;          
net.trainParam.show= 2;                
net.trainParam.lr = 0.06;          
net.trainParam.goal = 0.0001;          
net.trainParam.max_fail = 50;         

net = train(net, input_train, output_train);

% save fit net;
saveFileName = ['./00_result/net_fit_data/Tel' number 'fitting_net']; 
save(saveFileName, 'net'); 

% test
an_train = sim(net, input_train);
an_val = sim(net, input_val);
error_train = an_train - output_train;
error_val = an_val - output_val;

disp('train_error:');
disp(error_train);
disp('Validation_error:');
disp(error_val);

figure(1);
plot(input_train(8, :), output_train, 'bo'); 
hold on;
plot(input_train(8, :), an_train, 'r*'); 
plot(input_val(8, :), output_val, 'go'); 
plot(input_val(8, :), an_val, 'm*'); 
legend('Training Set Expected Values', 'Training Set Predicted Values', 'Validation Set Expected Values', 'Validation Set Predicted Values');
xlabel('Number of Data Groups');
ylabel('value');
title('');
hold off;