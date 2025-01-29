%% PINN ODE N degree of freedom traning
%


clc,clear all,close all,
OdeDir='ReportPINN-ODE_ND';
delete(OdeDir)
mkdir(OdeDir);
x0 = [3; 0];

N=10;
%% Scenario 1
% m=2*ones(N,1);
% k=1*ones(N,1);
% c=0.5*ones(N,1);
%% Scenario 2
sm=2;m0=2;
sk=2;k0=1;
sc=1;c0=2;
m=[m0:sm:m0+(N-1)*sm];
k=[k0:sk:k0+(N-1)*sk];
c=[c0:sc:c0+(N-1)*sc]*0.1;


M=diag(m);
K=diag(k);
C=diag(c);

M_1=inv(M);

O=zeros(N);
I=eye(N);
A=[-M_1*C -M_1*K;I O];

% % A = [-0.1 -1; 1 -0.1];% Define dynamci model
% % m=1;k=1;c=0.5;
% % C=1/m*[-c 0;0 0]*1;
% % K=1/m*[0 -k;0 0];
% % S=[0 0;1 0];%Momentum derivatives association
% % A = S+C+K+[0 -0.1];% Define dynamc model as a linear system of second order of type MSD (Mass Spring Damper)
% % 
% x0=[3;2;1];
x0=randn(2*N,1);
% % % M=randn(3)*0;
% % % A = [0 -1 0;0 0 -1;1 0 0]+(M+M')/2;% Define dynamic model
trueModel = @(t,y) A*y;
% % 
% % %% Nonlienar ode
% % trueModel = @(t,y)[-c*y(1)*y(2)-c/10*y(1)-k*y(2);y(1)];

% trueModel = @(t,y)[+k*y(2)^2-k*y(2)*atan(y(2))-c*y(1);y(1)];
% trueModel = @(t,y)[+k*y(2)^2-k*y(2)-c*y(1);y(1)];

  % trueModel = @(t,y)[k*y(2,:).^2-k*y(2,:)-c*y(1,:);y(1,:)];

numTimeSteps = 100000;
T =185;
odeOptions = odeset(RelTol=1.e-7);
t = linspace(0, T, numTimeSteps);
[~, xTrain] = ode45(trueModel, t, x0, odeOptions)
% Use simualtion data for training vector definition
xTrain = xTrain';
% Visualize the training data in a plot.

figure('Name','Phase plot of system dynamics')
tiledlayout('Flow');%Use a tailed layout
if size(xTrain,1)>1
for i=2:(N*2)
plot(xTrain(1,:),xTrain(i,:))
title("Ground Truth Dynamics") 
xlabel("x(1)") 
ylabel(["x(" num2str(i) ")"])
grid on; axis equal;
end
elseif size(xTrain,1)==3
% plot3(xTrain(1,:),xTrain(2,:),xTrain(3,:))
% title("Ground Truth Dynamics") 
% xlabel("x(1)") 
% ylabel("x(2)")
% zlabel("x(3)")
% grid onplot3(xTrain(1,:),xTrain(2,:),xTrain(3,:))
% title("Ground Truth Dynamics") 
% xlabel("x(1)") 
% ylabel("x(2)")
% zlabel("x(3)")
% grid on
end
neuralOdeTimesteps = 40;
dt = t(2);
timesteps = (0:neuralOdeTimesteps)*dt;

neuralOdeParameters = struct;

stateSize = size(xTrain,1);
hiddenSize = 20;
%% Extract neural network parametrs
neuralOdeParameters.fc1 = struct;
sz = [hiddenSize stateSize];
neuralOdeParameters.fc1.Weights = initializeGlorot(sz, hiddenSize, stateSize);
neuralOdeParameters.fc1.Bias = initializeZeros([hiddenSize 1]);

neuralOdeParameters.fc2 = struct;
sz = [stateSize hiddenSize];
neuralOdeParameters.fc2.Weights = initializeGlorot(sz, stateSize, hiddenSize);
neuralOdeParameters.fc2.Bias = initializeZeros([stateSize 1]);

% neuralOdeParameters.fc1 = struct;
% sz = [hiddenSize stateSize];
% neuralOdeParameters.fc1.Weights = initializeGlorot(sz, hiddenSize, stateSize);
% neuralOdeParameters.fc1.Bias = initializeZeros([hiddenSize 1]);
% 
% neuralOdeParameters.fc2 = struct;
% sz = [stateSize stateSize];
% neuralOdeParameters.fc2.Weights = initializeGlorot(sz, stateSize, stateSize);
% neuralOdeParameters.fc2.Bias = initializeZeros([stateSize 1]);
% 
% neuralOdeParameters.fc3 = struct;
% sz = [stateSize stateSize];
% neuralOdeParameters.fc3.Weights = initializeGlorot(sz, stateSize, stateSize);
% neuralOdeParameters.fc3.Bias = initializeZeros([stateSize 1]);
% 
% neuralOdeParameters.fc4 = struct;
% sz = [stateSize stateSize];
% neuralOdeParameters.fc4.Weights = initializeGlorot(sz, stateSize, stateSize);
% neuralOdeParameters.fc4.Bias = initializeZeros([stateSize 1]);
% 
% neuralOdeParameters.fc5 = struct;
% sz = [stateSize hiddenSize];
% neuralOdeParameters.fc5.Weights = initializeGlorot(sz, stateSize, hiddenSize);
% neuralOdeParameters.fc5.Bias = initializeZeros([stateSize 1]);


neuralOdeParameters.fc1
neuralOdeParameters.fc2

% neuralOdeParameters.fc3
% neuralOdeParameters.fc4
% neuralOdeParameters.fc5
% 
gradDecay = 0.9;
sqGradDecay = 0.999;
learnRate = 0.002;

numIter = 5000;
miniBatchSize = 200;

plotFrequency = 500;
% Initialize the averageGrad and averageSqGrad parameters for the Adam solver.
averageGrad = [];
averageSqGrad = [];
% Initialize the TrainingProgressMonitor object. Because the timer starts when you create the monitor object, make sure that you create the object close to the training loop.
monitor = trainingProgressMonitor(Metrics="Loss",Info=["Iteration","LearnRate"],XLabel="Iteration");

%% Training
numTrainingTimesteps = numTimeSteps;
trainingTimesteps = 1:numTrainingTimesteps;
plottingTimesteps = 2:numTimeSteps;

iteration = 0;

while iteration < numIter && ~monitor.Stop
    iteration = iteration + 1;

    % Create batch
    [X, targets] = createMiniBatch(numTrainingTimesteps, neuralOdeTimesteps, miniBatchSize, xTrain);

    % Evaluate network and compute loss and gradients
    [loss,gradients] = dlfeval(@modelLoss,timesteps,X,neuralOdeParameters,targets);

    % Update network
    [neuralOdeParameters,averageGrad,averageSqGrad] = adamupdate(neuralOdeParameters,gradients,averageGrad,averageSqGrad,iteration,...
        learnRate,gradDecay,sqGradDecay);

    % Plot loss
    recordMetrics(monitor,iteration,Loss=loss);

    % Plot predicted vs. real dynamics (updated graph)
    if mod(iteration,plotFrequency) == 0  || iteration == 1
        % Use ode45 to compute the solution 
        %% 
        y = dlode45(@odeModel,t,dlarray(x0),neuralOdeParameters,DataFormat="CB");
        %nexttile; 
       for i=1:N*2 
        subplot(round(N),2,i)
        plot(xTrain(1,plottingTimesteps),xTrain(i,plottingTimesteps),"r--")
        hold on
        plot(y(1,:),y(i,:),"b-")
        hold off
        xlabel("x(1)")
        ylabel(['x(' num2str(i) ')'])
        title(["Predicted vs. Real Dynamics: iteration " num2str(iteration)])
        legend("Training Ground Truth", "Predicted");
        grid on;axis auto;
        drawnow
       end
     saveas(gcf,[OdeDir '\PahsePlot#'  num2str(iteration)  '.fig' ]);
    end
    updateInfo(monitor,Iteration=iteration,LearnRate=learnRate);
    monitor.Progress = 100*iteration/numIter;
end
%% Post training validation
N_val=50;
ScaleDeviation=ones(N_val)*0.1;
for h=1:N_val
tPred = t;
% x0Pred1 = sqrt([2;2]);
% x0Pred2 = [-1;-1.5];
% x0Pred3 = [0;2];
% x0Pred4 = [-2;0];

x0Pred1 = x0+(ScaleDeviation(h)*rand(2*N,1));
% x0Pred2 = x0+(ScaleDeviation(h)*randn(2*N,1));
% x0Pred3 = x0+(ScaleDeviation(h)*randn(2*N,1));
% x0Pred4 = x0+(ScaleDeviation(h)*randn(2*N,1));

[~, xTrue1] = ode45(trueModel, tPred, x0Pred1, odeOptions);
% [~, xTrue2] = ode45(trueModel, tPred, x0Pred2, odeOptions);
% [~, xTrue3] = ode45(trueModel, tPred, x0Pred3, odeOptions);
% [~, xTrue4] = ode45(trueModel, tPred, x0Pred4, odeOptions);

xPred1 = dlode45(@odeModel,tPred,dlarray(x0Pred1),neuralOdeParameters,DataFormat="CB");
% xPred2 = dlode45(@odeModel,tPred,dlarray(x0Pred2),neuralOdeParameters,DataFormat="CB");
% xPred3 = dlode45(@odeModel,tPred,dlarray(x0Pred3),neuralOdeParameters,DataFormat="CB");
% xPred4 = dlode45(@odeModel,tPred,dlarray(x0Pred4),neuralOdeParameters,DataFormat="CB");

%% Plotting
figure('Name',['Predicted vs expected:  ' num2str(h)])
% subplot(2,2,1)
% plotTrueAndPredictedSolutions(xTrue1, xPred1);grid on;
% subplot(2,2,2)
% plotTrueAndPredictedSolutions(xTrue2, xPred2);grid on;
% subplot(2,2,3)
% plotTrueAndPredictedSolutions(xTrue3, xPred3);grid on;
% subplot(2,2,4)
% plotTrueAndPredictedSolutions(xTrue4, xPred4);grid on;
%%
plotTrueAndPredictedSolutionsN(xTrue1, xPred1);grid on;
end 
%% Helper functions
function X = model(tspan,X0,neuralOdeParameters)
X = dlode45(@odeModel,tspan,X0,neuralOdeParameters,DataFormat="CB");
end

function y = odeModel(~,y,theta)
y = tanh(theta.fc1.Weights*y + theta.fc1.Bias);
y = theta.fc2.Weights*y + theta.fc2.Bias;
end
%% Model loss function
function [loss,gradients] = modelLoss(tspan,X0,neuralOdeParameters,targets)

% Compute predictions.
X = model(tspan,X0,neuralOdeParameters);

% Compute L1 loss.
loss = l1loss(X,targets,NormalizationFactor="all-elements",DataFormat="CBT");

% Compute gradients.
gradients = dlgradient(loss,neuralOdeParameters);

end
%% Batch function
function [x0, targets] = createMiniBatch(numTimesteps,numTimesPerObs,miniBatchSize,X)

% Create batches of trajectories.
s = randperm(numTimesteps - numTimesPerObs, miniBatchSize);

x0 = dlarray(X(:, s));
targets = zeros([size(X,1) miniBatchSize numTimesPerObs]);

for i = 1:miniBatchSize
    targets(:, i, 1:numTimesPerObs) = X(:, s(i) + 1:(s(i) + numTimesPerObs));
end

end
%%
function plotTrueAndPredictedSolutions(xTrue,xPred)

xPred = squeeze(xPred)';

err = mean(abs(xTrue(2:end,:) - xPred), "all");

plot(xTrue(:,1),xTrue(:,2),"r--",xPred(:,1),xPred(:,2),"b-",LineWidth=1)

title("Absolute Error = " + num2str(err,"%.4f"))
xlabel("x(1)")
ylabel("x(2)")

xlim([-2 3])
ylim([-2 3])

legend("Ground Truth","Predicted")

end
%%
function plotTrueAndPredictedSolutionsN(xTrue,xPred)
N=size(xTrue,2);
xPred = squeeze(xPred)';
err = mean(abs(xTrue(2:end,:) - xPred), "all");
tiledlayout('flow');
for k=1:N
nexttile;
plot(xTrue(:,1),xTrue(:,k),"r--",xPred(:,1),xPred(:,k),"b-",LineWidth=1);grid on;
title("Absolute Error = " + num2str(err,"%.4f"))
xlabel("x(1)")
ylabel("x(2)")

xlim([-2 3])
ylim([-2 3])

legend("Ground Truth","Predicted")
end
end
%%
function parameter = initializeZeros(sz,className)

arguments
    sz
    className = 'single'
end

parameter = zeros(sz,className);
parameter = dlarray(parameter);

end

function weights = initializeGlorot(sz,numOut,numIn,className)

arguments
    sz
    numOut
    numIn
    className = 'single'
end

Z = 2*rand(sz,className) - 1;
bound = sqrt(6 / (numIn + numOut));

weights = bound * Z;
weights = dlarray(weights);

end