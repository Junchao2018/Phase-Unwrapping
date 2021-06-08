%% This Matlab code is for generating simulated data (Wrapped and Unwrapped Phases).
%%
TraingingdataName      = 'TrainingData_Phase';
TestdataName      = 'TestData_Phase';
H = 400;
W = 400;
CHANNEL = 1;
batchSize = 3;
TestNum = 500;
Num = 10000;
%%
if ~exist(TraingingdataName,'file')
    mkdir(TraingingdataName);
end
if ~exist(TestdataName,'file')
    mkdir(TestdataName);
end
%%
TrainingNum = Num - TestNum;
numPatches = ceil(TrainingNum/batchSize)*batchSize;
disp([numPatches,batchSize,numPatches/batchSize]);
training_inputs  = zeros(numPatches,H, W, CHANNEL, 'single'); % this is fast
training_outputs  = zeros(numPatches,H, W, 1, 'int8'); % this is fast
training_contours  = zeros(numPatches,H, W, 1, 'int8'); % this is fast


test_inputs  = zeros(TestNum,H, W, CHANNEL, 'single'); % this is fast
test_outputs  = zeros(TestNum,H, W, 1, 'single'); % this is fast
test_contours  = zeros(TestNum,H, W, 1, 'single'); % this is fast


c_train = 0;
c_test = 0;
count = 0;

while(count<Num)
    zterm=(rand(1,45)-0.5)*5;
    [DC,coeff]=ZernikeCalc([1:45],zterm',[],'square');
    DC=DC(101:500,101:500);
    unwrapped=DC;
    I1=1+cos(unwrapped);
    I2=1+cos(unwrapped+pi/2);
    I3=1+cos(unwrapped+pi);
    I4=1+cos(unwrapped+pi*1.5);   
    %%
    phase=atan2(I4-I2,I1-I3);
    %%
    Tmp = (unwrapped - phase)/pi/2;
    Tmp = max(max(abs(Tmp)));
    if Tmp<9
        count = count + 1;
        Unwrapped(:,:,count)=unwrapped;
        Wrapped(:,:,count)=phase;
        
        
        unwrappedimg_Tra = (unwrapped - phase)/pi/2;
        ime = edge(unwrappedimg_Tra,'Canny');
        ime = int8(ime);
        if count > TrainingNum
            c_test = c_test + 1;
            test_inputs(c_test,:,:,:) = phase;
            test_outputs(c_test,:,:,:) = unwrapped;
            test_contours(c_test,:,:,:) = ime;
        else
            c_train = c_train + 1;
            training_inputs(c_train,:,:,:)   = phase;
            training_outputs(c_train,:, :, :)   = int8(unwrappedimg_Tra);
            training_contours(c_train,:, :, :)   = ime;
        end
    end 
end
TT = (Unwrapped - Wrapped)/pi;
a = max(TT(:))
b = min(TT(:))
training_outputs = training_outputs -b;

save(fullfile(TraingingdataName,'imtrain'), 'training_inputs','training_outputs','training_contours','-v7.3');
save(fullfile(TestdataName,'imtest'), 'test_inputs','test_outputs','test_contours','-v7.3');

