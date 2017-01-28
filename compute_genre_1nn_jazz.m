clc;
clear all;
close all;

% Jazz Training
fclose('all');

if exist('./G1C_Features_jazz_training.mat', 'file')
    delete('./G1C_Features_jazz_training.mat');
end

% ma_g1c_FeatureExtraction('./training_jazz.txt', './', 'G1C_Features_jazz_training.mat');

ma_g1c_FeatureExtraction('./training_5fold_jazz.txt', './', 'G1C_Features_jazz_training.mat');

load('./G1C_Features_jazz_training.mat');

for i=1:length(data.filenames)-1
    for j=i+1:length(data.filenames),
        if all([data.feat.g1c.max_ico(i),data.feat.g1c.max_ico(j)]<10^10),
                tmp = squeeze(data.feat.g1.m(i,:,:))-squeeze(data.feat.g1.m(j,:,:));
                d_g1 = ... %% kl distance
                    trace(squeeze(data.feat.g1.co(i,:,:))*squeeze(data.feat.g1.ico(j,:,:))) + ...
                    trace(squeeze(data.feat.g1.co(j,:,:))*squeeze(data.feat.g1.ico(i,:,:))) + ...
                    trace((squeeze(data.feat.g1.ico(i,:,:))+squeeze(data.feat.g1.ico(j,:,:)))*tmp'*tmp);
                d_g1_computed = true;
        end
        
        d_fp = norm(data.feat.fp(i,:)-data.feat.fp(j,:));
        d_fpg = abs(data.feat.fpg(i)-data.feat.fpg(j));
        d_fpb = abs(data.feat.fp_bass(i)-data.feat.fp_bass(j));
    end    
    
    jazz_training(i,:) = [d_fp d_fpb d_fpg d_g1 1];
end

% Non Jazz Training
fclose('all');

if exist('./G1C_Features_non_jazz_training.mat', 'file')
    delete('./G1C_Features_non_jazz_training.mat');
end

% ma_g1c_FeatureExtraction('./training_nonjazz.txt', './', 'G1C_Features_non_jazz_training.mat');

ma_g1c_FeatureExtraction('./training_5fold_nonjazz.txt', './', 'G1C_Features_non_jazz_training.mat');

load('./G1C_Features_non_jazz_training.mat');

for i=1:length(data.filenames)-1
    for j=i+1:length(data.filenames),
        if all([data.feat.g1c.max_ico(i),data.feat.g1c.max_ico(j)]<10^10),
                tmp = squeeze(data.feat.g1.m(i,:,:))-squeeze(data.feat.g1.m(j,:,:));
                d_g1 = ... %% kl distance
                    trace(squeeze(data.feat.g1.co(i,:,:))*squeeze(data.feat.g1.ico(j,:,:))) + ...
                    trace(squeeze(data.feat.g1.co(j,:,:))*squeeze(data.feat.g1.ico(i,:,:))) + ...
                    trace((squeeze(data.feat.g1.ico(i,:,:))+squeeze(data.feat.g1.ico(j,:,:)))*tmp'*tmp);
                d_g1_computed = true;
        end
        
        d_fp = norm(data.feat.fp(i,:)-data.feat.fp(j,:));
        d_fpg = abs(data.feat.fpg(i)-data.feat.fpg(j));
        d_fpb = abs(data.feat.fp_bass(i)-data.feat.fp_bass(j));        
    end    
    
    non_jazz_training(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 

end

jazz_count = 0;
non_jazz_count = 0;
false_positive = 0;

% Jazz Testing
fclose('all');

if exist('./G1C_Features_jazz_testing.mat', 'file')
    delete('./G1C_Features_jazz_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_jazz.txt', './', 'G1C_Features_jazz_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_jazz.txt', './', 'G1C_Features_jazz_testing.mat');

load('./G1C_Features_jazz_testing.mat');

for i=1:length(data.filenames)-1
    for j=i+1:length(data.filenames),
        if all([data.feat.g1c.max_ico(i),data.feat.g1c.max_ico(j)]<10^10),
                tmp = squeeze(data.feat.g1.m(i,:,:))-squeeze(data.feat.g1.m(j,:,:));
                d_g1 = ... %% kl distance
                    trace(squeeze(data.feat.g1.co(i,:,:))*squeeze(data.feat.g1.ico(j,:,:))) + ...
                    trace(squeeze(data.feat.g1.co(j,:,:))*squeeze(data.feat.g1.ico(i,:,:))) + ...
                    trace((squeeze(data.feat.g1.ico(i,:,:))+squeeze(data.feat.g1.ico(j,:,:)))*tmp'*tmp);
                d_g1_computed = true;
        end
        
        d_fp = norm(data.feat.fp(i,:)-data.feat.fp(j,:));
        d_fpg = abs(data.feat.fpg(i)-data.feat.fpg(j));
        d_fpb = abs(data.feat.fp_bass(i)-data.feat.fp_bass(j));
    end
    
    jazz_testing(i,:) = [d_fp d_fpb d_fpg d_g1 1]; 
    
    dist(i,1) = sqrt((jazz_testing(i,1) - 299.4943).^2 + ...
        (jazz_testing(i,2) - 164.7449).^2 + ...
        (jazz_testing(i,3) - 0.7772).^2 + ...
        (jazz_testing(i,4) - 105.2591).^2);

    dist(i,2) = sqrt((jazz_testing(i,1) - 336.0751).^2 + ...
        (jazz_testing(i,2) - 183.9396).^2 + ...
        (jazz_testing(i,3) - 0.6920).^2 + ...
        (jazz_testing(i,4) - 89.7987).^2);

    [min_val, idx] = min(abs(dist(i,:)));
    if idx == 1
        jazz_count = jazz_count + 1;
    else
        false_positive = false_positive + 1;
    end
end

% Non Jazz Testing
fclose('all');

if exist('./G1C_Features_non_jazz_testing.mat', 'file')
    delete('./G1C_Features_non_jazz_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_nonjazz.txt', './', 'G1C_Features_non_jazz_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_nonjazz.txt', './', 'G1C_Features_non_jazz_testing.mat');

load('./G1C_Features_non_jazz_testing.mat');

for i=1:length(data.filenames)-1
    for j=i+1:length(data.filenames),
        if all([data.feat.g1c.max_ico(i),data.feat.g1c.max_ico(j)]<10^10),
                tmp = squeeze(data.feat.g1.m(i,:,:))-squeeze(data.feat.g1.m(j,:,:));
                d_g1 = ... %% kl distance
                    trace(squeeze(data.feat.g1.co(i,:,:))*squeeze(data.feat.g1.ico(j,:,:))) + ...
                    trace(squeeze(data.feat.g1.co(j,:,:))*squeeze(data.feat.g1.ico(i,:,:))) + ...
                    trace((squeeze(data.feat.g1.ico(i,:,:))+squeeze(data.feat.g1.ico(j,:,:)))*tmp'*tmp);
                d_g1_computed = true;
        end
        
        d_fp = norm(data.feat.fp(i,:)-data.feat.fp(j,:));
        d_fpg = abs(data.feat.fpg(i)-data.feat.fpg(j));
        d_fpb = abs(data.feat.fp_bass(i)-data.feat.fp_bass(j));
    end
    
    non_jazz_testing(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 
    
    dist(i,1) = sqrt((non_jazz_testing(i,1) - 299.4943).^2 + ...
        (non_jazz_testing(i,2) - 164.7449).^2 + ...
        (non_jazz_testing(i,3) - 0.7772).^2 + ...
        (non_jazz_testing(i,4) - 105.2591).^2);

    dist(i,2) = sqrt((non_jazz_testing(i,1) - 336.0751).^2 + ...
        (non_jazz_testing(i,2) - 183.9396).^2 + ...
        (non_jazz_testing(i,3) - 0.6920).^2 + ...
        (non_jazz_testing(i,4) - 89.7987).^2);
    
    [min_val, idx] = min(abs(dist));
    if idx == 1
        false_positive = false_positive + 1;
    else
        non_jazz_count = non_jazz_count + 1;
    end
end