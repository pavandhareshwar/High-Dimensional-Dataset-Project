clc;
clear all;
close all;

% Metal Training
fclose('all');

if exist('./G1C_Features_metal_training.mat', 'file')
    delete('./G1C_Features_metal_training.mat');
end

% ma_g1c_FeatureExtraction('./training_metal.txt', './', 'G1C_Features_metal_training.mat');

ma_g1c_FeatureExtraction('./training_5fold_metal.txt', './', 'G1C_Features_metal_training.mat');

load('./G1C_Features_metal_training.mat');

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
    
    metal_training(i,:) = [d_fp d_fpb d_fpg d_g1 1];
end

% Non Metal Training
fclose('all');

if exist('./G1C_Features_non_metal_training.mat', 'file')
    delete('./G1C_Features_non_metal_training.mat');
end

% ma_g1c_FeatureExtraction('./training_nonmetal.txt', './', 'G1C_Features_non_metal_training.mat');

ma_g1c_FeatureExtraction('./training_5fold_nonmetal.txt', './', 'G1C_Features_non_metal_training.mat');

load('./G1C_Features_non_metal_training.mat');

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
    
    non_metal_training(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 

end

metal_count = 0;
non_metal_count = 0;
false_positive = 0;

% Metal Testing
fclose('all');

if exist('./G1C_Features_metal_testing.mat', 'file')
    delete('./G1C_Features_metal_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_metal.txt', './', 'G1C_Features_metal_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_metal.txt', './', 'G1C_Features_metal_testing.mat');

load('./G1C_Features_metal_testing.mat');

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
    
    metal_testing(i,:) = [d_fp d_fpb d_fpg d_g1 1]; 
    
    dist(i,1) = sqrt((metal_testing(i,1) - 136.4170).^2 + ...
        (metal_testing(i,2) - 178.4570).^2 + ...
        (metal_testing(i,3) - 0.7572).^2 + ...
        (metal_testing(i,4) - 61.4877).^2);

    dist(i,2) = sqrt((metal_testing(i,1) - 320.6767).^2 + ...
        (metal_testing(i,2) - 200.4132).^2 + ...
        (metal_testing(i,3) - 0.623).^2 + ...
        (metal_testing(i,4) - 88.1832).^2);

    [min_val, idx] = min(abs(dist(i,:)));
    if idx == 1
        metal_count = metal_count + 1;
    else
        false_positive = false_positive + 1;
    end
end

% Non Metal Testing
fclose('all');

if exist('./G1C_Features_non_metal_testing.mat', 'file')
    delete('./G1C_Features_non_metal_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_nonmetal.txt', './', 'G1C_Features_non_metal_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_nonmetal.txt', './', 'G1C_Features_non_metal_testing.mat');

load('./G1C_Features_non_metal_testing.mat');

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
    
    non_metal_testing(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 
    
    dist(i,1) = sqrt((non_metal_testing(i,1) - 136.4170).^2 + ...
        (non_metal_testing(i,2) - 178.4570).^2 + ...
        (non_metal_testing(i,3) - 0.7572).^2 + ...
        (non_metal_testing(i,4) - 61.4877).^2);

    dist(i,2) = sqrt((non_metal_testing(i,1) - 320.6767).^2 + ...
        (non_metal_testing(i,2) - 200.4132).^2 + ...
        (non_metal_testing(i,3) - 0.623).^2 + ...
        (non_metal_testing(i,4) - 88.1832).^2);
    
    [min_val, idx] = min(abs(dist));
    if idx == 1
        false_positive = false_positive + 1;
    else
        non_metal_count = non_metal_count + 1;
    end
end