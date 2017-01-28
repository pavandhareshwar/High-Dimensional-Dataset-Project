clc;
clear all;
close all;

% Electronic Training
% fclose('all');
% 
% if exist('./G1C_Features_electronic_training.mat', 'file')
%     delete('./G1C_Features_electronic_training.mat');
% end
% 
% % ma_g1c_FeatureExtraction('./training_electronic.txt', './', 'G1C_Features_electronic_training.mat');
% 
% ma_g1c_FeatureExtraction('./training_5fold_electronic.txt', './', 'G1C_Features_electronic_training.mat');

load('./G1C_Features_electronic_training.mat');

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
    
    electronic_training(i,:) = [d_fp d_fpb d_fpg d_g1 1];
end

% Non Electronic Training
% fclose('all');
% 
% if exist('./G1C_Features_non_electronic_training.mat', 'file')
%     delete('./G1C_Features_non_electronic_training.mat');
% end
% 
% % ma_g1c_FeatureExtraction('./training_nonelectronic.txt', './', 'G1C_Features_non_electronic_training.mat');
% 
% ma_g1c_FeatureExtraction('./training_5fold_nonelectronic.txt', './', 'G1C_Features_non_electronic_training.mat');

load('./G1C_Features_non_electronic_training.mat');

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
    
    non_electronic_training(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 

end

electronic_count = 0;
non_electronic_count = 0;
false_positive = 0;

% Electronic Testing
fclose('all');

if exist('./G1C_Features_electronic_testing.mat', 'file')
    delete('./G1C_Features_electronic_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_electronic.txt', './', 'G1C_Features_electronic_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_electronic.txt', './', 'G1C_Features_electronic_testing.mat');

load('./G1C_Features_electronic_testing.mat');

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
    
    electronic_testing(i,:) = [d_fp d_fpb d_fpg d_g1 1]; 
    
    dist(i,1) = sqrt((electronic_testing(i,1) - 569.94).^2 + ...
        (electronic_testing(i,2) - 380.64).^2 + ...
        (electronic_testing(i,3) - 0.6874).^2 + ...
        (electronic_testing(i,4) - 88.1978).^2);

    dist(i,2) = sqrt((electronic_testing(i,1) - 122.27).^2 + ...
        (electronic_testing(i,2) - 58.704).^2 + ...
        (electronic_testing(i,3) - 0.2862).^2 + ...
        (electronic_testing(i,4) - 34.086).^2);

    [min_val, idx] = min(abs(dist(i,:)));
    if idx == 1
        electronic_count = electronic_count + 1;
    else
        false_positive = false_positive + 1;
    end
end

% Non Electronic Testing
fclose('all');

if exist('./G1C_Features_non_electronic_testing.mat', 'file')
    delete('./G1C_Features_non_electronic_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_nonelectronic.txt', './', 'G1C_Features_non_electronic_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_nonelectronic.txt', './', 'G1C_Features_non_electronic_testing.mat');

load('./G1C_Features_non_electronic_testing.mat');

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
    
    non_electronic_testing(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 
    
    dist(i,1) = sqrt((non_electronic_testing(i,1) - 569.94).^2 + ...
        (non_electronic_testing(i,2) - 380.64).^2 + ...
        (non_electronic_testing(i,3) - 0.6874).^2 + ...
        (non_electronic_testing(i,4) - 88.1978).^2);

    dist(i,2) = sqrt((non_electronic_testing(i,1) - 122.27).^2 + ...
        (non_electronic_testing(i,2) - 58.704).^2 + ...
        (non_electronic_testing(i,3) - 0.2862).^2 + ...
        (non_electronic_testing(i,4) - 34.086).^2);
    
    [min_val, idx] = min(abs(dist));
    if idx == 1
        false_positive = false_positive + 1;
    else
        non_electronic_count = non_electronic_count + 1;
    end
end