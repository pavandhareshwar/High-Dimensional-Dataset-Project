clc;
clear all;
close all;

% Classical Training
fclose('all');

if exist('./G1C_Features_classical_training.mat', 'file')
    delete('./G1C_Features_classical_training.mat');
end

% ma_g1c_FeatureExtraction('./training_classical.txt', './', 'G1C_Features_classical_training.mat');

ma_g1c_FeatureExtraction('./training_5fold_classical.txt', './', 'G1C_Features_classical_training.mat');

load('./G1C_Features_classical_training.mat');

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
    
    classical_training(i,:) = [d_fp d_fpb d_fpg d_g1 1];
end

% Non Classical Training
fclose('all');

if exist('./G1C_Features_non_classical_training.mat', 'file')
    delete('./G1C_Features_non_classical_training.mat');
end

% ma_g1c_FeatureExtraction('./training_nonclassical.txt', './', 'G1C_Features_non_classical_training.mat');

ma_g1c_FeatureExtraction('./training_5fold_nonclassical.txt', './', 'G1C_Features_non_classical_training.mat');

load('./G1C_Features_non_classical_training.mat');

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
    
    non_classical_training(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 

end

classical_count = 0;
non_classical_count = 0;
false_positive = 0;

% Classical Testing
fclose('all');

if exist('./G1C_Features_classical_testing.mat', 'file')
    delete('./G1C_Features_classical_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_classical.txt', './', 'G1C_Features_classical_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_classical.txt', './', 'G1C_Features_classical_testing.mat');

load('./G1C_Features_classical_testing.mat');

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
    
    classical_testing(i,:) = [d_fp d_fpb d_fpg d_g1 1]; 
    
%     dist(i,1) = sqrt((classical_testing(i,1) - classical_training(2,1)).^2 + ...
%         (classical_testing(i,2) - classical_training(2,2)).^2 + ...
%         (classical_testing(i,3) - classical_training(2,3)).^2); % + ...
% %         (classical_testing(i,4) - classical_training(2,4)).^2);
% 
%     dist(i,2) = sqrt((classical_testing(i,1) - non_classical_training(2,1)).^2 + ...
%         (classical_testing(i,2) - non_classical_training(2,2)).^2 + ...
%         (classical_testing(i,3) - non_classical_training(2,3)).^2); % + ...
% %         (classical_testing(i,4) - non_classical_training(2,4)).^2);

    dist(i,1) = sqrt((classical_testing(i,1) - 239.14).^2 + ...
        (classical_testing(i,2) - 116.62).^2 + ...
        (classical_testing(i,3) - 0.505).^2); %+ ...
        % (classical_testing(i,4) - 749.47).^2);

    dist(i,2) = sqrt((classical_testing(i,1) - 311.035).^2 + ...
        (classical_testing(i,2) - 181.37).^2 + ...
        (classical_testing(i,3) - 0.578).^2); %+ ...
%         (classical_testing(i,4) - 65).^2);

    [min_val, idx] = min(abs(dist(i,:)));
    if idx == 1
        classical_count = classical_count + 1;
    else
        false_positive = false_positive + 1;
    end
end

% Non Classical Testing
fclose('all');

if exist('./G1C_Features_non_classical_testing.mat', 'file')
    delete('./G1C_Features_non_classical_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_nonclassical.txt', './', 'G1C_Features_non_classical_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_nonclassical.txt', './', 'G1C_Features_non_classical_testing.mat');

load('./G1C_Features_non_classical_testing.mat');

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
    
    non_classical_testing(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 
    
%     dist(i,1) = sqrt((non_classical_testing(i,1) - classical_training(2,1))^2 + ...
%         (non_classical_testing(i,2) - classical_training(2,2))^2 + ...
%         (non_classical_testing(i,3) - classical_training(2,3))^2); %+ ...
% %         (non_classical_testing(i,4) - classical_training(2,4))^2);
% 
%     dist(i,2) = sqrt((non_classical_testing(i,1) - non_classical_training(2,1))^2 + ...
%         (non_classical_testing(i,2) - non_classical_training(2,2))^2 + ...
%         (non_classical_testing(i,3) - non_classical_training(2,3))^2); % + ...
% %         (non_classical_testing(i,4) - non_classical_training(2,4))^2);

    dist(i,1) = sqrt((non_classical_testing(i,1) - 239.14).^2 + ...
        (non_classical_testing(i,2) - 116.62).^2 + ...
        (non_classical_testing(i,3) - 0.505).^2); % + ...
%         (non_classical_testing(i,4) - 749.47).^2);

    dist(i,2) = sqrt((non_classical_testing(i,1) - 311.035).^2 + ...
        (non_classical_testing(i,2) - 181.37).^2 + ...
        (non_classical_testing(i,3) - 0.578).^2); %+ ...
%         (non_classical_testing(i,4) - 65).^2);
    
    [min_val, idx] = min(abs(dist));
    if idx == 1
        false_positive = false_positive + 1;
    else
        non_classical_count = non_classical_count + 1;
    end
end