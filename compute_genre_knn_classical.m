clc;
clear all;
close all;

k = 3;

% Classical Training
fclose('all');

if exist('./G1C_Features_classical_training.mat', 'file')
    delete('./G1C_Features_classical_training.mat');
end

ma_g1c_FeatureExtraction('./training_classical.txt', './', 'G1C_Features_classical_training.mat');

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

ma_g1c_FeatureExtraction('./training_nonclassical.txt', './', 'G1C_Features_non_classical_training.mat');

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

ma_g1c_FeatureExtraction('./testing_classical.txt', './', 'G1C_Features_classical_testing.mat');

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
    
    for index = 1:min(size(classical_training,1),size(non_classical_training,1))
        dist(index,1) = sqrt((classical_testing(i,1) - classical_training(index,1)).^2 + ...
            (classical_testing(i,2) - classical_training(index,2)).^2 + ...
            (classical_testing(i,3) - classical_training(index,3)).^2 + ...
            (classical_testing(i,4) - classical_training(index,4)).^2);

        dist(index,2) = sqrt((classical_testing(i,1) - non_classical_training(index,1)).^2 + ...
            (classical_testing(i,2) - non_classical_training(index,2)).^2 + ...
            (classical_testing(i,3) - non_classical_training(index,3)).^2 + ...
            (classical_testing(i,4) - non_classical_training(index,4)).^2);
    end

    l = 1;
    for idx = 1:size(dist,1)
        distance(l) = dist(idx,1);
        distance(l+1) = dist(idx,2);
        l = l+2;
    end
    
    [dist_sorted, dist_idx] = sort(distance);
    smallestKentries = dist_sorted(1:k);
    smallestKindices = dist_idx(1:k);

    classical_nearest_match = 0;
    non_classical_nearest_match = 0;
    
    for p = 1:length(smallestKindices)
        if mod(smallestKindices(p),2) ~= 0
            classical_nearest_match = classical_nearest_match + 1;
        else
            non_classical_nearest_match = non_classical_nearest_match + 1;
        end
    end
        
    if (classical_nearest_match >= non_classical_nearest_match)
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

ma_g1c_FeatureExtraction('./testing_nonclassical.txt', './', 'G1C_Features_non_classical_testing.mat');

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

    for index = 1:min(size(classical_training,1),size(non_classical_training,1))
        dist(index,1) = sqrt((non_classical_testing(i,1) - non_classical_training(index,1)).^2 + ...
            (non_classical_testing(i,2) - non_classical_training(index,2)).^2 + ...
            (non_classical_testing(i,3) - non_classical_training(index,3)).^2 + ...
            (non_classical_testing(i,4) - non_classical_training(index,4)).^2);

        dist(index,2) = sqrt((non_classical_testing(i,1) - classical_training(index,1)).^2 + ...
            (non_classical_testing(i,2) - classical_training(index,2)).^2 + ...
            (non_classical_testing(i,3) - classical_training(index,3)).^2 + ...
            (non_classical_testing(i,4) - classical_training(index,4)).^2);
    end

    l = 1;
    for idx = 1:size(dist,1)
        distance(l) = dist(idx,1);
        distance(l+1) = dist(idx,2);
        l = l+2;
    end
    
    [dist_sorted, dist_idx] = sort(distance);
    smallestKentries = dist_sorted(1:k);
    smallestKindices = dist_idx(1:k);
    
    classical_nearest_match = 0;
    non_classical_nearest_match = 0;

    for p = 1:length(smallestKindices)
        if mod(smallestKindices(p),2) ~= 0
            non_classical_nearest_match = non_classical_nearest_match + 1;
        else
            classical_nearest_match = classical_nearest_match + 1;
        end
    end
        
    if (non_classical_nearest_match >= classical_nearest_match)
        non_classical_count = non_classical_count + 1;
    else
        false_positive = false_positive + 1;
    end
end