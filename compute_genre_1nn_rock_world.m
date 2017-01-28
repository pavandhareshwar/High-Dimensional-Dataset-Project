clc;
clear all;
close all;

% Rock Training
fclose('all');

if exist('./G1C_Features_rock_training.mat', 'file')
    delete('./G1C_Features_rock_training.mat');
end

% ma_g1c_FeatureExtraction('./training_rock.txt', './', 'G1C_Features_rock_training.mat');

ma_g1c_FeatureExtraction('./training_5fold_rock.txt', './', 'G1C_Features_rock_training.mat');

load('./G1C_Features_rock_training.mat');

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
    
    rock_training(i,:) = [d_fp d_fpb d_fpg d_g1 1];
end

% Non Rock Training
fclose('all');

if exist('./G1C_Features_world_training.mat', 'file')
    delete('./G1C_Features_world_training.mat');
end

% ma_g1c_FeatureExtraction('./training_world.txt', './', 'G1C_Features_world_training.mat');

ma_g1c_FeatureExtraction('./training_5fold_world.txt', './', 'G1C_Features_world_training.mat');

load('./G1C_Features_world_training.mat');

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
    
    world_training(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 

end

rock_count = 0;
world_count = 0;
false_positive = 0;

% Rock Testing
fclose('all');

if exist('./G1C_Features_rock_testing.mat', 'file')
    delete('./G1C_Features_rock_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_rock.txt', './', 'G1C_Features_rock_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_rock.txt', './', 'G1C_Features_rock_testing.mat');

load('./G1C_Features_rock_testing.mat');

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
    
    rock_testing(i,:) = [d_fp d_fpb d_fpg d_g1 1]; 
    
    dist(i,1) = sqrt((rock_testing(i,1) - 634.0594).^2 + ...
        (rock_testing(i,2) - 715.3273).^2 + ...
        (rock_testing(i,3) - 0.7508).^2 + ...
        (rock_testing(i,4) - 72.6104).^2);

    dist(i,2) = sqrt((rock_testing(i,1) - 319.7992).^2 + ...
        (rock_testing(i,2) - 202.8354).^2 + ...
        (rock_testing(i,3) - 0.5823).^2 + ...
        (rock_testing(i,4) - 99.2944).^2);

    [min_val, idx] = min(abs(dist(i,:)));
    if idx == 1
        rock_count = rock_count + 1;
    else
        false_positive = false_positive + 1;
    end
end

% Non Classical Testing
fclose('all');

if exist('./G1C_Features_world_testing.mat', 'file')
    delete('./G1C_Features_world_testing.mat');
end

% ma_g1c_FeatureExtraction('./testing_world.txt', './', 'G1C_Features_world_testing.mat');

ma_g1c_FeatureExtraction('./testing_5fold_world.txt', './', 'G1C_Features_world_testing.mat');

load('./G1C_Features_world_testing.mat');

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
    
    world_testing(i,:) = [d_fp d_fpb d_fpg d_g1 2]; 
    
    dist(i,1) = sqrt((world_testing(i,1) - 634.0594).^2 + ...
        (world_testing(i,2) - 715.3273).^2 + ...
        (world_testing(i,3) - 0.7508).^2 + ...
        (world_testing(i,4) - 72.6104).^2);

    dist(i,2) = sqrt((world_testing(i,1) - 319.7992).^2 + ...
        (world_testing(i,2) - 202.8354).^2 + ...
        (world_testing(i,3) - 0.5823).^2 + ...
        (world_testing(i,4) - 99.2944).^2);
    
    [min_val, idx] = min(abs(dist));
    if idx == 1
        false_positive = false_positive + 1;
    else
        world_count = world_count + 1;
    end
end