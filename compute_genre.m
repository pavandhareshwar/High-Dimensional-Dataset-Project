clc;
clear all;
close all;

false_positives = 0;

% fclose('all');
% 
% if exist('./G1C_Features.mat', 'file')
%     delete('./G1C_Features.mat');
% end
% 
% ma_g1c_FeatureExtraction('./tracks.txt', './', 'G1C_Features.mat');

load('./G1C_Features.mat');

% Dividing the total 729 songs into classical and non-classical songs

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

    count = 1;
    
    dist(1,count) = 0.05*(d_g1-123.6654) + 0.45*(d_fp-354.0740) + ...
        0.05*(d_fpg-0.7585) + 0.45*(d_fpb-151.6057);
                
    dist(2,count) = 0.05*(d_g1-92.4866) + 0.45*(d_fp-601.3099) + ...
        0.05*(d_fpg-0.8668) + 0.45*(d_fpb-356.5221);

    [min_dist, idx] = min(abs(dist));
    
    final_result_classical(i) = idx;
end

class_count = 0;
non_class_count = 0;
non_matched_class_count = 0;

for index=1:length(data.filenames)-1
    if index >= 1 && index <= 320
        if final_result_classical(index) == 1
            class_count = class_count + 1;
            class_songs(class_count,1) = data.filenames(index);
        else
%             non_matched_class_count = non_matched_class_count + 1;
%             non_matched_class_songs(non_matched_class_count,1) = data.filenames(index);
            non_class_count = non_class_count + 1;
            non_class_songs(non_class_count,1) = data.filenames(index);
            false_positives = false_positives + 1;
        end
    else
        if final_result_classical(index) == 2
            non_class_count = non_class_count + 1;
            non_class_songs(non_class_count,1) = data.filenames(index);
        else
%             non_matched_class_count = non_matched_class_count + 1;
%             non_matched_class_songs(non_matched_class_count,1) = data.filenames(index);
            non_class_count = non_class_count + 1;
            non_class_songs(non_class_count,1) = data.filenames(index);
            false_positives = false_positives + 1;
        end
    end
end

fclose('all');

if exist('./non_classical_songs.txt', 'file')
    delete('./non_classical_songs.txt');
end

fileID = fopen('./non_classical_songs.txt', 'wt');
for i=1:length(non_class_songs)
    fprintf(fileID, '%s\n', non_class_songs{i});
end
fclose(fileID);

% if exist('./G1C_Features_non_classical.mat', 'file')
%     delete('./G1C_Features_non_classical.mat');
% end
% 
% ma_g1c_FeatureExtraction('./non_classical_songs.txt', './', 'G1C_Features_non_classical.mat');

load('./G1C_Features_non_classical.mat');

% Non Classical (Split into Jazz and Non-Jazz)

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

    count = 1;

    dist(1,count) = 0.33*(d_g1-134.5646) + 0.33*(d_fp-353.2240) + ...
        0.01*(d_fpg-0.7816) + 0.33*(d_fpb-321.7547);
                
    dist(2,count) = 0.33*(d_g1-77.2302) + 0.33*(d_fp-687.4780) + ...
        0.01*(d_fpg-0.8481) + 0.33*(d_fpb-374.7253);

    [min_dist, idx] = min(abs(dist));
    
    final_result_jazz(i) = idx;
end

jazz_count = 0;
non_class_jazz_count = 0;
non_matched_jazz_count = 0;

for index=1:length(data.filenames)-1
    if index >= 117 && index <= 142
        if final_result_jazz(index) == 1
            jazz_count = jazz_count + 1;
            jazz_songs(jazz_count,1) = data.filenames(index);
        else
%             non_matched_jazz_count = non_matched_jazz_count + 1;
%             non_matched_jazz_songs(non_matched_jazz_count,1) = data.filenames(index);
            non_class_jazz_count = non_class_jazz_count + 1;
            non_class_jazz_songs(non_class_jazz_count,1) = data.filenames(index);
            false_positives = false_positives + 1;
        end
    else
        if final_result_jazz(index) == 2
            non_class_jazz_count = non_class_jazz_count + 1;
            non_class_jazz_songs(non_class_jazz_count,1) = data.filenames(index);
        else
%             non_matched_jazz_count = non_matched_jazz_count + 1;
%             non_matched_jazz_songs(non_matched_jazz_count,1) = data.filenames(index);
            non_class_jazz_count = non_class_jazz_count + 1;
            non_class_jazz_songs(non_class_jazz_count,1) = data.filenames(index);
            false_positives = false_positives + 1;
        end
    end
end

fclose('all');

if exist('./non_class_jazz_songs.txt', 'file')
    delete('./non_class_jazz_songs.txt');
end

fileID = fopen('./non_class_jazz_songs.txt', 'wt');
for i=1:length(non_class_jazz_songs)
    fprintf(fileID, '%s\n', non_class_jazz_songs{i});
end
fclose(fileID);

% if exist('./G1C_Features_non_class_jazz.mat', 'file')
%     delete('./G1C_Features_non_class_jazz.mat');
% end
% 
% ma_g1c_FeatureExtraction('./non_class_jazz_songs.txt', './', 'G1C_Features_non_class_jazz.mat');

load('./G1C_Features_non_class_jazz.mat');

% Non Classical Non Jazz (Split into Metal and Non-Metal)

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

    count = 1;

    dist(1,count) = 0.25*(d_g1-75.0068) + 0.25*(d_fp-255.9069) + ...
        0.25*(d_fpg-0.6160) + 0.25*(d_fpb-319.7224);
                
    dist(2,count) = 0.25*(d_g1-91.3266) + 0.25*(d_fp-647.9263) + ...
        0.25*(d_fpg-0.8596) + 0.25*(d_fpb-386.2445);

    [min_dist, idx] = min(abs(dist));
    
    final_result_metal(i) = idx;
end

metal_count = 0;
non_class_jazz_metal_count = 0;
non_matched_metal_count = 0;

for index=1:length(data.filenames)-1
    if index >= 118 && index <= 162
        if final_result_metal(index) == 1
            metal_count = metal_count + 1;
            metal_songs(metal_count,1) = data.filenames(index);
        else
%             non_matched_metal_count = non_matched_metal_count + 1;
%             non_matched_metal_songs(non_matched_metal_count,1) = data.filenames(index);
            non_class_jazz_metal_count = non_class_jazz_metal_count + 1;
            non_class_jazz_metal_songs(non_class_jazz_metal_count,1) = data.filenames(index);
            false_positives = false_positives + 1;
        end
    else
        if final_result_metal(index) == 2
            non_class_jazz_metal_count = non_class_jazz_metal_count + 1;
            non_class_jazz_metal_songs(non_class_jazz_metal_count,1) = data.filenames(index);
        else
%             non_matched_metal_count = non_matched_metal_count + 1;
%             non_matched_metal_songs(non_matched_metal_count,1) = data.filenames(index);
            non_class_jazz_metal_count = non_class_jazz_metal_count + 1;
            non_class_jazz_metal_songs(non_class_jazz_metal_count,1) = data.filenames(index);
            false_positives = false_positives + 1;
        end
    end
end

fclose('all');

if exist('./non_class_jazz_metal_songs.txt', 'file')
    delete('./non_class_jazz_metal_songs.txt');
end

fileID = fopen('./non_class_jazz_metal_songs.txt', 'wt');
for i=1:length(non_class_jazz_metal_songs)
    fprintf(fileID, '%s\n', non_class_jazz_metal_songs{i});
end
fclose(fileID);