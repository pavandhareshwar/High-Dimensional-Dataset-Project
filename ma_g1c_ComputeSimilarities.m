function ComputeSimilarities(in_dir,out_file)
%%
%%         G1C, implementation as described in thesis
%%         for mirex, "music audio search" (was audio music similarity)
%%
%% USAGE EXAMPLE
%%       FeatureExtraction(in_file,out_dir)
%%       ComputeSimilarities(out_dir,'somepath/distance_matrix.txt')
%%
%% INPUT ARGUMENTS
%%
%% in_dir:   directory to which function "FeatureExtraction" was writing to.
%%           a log file will be created in this directory
%% out_file: whole distance matrix in the following format
%%
%%   <start file, exclude this line>
%%     <SUBMISSION_NAME>
%%     1   path\file1.wav
%%     2   path\file2.wav
%%     ...
%%     N  path\fileN.wav
%%     Q\R   1       2      ...      N
%%     1     0.000  10.234  ...   123.32
%%     2    10.234   0.000  ...    23.45
%%     .     ...     ...    0.000   ...
%%     N     4.1234  6.345  ...     0.0
%%   <end file, exclude this line>
%%
%%   delimiter: tabulator space, number format: float

%% HARDCODED PARAMETERS
exit_on_error = 0; %% set to 1 for final submission (use 0 for testing)
data.submission_name = 'G1C';
%% 

if nargin~=2,
    error('Number of input arguments is not 2. (try "help ComputeSimilarities")')
end

disp('--> G1C (Elias Pampalk, MIREX''06)')

if in_dir(end)~='/' && in_dir(end)~='\',
    in_dir(end+1)='/';
end
in_file = [in_dir,'G1C_features.mat'];

%% TEST WRITE ACCESS TO OUTPUT FILE
fid = fopen(out_file,'w');
if fid==-1, error('cannot write to output file (distance matrix)'); end
fprintf(fid,'%s','testwrite'); fclose(fid);              
delete(out_file);

%% TEST INPUT FILE
fid = fopen(in_file,'r');
if fid==-1, error('cannot find file with extracted features (path problem?)'); end
fclose(fid);
load(in_file)

%% START LOGFILE
logfile = [in_dir,'ComputeSimilarities-',data.submission_name,'-logfile.txt'];
fid = fopen(logfile,'a');
if fid==-1, error('can''t append logfile'); end
fclose(fid);              

mydisp(logfile,datestr(now));
mydisp(logfile, '-> ComputeSimilarities called.')
mydisp(logfile,['   Input directory: ',in_dir])
mydisp(logfile,'         (writing log-file to input directory)')
mydisp(logfile,['   Output file: ',out_file])

try %% big try catch to catch every error, write it to logfile and exit
    mydisp(logfile,'start computing distances ...')

    D = zeros(length(data.filenames),length(data.filenames));
    t0 = cputime;
    t1 = cputime;
    num_computations = length(data.filenames)*(length(data.filenames)-1)/2;
    num_computations_sofar = 0;
    for i=1:length(data.filenames)-1,
        t2 = cputime;
        if t2-t1>2, %% output current status only every 5 seconds (cputime)
            tmp_avg_t_sofar = (t1-t0)/(num_computations_sofar-1);
            tmp_est_tot = (t1-t0)/(num_computations_sofar-1)*num_computations;
            tmp_est_rem = tmp_est_tot - (t1-t0);
            mydisp(logfile,[num2str(num_computations_sofar),'/',num2str(num_computations), ...
                ' CS est rem ',num2str(tmp_est_rem/60),'m, est tot ', ...
                num2str(tmp_est_tot/60),'m'])
            t1 = cputime;
        end
        
        for j=i+1:length(data.filenames),
            d_g1_computed = false;
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

%             if d_g1_computed,
%                 D(i,j) = ...
%                     0.7*(-exp(-1/450*d_g1)+0.7950)/0.1493 + ...
%                     0.1*(d_fp-1688.4)/878.23 + ...
%                     0.1*(d_fpg-1.2745)/1.1245 + ...
%                     0.1*(d_fpb-1064.8)/932.79 + ...
%                     10; %% ensure all values larger than one (+1.5 should be enough in most cases)
%             else %% not evaluated work around to deal with inv covariance problems
%                 D(i,j) = ...
%                     0.4*(d_fp-1688.4)/878.23 + ...
%                     0.3*(d_fpg-1.2745)/1.1245 + ...
%                     0.3*(d_fpb-1064.8)/932.79 + ...
%                     10;
%             end

            count = 1;
            if d_g1_computed,
            % Dividing into two groups (one representing classical and other representing non-classical)
                dist(1,count) = ...
                    0.1*(d_g1-123.6654)/87.5750 + ...
                    0.4*(d_fp-354.0740)/174.7861 + ...
                    0.1*(d_fpg-0.7585)/0.5701 + ...
                    0.4*(d_fpb-151.6057)/128.9545;
                
                dist(2,count) = ...
                    0.1*(d_g1-92.4866)/43.7956 + ...
                    0.4*(d_fp-601.3099)/518.6457 + ...
                    0.1*(d_fpg-0.8668)/0.6632 + ...
                    0.4*(d_fpb-356.5221)/400.6435;
                                
%                 dist(1,count+1) = ...
%                     0.2*(d_g1-91.0186)/36.8091 + ...
%                     0.3*(d_fp-716.4740)/336.5471 + ...
%                     0.2*(d_fpg-1.2773)/0.9348 + ...
%                     0.3*(d_fpb-637.4104)/787.3581;
                
%                 dist(1,count+2) = ...
%                     1*(d_g1-134.5646)/58.3529 + ...
%                     0*(d_fp-353.2240)/127.6085 + ...
%                     0*(d_fpg-0.7816)/0.5387 + ...
%                     0*(d_fpb-321.7547)/298.1665;
%                
%                 dist(1,count+3) = ...
%                     1*(d_g1-75.0068)/16.2591 + ...
%                     0*(d_fp-255.9069)/162.5775 + ...
%                     0*(d_fpg-0.6160)/0.4511 + ...
%                     0*(d_fpb-319.7224)/359.0688;
%                 
%                 dist(1,count+4) = ...
%                     1*(d_g1-82.1003)/28.5872 + ...
%                     0*(d_fp-489.4245)/353.7383 + ...
%                     0*(d_fpg-0.7735)/0.5766 + ...
%                     0*(d_fpb-501.8877)/654.9026;
%                 
%                 dist(1,count+5) = ...
%                     1*(d_g1-138.0229)/104.2621 + ...
%                     0*(d_fp-448.0210)/213.7543 + ...
%                     0*(d_fpg-1.0480)/0.7781 + ...
%                     0*(d_fpb-480.9510)/514.6315;
            else
                dist(1,count) = ...
                    0.4*(d_fp-354.0740)/174.7861 + ...
                    0.3*(d_fpg-0.7585)/0.5701 + ...
                    0.3*(d_fpb-151.6057)/128.9545;
                
                dist(1,count+1) = ...
                    0.4*(d_fp-716.4740)/336.5471 + ...
                    0.3*(d_fpg-1.2773)/0.9348 + ...
                    0.3*(d_fpb-637.4104)/787.3581;
                
                dist(1,count+2) = ...
                    0.4*(d_fp-353.2240)/127.6085 + ...
                    0.3*(d_fpg-0.7816)/0.5387 + ...
                    0.3*(d_fpb-321.7547)/298.1665;
                
                dist(1,count+3) = ...
                    0.4*(d_fp-255.9069)/162.5775 + ...
                    0.3*(d_fpg-0.6160)/0.4511 + ...
                    0.3*(d_fpb-319.7224)/359.0688;
                
                dist(1,count+4) = ...
                    0.4*(d_fp-489.4245)/353.7383 + ...
                    0.3*(d_fpg-0.7735)/0.5766 + ...
                    0.3*(d_fpb-501.8877)/654.9026;
                
                dist(1,count+5) = ...
                    0.4*(d_fp-448.0210)/213.7543 + ...
                    0.3*(d_fpg-1.0480)/0.7781 + ...
                    0.3*(d_fpb-480.9510)/514.6315;
            end
            [min_dist, idx] = min(abs(dist));
            
            final_result(i) = idx;
            num_computations_sofar = num_computations_sofar+1;
        end     
    end
        
    D = D+D';
    disp('done. start writing output ...')

    fid = fopen(out_file,'w');
    fprintf(fid,'%s\r\n',data.submission_name);
    for i=1:length(data.filenames),
        fprintf(fid,'%d\t%s\r\n',i,data.filenames{i});
    end
    fprintf(fid,'%s','Q\R');
    fprintf(fid,'\t%d',1:length(data.filenames));
    fprintf(fid,'\r\n');
    for i=1:length(data.filenames),
        fprintf(fid,'%d',i);
        fprintf(fid,'\t%d',D(i,:));
        fprintf(fid,'\r\n');
    end
    fclose(fid);    

    disp('output file created.')

    tot_time = cputime-t0;
    disp(['total CPU time [h] ',num2str(tot_time/60/60)])
    
    mydisp(logfile,'done. exiting ...')
    mydisp(logfile,datestr(now))
    if exit_on_error, 
        exit; 
    end
catch
    mydisp(logfile,'-- caught error!')
    mydisp(logfile,lasterr)
    mydisp(logfile,datestr(now))
    if exit_on_error, 
        exit; 
    else
        error('-- error')
    end    
end

function mydisp(logfile,str)
fid = fopen(logfile,'a'); disp(str); fprintf(fid,'%s\r\n',str); fclose(fid);
