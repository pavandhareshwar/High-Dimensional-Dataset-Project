The song that we want to analyze and categorize into the right genre is approximately 1.3 million data points. 
Its hard to justify analyzing such a huge dataset, especially when each data-point gives very little information 
about the song and we want to work with hundreds or thousands of songs. Instead, it is more practical to 
implement dimensionalty reduction techniques to reduce the number of sample points we are working with so we can 
use the entire set of songs to train and test the classifier we are building. In order to reduce the amount of data 
we are working with, while preserving as much information as possible, we need to apply some pre-processing techniques.

The invocation of ma_g1c_FeatureExtraction function(the usage instructions is listed below), part of the 
"Elias Pampalk - MA Toolbox" computes the features necessary for us to compute the distance metrics for evaluating 
spectral similarity. The ma_g1c_FeatureExtraction function internally invokes the ma_mfcc function, which computes 
the Mel Frequency Cepstral Coefficients (MFCCs) followed by Discrete Cosine transform (DCT) on the raw data. 
DCT is used here to play the role of K-L transform, as the transformation function for DCT is a fixed one but 
for K-L transform, it has to be computed for each song and this can an overhead. The MFCCs is an effective tool 
to extract information from an audio signal. To compute MFCCs, we use some very simple filters and transformations 
that roughly model some of the characteristics of the human auditory system. The important characteristics are:

1. The non-linear frequency resolution using the Mel frequency scale.

2. The non-linear perception of loudness using decibel.

3. Spectral masking effects using DCT.

With Mel Frequency Cepstral Coefficients (MFCC), the dimensionality of each song was reduced from a million data points 
to a several hundred thousand points, but this wasn’t good enough to go ahead with computation of feature vectors for 
classification. To further reduce the dimension of the song, we used fluctuation patterns (FP) in combination
with MFCC. The ma_g1c_FeatureExtraction function is used to extract the fluctuation patterns.

Once the fluctuation patterns are extracted, the distance metrics are computed. Similarity between individual spectra, 
or between clusters of similar spectra, can be mathematically analyzed using different distance metrics between points
in the feature space. So, given a query song, the objective of classification problem is to find out the songs that 
match the query song (in terms of clustering, identify the cluster to which the query song can be assigned to) according 
to a similarity metric. The smaller the distance between two vectors in the feature space, the greater the similarity.

One thing we tried in this project was to work with just one genre at a time and try classifying the song as either 
belonging to the genre or not (ex either the song belongs to classical or not). We used compute_genre_1nn_xxx.m 
where xxx is the different genre we tried to classify (ex. compute_genre_1nn_classical.m for classical). The
compute_genre_1nn_xxx.m file uses the features extracted by ma_g1c_FeatureExtraction function and uses them to 
initially train the classifier to recognize the song genre and then during the testing phase, to categorize the 
song into the right genre.

ma_g1c_FeatureExtraction -- extract the features necessary to compute the distance metrics used by the classifier to 
train/categorize the song into right genre.

%% USAGE EXAMPLES
%% FeatureExtraction(in_file,out_dir)
%%
%FeatureExtraction('mypath/myInputFile.txt','mypath/myOutputDirectory')
%% FeatureExtraction('mypath/myInputFile','mypath/myOutputDirectory/')
%%
%% INPUT ARGUMENTS
%%
%% in_file: path fo file containing list of wav files to extract
%features from
%% all files are 22050Hz mono wav (this function checks if this is true)
%% out_dir: directory to which this function can write output.
%(logfiles,
%% and extracted features), files written to out_dir must be read
%% by function "DistanceCumputation"
%%
%% in_file format (text):
%% path/to/audiofile/0000001.wav
%% path/to/audiofile/0000002.wav
%% path/to/audiofile/0000003.wav
%% ...
%% path/to/audiofile/9999999.wav


