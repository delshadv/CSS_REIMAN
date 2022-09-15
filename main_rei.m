
%restoredefaultpath
bwd = pwd;
wd  = fullfile(bwd,'MKL');
addpath(wd)

% Setup OSL
addpath(fullfile(bwd,'osl','osl-core'))
osl_startup
osl_check_installation

% BIDS and Processed directories
bidspth = fullfile(bwd,'BioFIND','MCIControls'); %BIDS Path
procpth = fullfile(bidspth,'derivatives','meg_derivatives'); % If want maxfiltered files

% Define participant variables
participants = spm_load(fullfile(wd,'participants-imputed.tsv'));
mri_num      = grp2idx(participants.sImaging);
% Remove noisy-MRIs and non-MRI subjects
mri_num([23 197]) = 2;
nsub = 324;

% Freq bands and modalities
freqbands = {[8 12],[30 86],[2 86]};
modal = {'MEGPLANAR','MEGMAG'};

%% Rei Covs without alignment

processed_pth = fullfile(bwd,'Processed_norm_new');


for k=1%:numel(modal)
    
    covariance = cell(1,numel(freqbands));
    variance = cell(1,numel(freqbands));   
    
    for ii = 1:length(freqbands)
        
        Cov = cell(1,nsub); Var = Cov;
        cm=[];
        parfor sub = 1:nsub
            try
                 
            infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'ffdspmeeg');
            D = spm_eeg_load(infile);
            
            % Remove bad badtrials
            chans = D.indchantype(modal{k},'GOOD'); % Retrieve all good channels
            g = D(chans,:);
            g = g(:,good_samples(D,chans)); % Keep only the good samples
            
            % Filter to desired freq band
            y1 = ft_preproc_bandpassfilter(g, D.fsample, freqbands{ii}, 5, 'but');
            
            % Despiking
            y = filloutliers(y1','clip','median','ThresholdFactor',3);
            
            % Remove 50th and above PCs
            %[coeff_1,score_1,~,~,explained_1,mu1] = pca(y','NumComponents',60);
            %idx = find(cumsum(explained_1)>95,1)
            %y_1 = score_1*coeff_1' + repmat(mu1,size(chans,1),1);

            
            % Calculate Covariance Matrix
            cm(:,:,sub) = cov(y);
            

            catch
            end
            
        end

        covariance_rei_na{ii} = cm(:,:,mri_num==1);
    end
    save(sprintf('%s',modal{k}), 'covariance_rei_na');
    
    
end

%% Rei Covs with alighnment

processed_pth = fullfile(bwd,'Processed_trans');

for k=1%:numel(modal)
    
    
    for ii = 1:length(freqbands)
        
        Cov = cell(1,nsub); Var = Cov;
        cm=[];
        parfor sub = 1:nsub
            try
                 
            infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'ffdspmeeg');
            D = spm_eeg_load(infile);
            
            % Remove bad badtrials
            chans = D.indchantype(modal{k},'GOOD'); % Retrieve all good channels
            g = D(chans,:);
            g = g(:,good_samples(D,chans)); % Keep only the good samples
            
            % Filter to desired freq band
            y1 = ft_preproc_bandpassfilter(g, D.fsample, freqbands{ii}, 5, 'but');
            
            % Despiking
            y = filloutliers(y1','clip','median','ThresholdFactor',3);
            
            % Remove 50th and above PCs
            %[coeff_1,score_1,~,~,explained_1,mu1] = pca(y','NumComponents',60);
            %idx = find(cumsum(explained_1)>95,1)
            %y_1 = score_1*coeff_1' + repmat(mu1,size(chans,1),1);

            
            % Calculate Covariance Matrix
            cm(:,:,sub) = cov(y);
            

            catch
            end
            
        end

        covariance_rei_wa{ii} = cm(:,:,mri_num==1);
    end
    save(sprintf('%s',modal{k}), 'covariance_rei_wa');
    
    
end

%% Sensor Covs without alignment
processed_pth = fullfile(bwd,'Processed_norm_new');

anat = {'MRI'};
ana=1;
Nanat = length(anat);

for k=1%:numel(modal)
    
    covariance = cell(1,numel(freqbands));
    variance = cell(1,numel(freqbands));
    
    for ii = 1:length(freqbands)
        
        Cov = cell(1,nsub); Var = Cov;
        cm=[];
        parfor sub = 1:nsub
            T1file = fullfile('/imaging/henson/users/dv01/defacing_project',anat{ana},sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
            
            if exist(T1file,'file')
                
                
                infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'ffdspmeeg');
                D = spm_eeg_load(infile);
                
                % Remove bad badtrials
                chans = D.indchantype(modal{k},'GOOD'); % Retrieve all good channels
                g = D(chans,:);
                g = g(:,good_samples(D,chans)); % Keep only the good samples
                
                % Filter to desired freq band
                y1 = ft_preproc_bandpassfilter(g, D.fsample, freqbands{ii}, 5, 'but');
                
                % Despiking
                y = filloutliers(y1','clip','median','ThresholdFactor',3);
                
                % Remove 50th and above PCs
                %[coeff_1,score_1,~,~,explained_1,mu1] = pca(y','NumComponents',60);
                %idx = find(cumsum(explained_1)>95,1)
                %y_1 = score_1*coeff_1' + repmat(mu1,size(chans,1),1);
                
                
                % Calculate Covariance Matrix
                cm = cov(y);
                Cov{sub} = cm(find(triu(cm,0)))';
                
                
            end
            
        end
        Cov = cat(1,Cov{:});
        
        covariance_sen_na{ii} = Cov(:,:);
    end
    save(sprintf('%s',modal{k}), 'covariance_sen_na');
    
    
end

%% Sensor Covs with alignment

processed_pth = fullfile(bwd,'Processed_trans');

for k=1%:numel(modal)
    
    covariance = cell(1,numel(freqbands));
    variance = cell(1,numel(freqbands));
    
    for ii = 1:length(freqbands)
        
        Cov = cell(1,nsub); Var = Cov;
        cm=[];
        parfor sub = 1:nsub
            T1file = fullfile('/imaging/henson/users/dv01/defacing_project',anat{ana},sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
            
            if exist(T1file,'file')
                
                
                infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'ffdspmeeg');
                D = spm_eeg_load(infile);
                
                % Remove bad badtrials
                chans = D.indchantype(modal{k},'GOOD'); % Retrieve all good channels
                g = D(chans,:);
                g = g(:,good_samples(D,chans)); % Keep only the good samples
                
                % Filter to desired freq band
                y1 = ft_preproc_bandpassfilter(g, D.fsample, freqbands{ii}, 5, 'but');
                
                % Despiking
                y = filloutliers(y1','clip','median','ThresholdFactor',3);
                
                % Remove 50th and above PCs
                %[coeff_1,score_1,~,~,explained_1,mu1] = pca(y','NumComponents',60);
                %idx = find(cumsum(explained_1)>95,1)
                %y_1 = score_1*coeff_1' + repmat(mu1,size(chans,1),1);
                
                
                % Calculate Covariance Matrix
                cm = cov(y);
                Cov{sub} = cm(find(triu(cm,0)))';
                
                
            end
            
        end
        
        Cov = cat(1,Cov{:});
        covariance_sen_wa{ii} = Cov(:,:);
    end
    save(sprintf('%s',modal{k}), 'covariance_sen_wa');
    
    
end

%% source Cov without PCA

for k=1;%:numel(modal)
    
    covariance = cell(1,numel(freqbands));
    
    for ii = 1:length(freqbands)
        
        Cov = cell(1,nsub); Var = Cov;
        
        parfor sub = 1:nsub
            
            T1file = fullfile('/imaging/henson/users/dv01/defacing_project',anat{ana},sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
            
            if exist(T1file,'file')
                infile = fullfile('/imaging/henson/users/dv01/Github/Processed_trans',(sprintf('sub-Sub%04d',sub)),'_ffdspmeeg_187MEGPLANAR');
                D = spm_eeg_load(infile);
                D = D.montage('switch',2);
                
                % Remove bad badtrials
                chans = D.indchantype('LFP','GOOD'); % MEG :MAG or MEGPLANAR : GRD

                g = D(:,:);
                g = g(chans,good_samples(D,chans));

                
                % Filter to desired freq band
                y1 = ft_preproc_bandpassfilter(g, D.fsample, freqbands{ii}, 4, 'but');
                
                % Despiking
                y = filloutliers(y1,'clip','median','ThresholdFactor',3);
                
                % Calculate Covariance Matrix
                cm = cov(y');
                
                Cov{sub} = cm(find(triu(cm,0)))';
                
            end
        end
        Cov = cat(1,Cov{:});
        covariance_sou_nsvd{ii} = Cov;
    end
    
    save(sprintf('%s',modal{k}), 'covariance_sou_nsvd');
    
end

%% Machine Learning
V = {{covariance_sen_na{1}},{covariance_sen_wa{1}},{covariance_rei_na{1}},{covariance_rei_wa{1}},{covariance_sou_nsvd{1}}}
labels = csvread('derived/labels.csv');

rng('default') % For reproducibility
[acc4,~] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'Nfold',5,'Nrun',100,'PCA_cut',0,'feat_norm',1,'ens',0);


titles = {'Sensor','Sensor with align','Reiman','Reiman with align','Source'};

pos_titles = {'Sensor with align > Sensor','Reiman with align > Reiman','Sensor with align > Source',...
    'Reiman with align > Source'};

        
%  define contrasts
c = [-1 1 0  0 0 ; 
     0 0 -1  1 0 ; 
     0 1 0  0 -1 ;
    
     0 0 0  1 -1 ];

f1 = plot_results(1,titles,acc3,pos_titles,c); % main figure 1
sgtitle('Gradiometer (Alpha)')
sgt.FontSize = 20;