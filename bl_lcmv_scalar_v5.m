function [wts,Hscalar,ori,P]=bl_lcmv_scalar_v5(H,H_idx,Href,RNcov,ref_idx,loc_flag)
%function [wts,Hscalar,ori,P]=bl_lcmv_scalar_v2(H,H_idx,Href,RNcov,ref_idx,loc_flag);
%
% INPUT:
%   H =  lead field  [channels x 3 orientations x voxels]
%   H_idx = indices of voxel to calculate --> must exclude the ref_idx indices from the H_idx .
%   Href = lead field for voxels to place nulls (can use wts.scalar_LeadFields from "BRANElab_LCMV_beamformer_singlesource_v3.m".).
%   ref_idx = indices of the mHref nulls -- to start with single-source search mHref=[]; mref=[];.
%   loc_flag =  (0) returns 'pseudoZ' for single-source or 'MPZ' for multisource.
%                  (1) returns 'pseudoZ_er' for single-source or 'MER' for multisource.
%                  (2) returns 'pseudoZ_rer' for single-source or 'rMER' for multisource.
%                  (3) returns 'activity_index' for single-source or 'MAI' for multisource.
%           Note: single-source localizer returned if Href=[]; or single-source localizer returned if Href=[channels x voxels].  
%
% OUTPUT:
%   wts = source weights [channels x voxels];
%   Hscalar = scalar lead-field matrix based on single- or mult-source
%   ori = orientation of voxels [voxels x (X,Y,Z)];
%   P.type = possible localizers -  formulae from Moiseev et al., 2011 .
%           pseudoZ       = (H'*Rinv*H)/(H'*Rinv*N*Rinv*H);
%           pseudoZ_er   = (H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*N*Rinv*H)
%           pseudoZ_rer   = (H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*H)
%           activity_index  = (H'*Ninv*H)/(H'*Rinv*H)
%           MPZ             = trace((H'*Rinv*H)/(H'*Rinv*N*Rinv*H))-n,    where n = number of multiple sources.
%           MER             =  trace((H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*N*Rinv*H))-n,    where n = number of multiple sources.
%           rMER             =  trace((H'*Rinv*Rbar*Rinv*H)/(HF'*Rinv*HF))-n,    where n = number of multiple sources.
%           MAI             =  trace((H'*Ninv*H)/(HF'*Rinv*HF))-n,    where n = number of multiple sources.
%   P.img = image for localizer values [voxels x 1]
%
%   written  by A. Herdman Jun 13, 2015
%       updates:
%           Aug 15, 2015 - added MAI and correct bug in "bl_lcmv_scalar.m" - now performs eigen solution for localizer speficied. 
warning('off','all');
[num_chans,num_ori,num_voxels]=size(H);
wts=nan(num_chans,num_voxels);    % initializing LCMV weights
Hscalar=nan(num_chans,num_voxels); % initializing scalar lead-field matrix
ori=nan(num_voxels,3); % initializing oreientaiton matrix
img=nan(1,num_voxels); 
P.type=''; P.img=zeros(num_voxels,1);
if size(H,2)~=3;
    fprintf('Terminated program because H is not in the proper format.\n H =  lead field  [channels x 3 orientations x voxels]\n'); 
    return;
end
R=RNcov.R; Rinv=RNcov.Rinv; Rbar=RNcov.Rbar; N=RNcov.N; Ninv=RNcov.Ninv; Nbar=RNcov.Nbar;
N2=N/size(Href,2);
vv=setdiff(H_idx,ref_idx); % exlcuding already found sources in ref_idx
%vv=1:size(H,3);
%   S=[];Skk=[]; Tkk=[]; D=[]; F=[]; Skr=[]; TkR=[]; SRR=[]; TRR=[];
parfor v=1:size(H,3); %1:length(H_idx); %1:size(H,3);
    if ismember(v,vv);
        if ~isempty(Href);  % Multi-source LCMV
            Hk = H(:,:,v);
            % for different localizers
            if loc_flag==0;     % MPZ
                Skk = (Hk'*Rinv*Hk); SkR = (Hk'*Rinv*Href); SRR = (Href'*Rinv*Href); SRk = (Href'*Rinv*Hk);     % S
                Tkk = (Hk'*Rinv*N*Rinv*Hk); TkR = (Hk'*Rinv*N*Rinv*Href); TRR = (Href'*Rinv*N*Rinv*Href); TRk = (Href'*Rinv*N*Rinv*Hk); % T
            elseif loc_flag==1;    % MER
                Skk = (Hk'*Rinv*Rbar*Rinv*Hk); SkR = (Hk'*Rinv*Rbar*Rinv*Href); SRR = (Href'*Rinv*Rbar*Rinv*Href); SRk = (Href'*Rinv*Rbar*Rinv*Hk);    % E
                Tkk = (Hk'*Rinv*N*Rinv*Hk); TkR = (Hk'*Rinv*N*Rinv*Href); TRR = (Href'*Rinv*N*Rinv*Href); TRk = (Href'*Rinv*N*Rinv*Hk);     % T
            elseif loc_flag==2;    % rMER
                Skk = (Hk'*Rinv*Rbar*Rinv*Hk); SkR = (Hk'*Rinv*Rbar*Rinv*Href); SRR = (Href'*Rinv*Rbar*Rinv*Href); SRk = (Href'*Rinv*Rbar*Rinv*Hk);    % E
                Tkk = (Hk'*Rinv*Hk); TkR = (Hk'*Rinv*Href); TRR = (Href'*Rinv*Href); TRk = (Href'*Rinv*Hk); % Sred
            elseif loc_flag==3;    % MAI
                Skk = (Hk'*Ninv*Hk); SkR = (Hk'*Ninv*Href); SRR = (Href'*Ninv*Href); SRk = (Href'*Ninv*Hk);    % G
                Tkk = (Hk'*Rinv*Hk); TkR = (Hk'*Rinv*Href); TRR = (Href'*Rinv*Href); TRk = (Href'*Rinv*Hk); % S
            end
            
            % for multisource nulling -- see Mosieev et al, 2011 (page 486 appendix.
            D = (TkR*inv(TRR)*SRR*inv(TRR)*TRk)-(TkR*inv(TRR)*SRk)-(SkR*inv(TRR)*TRk)+Skk;  % D = S
            F = Tkk-(TkR*inv(TRR)*TRk); % F = T
            
            %%% solving traditional eigenvalue problem: S*u=lamda_max[T*u*] , where lamda_max = maximum eigenvalue. --> [uvect,uval] = eigen(Sinv*T) .
            %[eigvect,eigval] = eig(F\D); % same as [eigvect,eigval] = eig(inv(F)*D);
%               try
%                 [eigvect,eigval] = eig(F\D); % same as [eigvect,eigval] = eig(inv(F)*D);
%             catch
                [eigvect,eigval] = eig(pinv(F)*D); % using pseudo-inverse becuase F is badly scaled;
%             end
            %         x=diag(eigval); % some eigval were complex thus only selecting real eigen values
            %         if ~isreal(x);
            %             for xt=1:length(x); if isreal(x(xt))==1; t(xt)=1; end; end;
            %             t_idx=find(t==1);
            %             [~,w]=max(x(t_idx));  e_idx=t_idx(w);
            %         else
            [~,e_idx]=max(diag(eigval));   % find maximum eigenvalue
            %         end
            ori(v,:)=real(eigvect(:,e_idx));  % source orientiation vectors - columns of eigvect correspond to eigenvectors
            
            %This converts the 3 vector dipoles to the single scalar dipole orientation.
            h = H(:,:,v)*ori(v,:)';     % lead field for dipole with best orientation (yields highest signal).
            h = h/norm(h);              % normalizes the leadfield across channels for each voxel (Needed to get proper weights)
            h = h-mean(h);              % average referrencing lead fields
            Hscalar(:,v)=h;             % lead-field for dipole with best orientation (i.e., with largest eigenvalue).
            
            wts(:,v)=(Rinv*h)/(h'*Rinv*h); % Alex's suggested script - returns same as nutMEG script
            HF=[h,Href];    % multi-source lead field
            
            % Multi-source Localizers
            if loc_flag==0;     % MPZ
                n=size(HF,2);
                S = (HF'*Rinv*HF);
                T = (HF'*Rinv*N*Rinv*HF);
                img(v)=trace(S/T)-n;
            elseif loc_flag==1;    % MER
                E = (HF'*Rinv*Rbar*Rinv*HF);
%                 T = (HF'*Rinv*N*Rinv*HF);
                T = (HF'*Rinv*N2*Rinv*HF);
                img(v)=trace(E/T);    % forward slash "/" means E\Tinv
            elseif loc_flag==2;    % rMER
                E = (HF'*Rinv*Rbar*Rinv*HF);
                Sred = (HF'*Rinv*HF); % for reduced event-related
                img(v)=trace(E/Sred);
            elseif loc_flag==3;    % MAI
                n=size(HF,2);
                G = (HF'*Ninv*HF);
                S = (HF'*Rinv*HF);
                img(v)=trace(G/S)-n;
            end
            
            
        else % Single-source LCMV
            if loc_flag==0; % pseudoZ
                S = (H(:,:,v)'*Rinv*H(:,:,v));                     % S
                T = (H(:,:,v)'*Rinv*N*Rinv*H(:,:,v));          % T
            elseif loc_flag==1; % pseudoZ_er
                S = (H(:,:,v)'*Rinv*Rbar*Rinv*H(:,:,v));     % E
                T = (H(:,:,v)'*Rinv*N*Rinv*H(:,:,v));          % T
            elseif loc_flag==2; % pseudoZ_rer
                S = (H(:,:,v)'*Rinv*Rbar*Rinv*H(:,:,v));      % E
                T = (H(:,:,v)'*Rinv*H(:,:,v));                      % Sred
            elseif loc_flag==3; % activity index - Van Veen et al., 1997
                S = (H(:,:,v)'*Ninv*H(:,:,v));                      % G
                T = (H(:,:,v)'*Rinv*H(:,:,v));                      % S
            end
            
            
            %%% solving trasitional eigenvalue problem: S*u=lamda_max[T*u*] , where lamda_max = maximum eigenvalue. --> [uvect,uval] = eigen(Sinv*T) .
%             [eigvect,eigval] = eig(T\S);    % same as [eigvect,eigval] = eig(Tinv*S);
%             try
%                 [eigvect,eigval] = eig(T\S);  % same as [eigvect,eigval] = eig(Tinv*S);
%             catch
                [eigvect,eigval] = eig(pinv(T)*S);    % using pseudo-inverse because inv(T) returns inf or nan
%             end
%             
            %         x=diag(eigval); % some eigval were complex thus only selecting real eigen values
            %              keyboard;
            %         if ~isreal(x);
            %             fprintf('error\n'); return;
            %             for xt=1:length(x); if isreal(x(xt))==1; t(xt)=1; end; end;
            %             t_idx=find(t==1);
            %             [~,w]=max(x(t_idx));  e_idx=t_idx(w);
            %         else
            [~,e_idx]=max(diag(eigval));   % find maximum eigenvalue
            %         end
            ori(v,:)=real(eigvect(:,e_idx));  % source orientiation vectors - columns of eigvect correspond to eigenvectors
            
            %This converts the 3 vector dipoles to is the single scalar dipole orientation.
            h = H(:,:,v)*ori(v,:)';     % lead field for dipole with best orientation (yields highest signal).
            %        Hscalar(:,v)=h;             % lead-field for dipole with best orientation (i.e., with largest eigenvalue).
            h = h/norm(h);              % normalizes the leadfield across channels for each voxel (Needed to get proper weights)
            h = h-mean(h);              % average referrencing lead fields
            Hscalar(:,v)=h;             % lead-field for dipole with best orientation (i.e., with largest eigenvalue).
            
            w=(Rinv*h)/(h'*Rinv*h); % Alex's suggested script - returns same as nutMEG script
            wts(:,v)=w; HF=h;
            % Single-source Localizers
            if loc_flag==0; % pseudoZ
                img(v)= (HF'*Rinv*HF)/(HF'*Rinv*N*Rinv*HF);    % equivalent to (w'*R*w)\(w'*N*w);
            elseif loc_flag==1; % pseudoZ_er
                E = (HF'*Rinv*Rbar*Rinv*HF);
                T = (HF'*Rinv*N*Rinv*HF);
                img(v)= E/T;
            elseif loc_flag==2; % pseudoZ_rer
                E = (HF'*Rinv*Rbar*Rinv*HF);
                Sred = (HF'*Rinv*HF); % for reduced event-related
                img(v)= E/Sred;   % similar to rMER
            elseif loc_flag==3; % activity index - Van Veen et al., 1997
                E = (HF'*Ninv*HF);
                Sred = (HF'*Rinv*HF); % for reduced event-related
                img(v)= E/Sred;   % similar to rMER
            end
        end
    end
end

% Multi-source Localizers
if loc_flag==0 && ~isempty(Href);     % MPZ
    P.type='MPZ';
elseif loc_flag==1 && ~isempty(Href);    % MER
    P.type='MER';
elseif loc_flag==2 && ~isempty(Href);    % rMER
    P.type='rMER';
elseif loc_flag==3 && ~isempty(Href);    % MAI
    P.type='MAI';
elseif loc_flag==0 && isempty(Href);     % pseudoZ
    P.type='pseudoZ';
elseif loc_flag==1 && isempty(Href);    % pseudoZ_er 
    P.type='pseudoZ_er';
elseif loc_flag==2 && isempty(Href);    % pseudoZ_rer
    P.type='pseudoZ_rer';
elseif loc_flag==3 && isempty(Href);    % Activity_Index
    P.type='Activity_Index';
end

P.img=img';
%[wts,Hscalar,ori,P]

