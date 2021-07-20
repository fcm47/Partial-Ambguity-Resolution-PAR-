function [sqnorm,Ps,nfixed,zfixed, iZt, delete_index_in_origin]=myparsearch(ahat, Qahat, record, zhat,Qzhat,Z,L,D,P0, ncands)
%
%   [zpar,sqnorm,Qzpar,Zpar,Ps,nfixed,zfixed]=parsearch(zhat,Qzhat,Z,L,D,P0,ncands)
%
%This routine performs an integer bootstrapping procedure for partial 
%ambiguity resolution (PAR) with user-defined success-rate P0
%
%INPUTS:
%   zhat: decorrelated float ambiguities (zhat, must be a column!)
%  Qzhat: variance-covariance matrix of decorrelated float ambiguities
%      Z: Z-matrix from decorrel
%    L,D: lower-triangular and diagonal matrix from LtDL-decomposition of Qzhat
%     P0: Minimum required sucess rate [DEFAULT=0.995]
% ncands: Number of requested integer candidate vectors [DEFAULT=2]
%
%OUTPUTS:
%    zpar: subset of fixed ambiguities (nfixed x ncands) 
%  sqnorm: squared norms corresponding to fixed subsets
%   Qzpar: variance-covariance matrix of float ambiguities for the
%          subset that is fixed
%    Zpar: Z-matrix corresponding to the fixed subset
%      Ps: Bootstrapped sucess rate of partial ambiguity resolution
%  nfixed: The number of  fixed ambiguities
%  zfixed: [OPTIONAL]: also give the complete 'fixed' ambiguity vector
%          where the remaining (non-fixed) ambiguities are adjusted 
%          according to their correlation with the fixed subset
%
%
% NOTE:
% This PAR algorithm should be applied to decorrelated ambiguities, which
% can be obtained from the original ahat and its variance-covariance matrix 
% Qahat as:
% 
%     [Qzhat,Z,L,D,zhat] = decorrel(Qahat,ahat);
%
% The fixed baseline solution can be obtained with (see documentation):
%
%     s      = length(zhat) - nfixed + 1;
%     Qbz    = Qba * Zpar ;   
%     bfixed = bhat - Qbz/Qzpar * (zhat(s:end)-zpar(:,1) );
%     Qbfixed= Qbhat - Qbz/Qzpar * Qbz';
%
% Hence, zfixed is not required. LAMBDA, however, does give the solution
% in terms of the full ambiguity vector.
%
%   *************************************************************
%   *                Author: Sandra Verhagen                    *
%   *                  Date: 04/APRIL/2012                      *
%   *                  GNSS Research Centre                     *
%   *  Dept.of Spatial Science, Curtin University of Technology *
%   *************************************************************

%============================START PROGRAM==========================%
if(nargin <4 )
    P0=0.999;
    warning(['user-defined success rate is necessary for PAR method',...
        'the default value 0.995 is used']);
end
if nargin<6
    ncands = 2;
end
global Thres_ratio
n      = size (Qzhat,1);
%% 
% bootstrapped success rate if all ambiguities would be fixed
delete_index_in_origin =[];
ahat_old = ahat;
for k=1:n-4
    if k==1
        [Qzhat,Z,L,D,zhat,iZt, record] = decorrel(Qahat,ahat);
    elseif k>1
        %         find the maximum in decorrelated Qvar
        [val, ind]=max(diag(Qzhat)); 
        %         find the maximum in  corresponding Qahat according to the
        %         record when decorrelaing
        Qahat_index = record(ind);
        %       delete corresponding elements in ahat and Qahat
        deleta_val = ahat(Qahat_index);
        ahat = [ahat(1:Qahat_index-1); ahat(Qahat_index+1:end)];
        Qahat = [Qahat(:, 1:Qahat_index-1), Qahat(:, Qahat_index+1:end)];
        Qahat = [Qahat(1:Qahat_index-1, :); Qahat(Qahat_index+1:end, :)];
        %         get its index in original ahat
        ahat_index_old = find(ahat_old==deleta_val);
        %         recording the index of deleted ambiguities in original ahat, 
        %         it would be useful while convert to fixed solution
        delete_index_in_origin = [delete_index_in_origin, ahat_index_old];
        [Qzhat,Z,L,D,zhat,iZt, record] = decorrel(Qahat,ahat);
    end
    
    Ps = prod ( 2 * normcdf(1./(2*sqrt(D(k:end)))) -1 ); 
    if(Ps<P0)
        continue;
    end

    
    % last n-k+1 ambiguities are fixed to integers with ILS
    [zfixed,sqnorm] = ssearch(zhat,L,D,ncands);
    if(sqnorm(2)/sqnorm(1)>Thres_ratio)
        
        nfixed = n-k+1;
        return;
    end 
end

sqnorm = [];
Ps     = NaN;
zfixed = zhat;
nfixed = 0;


return;
