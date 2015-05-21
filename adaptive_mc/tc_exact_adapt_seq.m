function [ M_hat , obs] = tc_exact_adapt_seq( M, m_array )
%TC_EXACT_ADAPT : tensor completion in Arti Singh 2013 NIPS paper
% exact version, sequential program 
%%%%%%%%%%%%%%%
% X: tensor 
% X_hat: sampled sparse tensor
% get the size of the tensor


T = ndims(M); 
sz_array = size(M);
M_hat = zeros(sz_array);
obs = zeros(sz_array);

% randomly draw entries for T's mode
n_neg_T =  prod(sz_array(1:T-1));
U = zeros(n_neg_T,1);

    % for each sub-tensor, fixing the T's mode
    for  i = 1: sz_array(end)
        Omega  = unique( randi(n_neg_T, [int32( m_array(T) * n_neg_T  ), 1]) );
        M_i = M (:,:,i);   
        M_i_Omega =  M_i( Omega ); %sub-tensor with sampled entries, turn into a vector


        U_Omega = U (Omega );
        M_i_project = orth_proj(U_Omega)* M_i_Omega;

        if norm( M_i_Omega - M_i_project)>0 % test whether M is in the subspace of U      
    %         whos('M_i');
            [M_hat_i, obs_i] =  mc_exact_adapt (M_i, m_array(i)) ;% recursive
            
          % vectorize M_hat_i
            U  = [ U,  M_hat_i(:)]; % vectorize
            U = orthogonalize(U);
%             disp(U);
%             U_null = null(U ,'r'); % null space of U    
%             M_hat_project = orth_proj(U_null) * M_hat_i;
%             U_i = M_hat_project/ norm(M_hat_project ); % new subspace tensor
            
        else

            M_hat_i = U / ( U_Omega' * U_Omega )* U_Omega * M_i_Omega;
        end

         M_hat(:,:,i) = reshape(M_hat_i, sz_array(1:end-1));
         obs(:,:,i) = reshape(obs_i, sz_array(1:end-1));
    end
end

function P = orth_proj(U)
    % construct an orthogonal projection matrix from basis 
    if all(U==0)
        P  = zeros(size(U,1));
        return
    end
    if(size(U,2)> size(U,1))
        P = U * pinv(U'*U) * U';

    else
        P = U / (U'*U) *U';
    end
end

function [U] = orthogonalize(U)
   [u,s,v] = svd(U);
   U = u(:,diag(s)>0.001);
end