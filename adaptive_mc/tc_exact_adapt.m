function [ X_hat ] = tc_exact_adapt( X )
%TC_EXACT_ADAPT : tensor completion in Arti Singh 2013 NIPS paper
% exact version
%%%%%%%%%%%%%%%
% X: tensor 
% X_hat: sampled sparse tensor

% get the size of the tensor
T = ndims (X);
m_array = zeros ( 1, T);
for t = 1: T
    m_array(t) = 0.3;
end
[X_hat ] = sequetial_tc(X,m_array);

end

function [M_hat ] = sequetial_tc(M,m_array)
% implemention of the qseudo code in sequential tensor completion
% m_t: sampling rate for the t's mode

if numel(M)==1
    M_hat = M;
    return;
end
U = [];
T = ndims(M); 

sz_array = size(M);
M_hat = zeros(sz_array);


if isvector(M) || numel(M)==1
 % fully observed 
  Omega = 1:sz_array(end);
else
   % randomly draw entries for T's mode
   n_neg_T =  prod(sz_array(1:T-1));
   Omega  = unique( randi(n_neg_T, [int32( m_array(T) * n_neg_T  ), 1]) );

end

% if isvector(M)
%     fprintf('T = 1\n');
% else
%     fprintf('T = %d \n',T );
% end

% for each sub-tensor, fixing the T's mode
for   i = 1: sz_array(end)
  
    if (T ==3)
         M_i = M (:,:,i);
    elseif isvector(M);
         M_i = M(i);     
    elseif ismatrix(M) % this contains the vector,
        M_i = M (:,i);
    end
    
   
    M_i_Omega =  M_i( Omega ); %sub-tensor with sampled entries
     
    
    if isempty (U) % initial state
       P_U_Omega  = 0;
%        U = M_i;
    else
       U_Omega = U (Omega );
       P_U_Omega = U_Omega / ( U_Omega' * U_Omega ) * U_Omega';         
    end
    
    M_i_project = P_U_Omega* M_i_Omega;
   
    if norm( M_i_Omega - M_i_project)>0 % test whether M is in the subspace of U      
%         whos('M_i');
        M_hat_i =  sequetial_tc (M_i, m_array) ;% recursive
%         disp(M_hat_i);
%         disp(M_i);
        U_null = null(U ,'r'); % null space of U
        P_U_null = U_null / ( U_null' * U_null )  * U_null';
        
%         whos('P_U_null')
%         whos('M_hat_i')
        if P_U_null~=0
            M_hat_project = P_U_null * M_hat_i;
            U_i = M_hat_project/ norm(M_hat_project ); % new subspace tensor
            U  = [ U; U_i(:)]; % vectorize
        end
    else
        
        M_hat_i = U / ( U_Omega' * U_Omega )* U_Omega * M_i_Omega;
    end
    
    if ( T == 3)
        M_hat(:,:,i) = M_hat_i;
    else
        M_hat(:,i)  = M_hat_i;
    end
    
end
    
end
