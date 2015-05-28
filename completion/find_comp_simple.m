function [K,C] =  find_comp_simple(G,r)

v_s = []; % set of removed vertices
degree = sum(G);
v = find(degree < r );
v_s=[v_s,v];
while( ~isempty(v))
    % remove vertex
    G(v,:) = 0;
    G(:,v) = 0; 
    % check degree
    degree = sum(G);
    v = setdiff(v_s, find(degree < r ));    
    v_s = [v_s, v];
%     disp(length(v_s));
end

if(length(v_s) ==size(G,1))
    fprintf('No completable submatrix\n');
end

%% find strongly-connected component

graph = sparse(G);
[K,C] = graphconncomp(graph); % BUG: too many singletons

end