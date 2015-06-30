%   Algorithm Linear Neighborhood Propagation by Fei Wang et. al.
%
%   Usage: owner = lnp(X,slabel,k,alpha,nclass,iter)
%
%   X = attributes vector (line = elements, columns = attributes)
%   slabel = vector with numerical labels (>0) of pre-labeled elements (0 for
%          unlabeled elements)
%   nclass = number of classes
%   iter = number of interations
%   alpha = in the range [ 0 1 ] , defines relative amount of information coming from neighbors and the initial information
%   k = number of nearest neighbors in the formation of graph

function owner = lnp(X,slabel,k,disttype,alpha,nclass,iter)
  
    if (nargin < 7) || isempty(iter),
        iter = 10000; 
    end
    
    if (nargin < 6) || isempty(nclass),
        nclass = max(slabel);
    end
    
    if (nargin < 5) || isempty(alpha),
        alpha = 0.99;
    end
    
    if (nargin < 4) || isempty(disttype),
        disttype = 'euclidean'; 
    end
    
    [qtnode,qtfeat] = size(X);  
    
    [~,IX] = sort(squareform(pdist(X,disttype).^2)); 
    
    IX = IX(2:k+1,:); 
    
    W = zeros(qtnode); 
    
    parfor i=1:qtnode
       N = X(IX(:,i),:)'; 
       Gp = ((X(i,:)' * ones(1,k)) - N);
       G = Gp' * Gp;
      
       if k>qtfeat
            G = G + 0.1 * trace(G) * eye(k); 
       end
      
       options = optimset('Algorithm','active-set','Display','off');
      
       w = quadprog(G,ones(1,k),[],[],ones(1,k),1,zeros(k,1),ones(k,1),[],options);
       
       wi = zeros(1,qtnode); 
       
       wi(IX(:,i)) = w;
       
       W(i,:) = wi;    
    
    end
    
    clear N Gp G w wi;
    
    Y = zeros(qtnode,nclass);
    
    labelednodes = find(slabel);
    
    Y(sub2ind(size(Y),labelednodes,slabel(labelednodes))) = 1;
    
    noch = 0; 

    F = Y;
    
    [~,owner] = max(F,[],2);
    
    alphaW = alpha * W;
    
    clear W;
    
    alphaY = (1 - alpha) * Y;
    
    clear Y;
    
    for i=1:iter 
         F = alphaW * F + alphaY;
         ownerbak = owner;
         [~,owner] = max(F,[],2);
         if sum(ownerbak~=owner)==0  
             noch = noch + 1;
             if noch>=100  
                 break; 
             end
         else
             noch = 0;
         end
   
    end

end