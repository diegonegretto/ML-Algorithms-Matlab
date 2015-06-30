%   Algorithm "Label Propagation", by Zhu and Ghahramani (2002)
%   
%   Usage: owner = labelprop(X,slabel,sigma,disttype,nclass,iter)
%   
%   X = attributes vector (line = elements, columns = attributes)
%   slabel = vector with numerical labels (>0) of pre-labeled elements (0 for
%          unlabeled elements)
%   nclass = number of classes
%   iter = number of interations
%   sigma = width of the Gaussian kernel


function owner = labelprop(X,slabel,sigma,disttype,nclass,iter)
    
    if (nargin < 6) || isempty(iter),
        iter = 10000; 
    end
    
    if (nargin < 5) || isempty(nclass),
        nclass = max(slabel); 
    end
    
    if (nargin < 4) || isempty(disttype),
        disttype = 'euclidean';
    end
    
    qtnode = size(X,1);  
    
    W = exp(-squareform(pdist(X,disttype).^2)/2*sigma^2); 
    
    D = sparse(diag(sum(W,2)));
    
    Y = zeros(qtnode,nclass); 
    
    labelednodes = find(slabel);
    
    Y(sub2ind(size(Y),labelednodes,slabel(labelednodes))) = 1;
    
    noch = 0; 
    
    [~,owner] = max(Y,[],2);
    
    DW = D\W;
    
    clear D W;
    
    YI = Y;
    
    YTP = repmat(slabel~=0,1,nclass);
    
    YTN = 1-YTP;
    
    for j=1:iter  
        
        Y = DW * Y;
        Y = YTP .* YI + YTN .* Y;
        ownerbak = owner;
        [~,owner] = max(Y,[],2);
        
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