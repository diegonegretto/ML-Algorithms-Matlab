% Script to create a confusion matrix
%
% Developed by: Diego Henrique Negretto
% Date: 07/04/2015
% 
% Script to create a confusion matrix from output vector machine learning algorithms
%
% Usage: confMatrix(lb, result);
% lb = vector with the original labels
% result = vector with the generated labels 

function confMatrix(lb, result)

    c = zeros(max(lb));

    for i=1:size(lb,1)
       
       c(lb(i),result(i)) = c(lb(i),result(i)) + 1;
      
    end
    
   
   disp('Confusion Matrix ');
   disp(c);

end