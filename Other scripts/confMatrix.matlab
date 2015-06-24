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
%
% Copyright 2015 Diego H. Negretto
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function confMatrix(lb, result)

    c = zeros(max(lb));

    for i=1:size(lb,1)
       
       c(lb(i),result(i)) = c(lb(i),result(i)) + 1;
      
    end
    
   
   disp('Confusion Matrix ');
   disp(c);

end