close all;clear;clc

%%% Charger une image
load('DTI_6dir.mat')
I = DTI_Data;
I0 = I(:,:,:,1);
% figure; imshow(I0(:,:,35),[]);

%-fanDTasia ToolBox------------------------------------------------------------------
% This Matlab script is part of the fanDTasia ToolBox: a Matlab library for Diffusion 
% Weighted MRI (DW-MRI) Processing, Diffusion Tensor (DTI) Estimation, High-order 
% Diffusion Tensor Analysis, Tensor ODF estimation, Visualization and more.
%
% A Matlab Tutorial on DW-MRI can be found in:
% http://www.cise.ufl.edu/~abarmpou/lab/fanDTasia/tutorial.php
%
%-CITATION---------------------------------------------------------------------------
% If you use this software please cite the following work:
% A. Barmpoutis and B.C. Vemuri, "A Unified Framework for Estimating Diffusion Tensors 
% of any order with Symmetric Positive-Definite Constraints", 
% In the Proceedings of ISBI, 2010
%
%-DESCRIPTION------------------------------------------------------------------------
% This demo script shows how to compute a field of Diffusion Tensors from a given DW-MRI
% (Diffusion-Weighted MRI) dataset. The method guarantees that the estimated tensors
% are positive-definite or at least positive semi-definite. Here the given demo dataset 
% is simulated using the DW-MRI simulator developed by Angelos Barmpoutis.
%
%-USE--------------------------------------------------------------------------------
% D=DEMO_DTI_Field_Estimation;
%
% D: is the computed 2D Diffusion Tensor field in the form of a 3x3xSizeXxSizeY matrix
%
%-DISCLAIMER-------------------------------------------------------------------------
% You can use this source code for non commercial research and educational purposes 
% only without licensing fees and is provided without guarantee or warrantee expressed
% or implied. You cannot repost this file without prior written permission from the 
% authors. If you use this software please cite the following work:
% A. Barmpoutis and B.C. Vemuri, "A Unified Framework for Estimating Diffusion Tensors 
% of any order with Symmetric Positive-Definite Constraints", In Proc. of ISBI, 2010.
%
%-AUTHOR-----------------------------------------------------------------------------
% Angelos Barmpoutis, PhD
% Computer and Information Science and Engineering Department
% University of Florida, Gainesville, FL 32611, USA
% abarmpou at cise dot ufl dot edu
%------------------------------------------------------------------------------------

% order=2;%In standard DTI the order is 2
% 
% %Here is a sample demo DW-MRI dataset of 1 voxel (21 gradient orientations)
% %This dataset was synthesized using the DW-MRI simulator developed by Angelos Barmpoutis.
% fprintf(1,'Simulating DW-MRI dataset...');
% UnitVectors;
% 
% GradientOrientations = [-1 0 1; 1 0 1; 0 1 1; 0 1 -1; -1 1 0; 1 1 0];
% b_value = ones(6,1);
% S = I(:,:,:,(2:7));
% 
% S = bsxfun(@rdivide, S(:,:,:,(1:6)), I0);
% 
% %Construct all possible monomials of a specific order
% G=constructMatrixOfMonomials(GradientOrientations, order); %computes G from section 5.1 (ISBI'10)
% %Construct set of polynomial coefficients C
% C=constructSetOf81Polynomials(order)'; %computes C from section 5.1 (ISBI'10)
% P=G*C;
% P=[-diag(b_value)*P ones(size(GradientOrientations,1),1)];
% 
% 
% %The next lines implement the core of the algorithm. 
% %It should be repeated for every voxel in the DW-MRI dataset (here there is only 1 voxel).
start_time=cputime;
% 
% % points 39,44,35
% 
% for k=1:size(S,3)
%     for i=1:size(S,1)
%         for j=1:size(S,2)
%         
%             y=squeeze(log(S(i,j,k,:)));
%             x=lsqnonneg(P, y);
%             UniqueTensorCoefficients = C * x([1:81]);
% 
%             %Put the result in the form of a 3x3 matrix
%             T=[UniqueTensorCoefficients(6) UniqueTensorCoefficients(5)/2 UniqueTensorCoefficients(4)/2
%                UniqueTensorCoefficients(5)/2 UniqueTensorCoefficients(3) UniqueTensorCoefficients(2)/2
%                UniqueTensorCoefficients(4)/2 UniqueTensorCoefficients(2)/2 UniqueTensorCoefficients(1)];
% 
%             D(:,:,i,j,k) = T;
%             TensorCoefficients(:,i,j,k) = UniqueTensorCoefficients;
%             
%             [V,Diag] = eig(T);
%             eigenValues(i,j,k,:) = [Diag(1,1) Diag(2,2) Diag(3,3)];
%             eigenVectors(i,j,k,:) = V(:,3);
%             
%         end
%     end
%     fprintf(1,'\n %.0f is done',k);
% end
% 
% % Save
% save('ellipsoideMatrix.mat', 'D');
% save('eigenValues.mat', 'eigenValues');
% save('eigenVectors.mat', 'eigenVectors');

load('ellipsoideMatrix.mat');
% load('eigenValues.mat');
% load('eigenVectors.mat');

end_time=cputime;
fprintf(1,'\nTotal estimation time: %.0f ms\n\n',(end_time-start_time)*1000);

 tractography(D,I0,10000);

% If you want to plot a tensor or a tensor field as an ellipsoid or a field of ellipsoids
% you have to download the plotDTI.m function developed by Angelos Barmpoutis, Ph.D.
% and then uncomment the following line.

%plotDTI(D);

% or if you want to plot a tensor or a tensor field as spherical functions
% you have to download the plotTensors.m function developed by Angelos Barmpoutis, Ph.D.
% and then uncomment the following line.

%plotTensors(TensorCoefficients,0.5,[321 1]);