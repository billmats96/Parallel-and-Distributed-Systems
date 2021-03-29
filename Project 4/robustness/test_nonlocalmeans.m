
%---------------test_nonlocalmeans.m Script-----------------
% This script generates random grayscale image for denoising. Then we
% execute demo_non_local_means.m with the random image. We compare matlab
% demo output with the ouput of our testing code using the same data as
% input. If the results satisfy identification prerequisit ( norm of
% difference of ouputs <0.1) then test is successful. Parametres: nLevel,
%  patchSize, filtSigma, patchSigma must be same with those in test code.
%
% 
% Authors:  Athanasiadis Christos   athanasc@ece.auth.gr   AEM 8401
%                 Matsoukas Vasileios     vmatsouk@ece.auth.gr  AEM 8743
%
% demo_non_local_means.m by Dimitrios Floros
%

close all
clear 

%Define parameters
im_rows=60;
im_cols=50;
test_image=rand(im_rows,im_cols);

% noise
noiseParams = {'gaussian', 0, 0.001}; 

% number of regions
nLevel = 6; 

% filter sigma value
filtSigma = 0.1;
patchSize = [5 5];
patchSigma = 5/3;

executable_file="serialnlm_adaptive"; %choose between serialnlm, serialnlm_adaptive, nlm, adaptive

save test_nlm test_image 

%Run matlab non local means script
run('demo_non_local_means.m')


%Save random input (noisy)  as a binary file
A=J; %J is the noisy image produced by the demo
A=A';
fileID = fopen('randInput.bin','w');
fwrite(fileID,A,'float');
fclose(fileID);

%Obtain matlab's output and write to a binary file
if (executable_file=="serialnlm" || executable_file=="nlm")
    Out=If;
elseif (executable_file=="serialnlm_adaptive" || executable_file=="adaptive")
    Out=Ia;
end
    
% Test the executable with the random binary input produced above. 
% Compare results with output produced by matlab script.
fprintf("Testing dataset on: %s\n",executable_file);

im_rows=num2str(im_rows);
im_cols=num2str(im_cols);
binaryInput="randInput.bin";

command=['./',char(executable_file),' ',char(binaryInput),' ','outputTEST.bin ',' ',char(im_rows),' ',char(im_cols)];
[status,cmdout] = system(command);
cmdout   
 
im_rows=str2double(im_rows);
im_cols=str2double(im_cols);

%Get output of testing code and compare with matlab results. 
fileID = fopen('outputTEST.bin','r');
B=fread(fileID,[im_cols im_rows],'float');
B=B';
fclose(fileID);

C=B-Out; %C matrix holds the difference of matlab output and testing code output
fprintf('Norm of comparison matrix C is: %f\n',norm(C));

if(norm(C)>0.1)
    fprintf('Considerable large value of norm(C)\n');
else
    fprintf('Accepted norm value\n'); 
end

fprintf('End of test\n');







