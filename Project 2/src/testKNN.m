%---------------testKNN.m Script-----------------
% This script generates datasets containig random real positive
% numbers. We execute knn using knn_example.m and compare matlab output 
% with our testing code. 
%Author: Matsoukas Vasileios vmatsouk@auth.gr

%knn_example.m by Dimitrios Floros 

%%
%Clean-up
clear
close all

%Define parameters
dataset_rows=1000;
dataset_cols=20;
kNbr=3;
processes="8"; 
executable_file="mpi_blocking"; %use one of the following "serialKNN",
                                                 %"mpi_blocking",
                                                 %"mpi_nonblocking",
                                                 %"mpi_threaded_blocking",
                                                 %"mpi_threaded_nonblocking"


%Create a random sparse dataset and save as -> train_testKNN.mat
%A=abs(rand(dataset_rows,dataset_cols));
A=abs(sprand(dataset_rows,dataset_cols,0.15)); %Sparse arrays are mainly used in knn application. 
A=full(A);
labels=ones(dataset_rows,1);
train_labels=labels;
train_X=A;
save train_testKNN train_X train_labels


%save random input as a binary file
A=A';
fileID = fopen('randInput.bin','w');
fwrite(fileID,A,'float');
fclose(fileID);

%run matlab knn script
run('knn_example.m')

%obtain matlab's output distances and write to a binary file
DIST=single(DIST);
DIST=DIST.^2;
DIST=DIST';
fileID = fopen('outputDISTtestKNN.bin','w');
fwrite(fileID,DIST,'float');
fclose(fileID);

%obtain matlab's output labels and write to a binary file 
IDX=IDX';
fileID=fopen('outputIDXtestKNN.bin','w');
fwrite(fileID,IDX,'int32');
fclose(fileID);

%Test the executable with the random binary input produced above and compare 
%results with output binary files produced by matlab script
%execute and check output in matlab command window
fprintf("Testing dataset on: %s\n",executable_file);
k=num2str(kNbr);
dataset_cols=num2str(dataset_cols);
if executable_file=="serialKNN"
    %%if you want to test serialknn use this section of code
    command=['./serialKNN ', char(k),' ',char(dataset_cols)];
    
else
    %%use this section of code for mpi versions.
    system('export I_MPI_SHM_LMT=shm');
    command = ['export I_MPI_SHM_LMT=shm ; mpiexec -np ',char(processes),' ./',char(executable_file),' ',char(k),' ',char(dataset_cols)];
end 
[status,cmdout] = system(command); 

cmdout   
    
fprintf("End of test\n");

