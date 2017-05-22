clc;  clear all;


datadir = 'Hopkins155';
seqs = dir(datadir);
seq3 = seqs(3:end);

%% Load the Data
data = struct('X',{},'name',{},'ids',{});
for i=1:length(seq3)
    fname = seq3(i).name;
    fdir = [datadir '/' fname];
    if isdir(fdir)
        datai = load([fdir '/' fname '_truth.mat']);
        id = length(data)+1;
        data(id).ids = datai.s;
        data(id).name = lower(fname);
        X = reshape(permute(datai.x(1:2,:,:),[1 3 2]),2*datai.frames,datai.points);
        data(id).X = [X;ones(1,size(X,2))];
    end
end


p = 4;  % which task


data1 = data(p).X;
label1 = data(p).ids;
label1 = find(label1==1);
n1 = length(label1);
data1 = data1(:, label1);
data2 = data(p).X;
label2 = data(p).ids;
label2 = find(label2==2);
n2 = length(label2);
data2 = data2(:, label2);
data3 = data(p).X;
label3 = data(p).ids;
label3 = find(label3==3);
n3 = length(label3);
data3 = data3(:, label3);

%% Process Data
M = [data1, data2, data3];  % delete 3
label = [ones(1, n1), 2*ones(1, n2), 3*ones(1, n3)]';  % delete 3
[m, n] = size(M);
%[U, S, V] = svd(M);
%L = U(:, 1:8)*S(1:8, 1:8)*V(:, 1:8)';
L = M;

%% Life-long Matrix Completion
[L_hat, success] = mc_sp(L, ceil(0.95*m), 8);
precision = norm(L_hat-L, 'fro')/norm(L, 'fro')




