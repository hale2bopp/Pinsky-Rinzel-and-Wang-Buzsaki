clc
close all;
clear all;
global N I g E  adj mat num_diff
N = 10; f = 0.5; C = 0.5;
size_wb = 3; size_pr = 8;

tspan = [1 1000];
tic
num_diff = ceil(f*N)*size_wb+ceil((1-f)*N)*size_pr+N*(N-1);
% fraction of the neurons that are excited
stim = 1; %fraction of stimulated neurons
stimmag=1.0; %magnitude of stimulation
I=zeros(N,1); %stimulation current
k=randperm(N);kk=k(1:ceil(stim*N));
I(kk)=stimmag
I
%for iter = 1:10
%for C = 0:0.1:1
adj=0.5+0.5*sign(C-rand(N));
for i=1:N,adj(i,i)=0;end;

gsyn = 0.1;
g = gsyn*bsxfun(@rdivide,adj,(1+sum(adj,2)))
%Synaptic reversal potentials
E=5*ones(1,N); %excitatory
k=randperm(N);kk=k(1:ceil(f*N));
E(kk)=-5

[T,Y] = ode15s(@oscnetwork_pr_wb,tspan,zeros(1,num_diff));
size(Y)
list  = 1:N;
col=hsv(N);
grid on
figure(1);
count = 0;
for i = mat
    count =count+1;
    grid on
    if E(count)<0
        subplot(2,1,1)
        ax1 = plot(T,Y(:,i),'color',col(count,:));
        hold on

    elseif E(count)>0
        subplot(2,1,2) 
        ax2 = plot(T,Y(:,i),'color',col(count,:));
        hold on
    end

    xlabel('Time(s)');
    ylabel('Voltage (mV)');
    title('Voltage vs Time-for N cells');       
end



toc
