%%%%%modified henon map
%       Plots semi-stable values of
%       x(n+1) = r*x(n)*(1-x(n)) as r increases to 4.
%
% FARSAN F J, LBS CENTRE FOR SCIENCE AND TECHNOLOGY, INDIA, 2019

clear
scale = 1000; % determines the level of rounding
maxpoints = 200; % determines maximum values to plot
N = 50; % number of "r" values to simulate
a = 0.0; % starting value of "r"
b = 12.0; % final value of "r"... anything higher diverges.
rs = linspace(a,b,N); % vector of "r" values
M = 500; % number of iterations of henon equation

% Loop through the "r" values
for j = 1:length(rs)
    
    r=rs(j); % get current "r"
    x=zeros(M,1);
     y=zeros(M,1);               % allocate memory
    x(1) = 0.5; % initial condition (can be anything from 0 to 1)
    y(1)=0.07;
    for i = 2:M, % iterate
        x(i)=1-r*cos(x(i-1))+0.3*y(i);
        y(i)=x(i);
    end
    % only save those unique, semi-stable values6
    out{j} = unique(round(scale*x(end-maxpoints:end)));
end

% Rearrange cell array into a large n-by-2 vector for plotting
data = [];
for k = 1:length(rs)
    n = length(out{k});
    data = [data;  rs(k)*ones(n,1),out{k}];
end

% Plot the data
figure(97);clf

 plot(data(:,1),data(:,2)/scale,'k.','markersize',0.0001);
% set(h,'markersize',0.1)
% axis tight
% set(gca,'units','normalized','position',[0 0 1 1])
% set(gcf,'color','white')
%axis on