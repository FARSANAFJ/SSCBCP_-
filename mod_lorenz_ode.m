function dx=lotenz_utube(~,x,beta)
dx=[beta(1)*(x(2)-x(1));x(1)*(beta(2)-x(3))-x(2);x(1)*x(2)-beta(3)*x(3);
    beta(4)*beta(1)*(x(2)-x(1)+beta(5)*x(1))
];
end



%%%%(w) ?=k_p a(y-x)+k_i x 
%kp=beta(4)
%ki=beta(5)