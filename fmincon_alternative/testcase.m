x1=[0.8147;0.9058];
x2=[0.1270; 0.9134];
l=[0.6324; 0.0975]; % variables

lb=-1*ones(2,1);
ub=1*ones(2,1);

%%%%% trust region %%%%%
[newl,newq]=fmincon(@(l) testfmincon(l,x1,x2),l,[],...
    [],[],[],lb,ub,[],optimset('Algorithm','trust-region-reflective',...
    'GradObj','on','Hessian','on','TolX',1e-5,'MaxIter',20,'Display','iter'));
fprintf(1, '\n\n#### trust region #####\n');
display(newl);
display(newq);

%%%%% interior-point %%%%%
[newl, newq] = fmincon(@(l) testfmincon2(l,x1,x2),l,[],[],[],[],lb,ub,[],...
    optimset('Algorithm','interior-point','GradObj','on','Hessian', ...
    'user-supplied','HessFcn',@test_hessian,'TolX',1e-5,'MaxIter',20,'Display','iter'));
    
fprintf(1, '\n\n#### interior-point #####\n');
display(newl);
display(newq);