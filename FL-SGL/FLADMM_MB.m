function [theta] = FLADMM_MB(XT,YT,s1,s2,s3,opts)

    Print = opts.Print;
    Times = opts.Times;
    Dimen = opts.Dimen;
    MaxIter = opts.MaxIter;
    rho = opts.rho;
    sigma = opts.sigma;
    rho_1 = 3*rho;
    Tol = opts.Tol;
   
    for i = 1:Times
        for j = 1:Times
            D(i,j) = exp(-((i-j)^2)/(2*sigma^2));
        end
    end
    DT = sum(D) - 1;
    for i = 1:Times
        D(:,i) = -D(:,i)/DT(i);
        D(i,i) = 1;
    end 
    DD = -D + eye(Times);
    
    bv = [];
    I = eye(Dimen);    
    for t = 1:Times
        A = XT{1,t};
        y = YT{1,t};
        bv = cat(2,bv,A'*y);   % X'*y
        CC{1,t} = A'*A + (2*rho+rho_1)*I;
        Ch_R{1,t} = dpotrf_(CC{1,t});   % Ch_R = U'*U
    end
    
    mu = zeros(Dimen,Times);
    u = zeros(Dimen,Times);
    v = zeros(Dimen,Times);
    q = zeros(Dimen,Times);
    theta = zeros(Dimen,Times);
    gamma = zeros(Dimen,Times);
    pi = zeros(Dimen,Times);
    theta_p = zeros(Dimen,Times);
    
    if Print
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\n','iter','funv', 'Tol tv', 'td', 'Tol tz');
    end
      
    for iter = 1:MaxIter
        % x-update
        TTtheta = rho*theta_p*DD; 
        thetat = bv-mu-u+rho*q+TTtheta+rho*gamma+rho_1*theta_p;
        for t = 1:Times
            theta(:,t) = dpotrs_(Ch_R{1,t},thetat(:,t));
        end
        
        Ttheta = theta*D;
        gamma = pi + Ttheta + (u - v)/rho;
        
        q = theta + mu/rho;
        q = shrinkage(q,s1/rho);
        q = shrinkage_21(q,s2/rho);
        
        pi = gamma + v/rho;
        pi = shrinkage(pi,s3/rho);
        
        mu = mu + rho*(theta - q);
        u = u + rho*(Ttheta - gamma);
        v = v + rho*(gamma - pi);
        
        funv(iter) = norm(theta-q);
        td = norm(theta-theta_p);
        
        if iter > 1
            tv = max(1,funv(iter-1));
            tz = max(1,td);
            if Print
                fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n',iter,funv(iter), Tol*tv,td,Tol*tz);
             end
            if funv(iter) < Tol*tv && td < Tol*tz
                break;
            end
        end
        
        theta_p = theta;
        
    end
    theta = q;
             
function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
    
function w = shrinkage_21(W,lambda_3)

    for i = 1:size(W,1)
        w_1 = W(i,:);
        nm = norm(w_1, 2);
        if nm == 0
            w_2 = zeros(size(w_1));
        else
            w_2 = max(nm - lambda_3, 0)/nm * w_1;
        end
        w(i,:) = w_2;
    end

    