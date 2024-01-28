% the program calculate the thickness of each film in multilayers system
%  xinit--initial guess
%  tol--tolerance for termination
%  max_feval--maximum number of function evaluations 
%  myfunction--objective function. (In myfunction(x), x is a ROW vector) 
function d_fitting_list=fitting_regresion_analysis_simplex(myfunction,xinit,tol,max_feval)
x0=xinit(:)';  %row vector.
myfunction = fcnchk(myfunction);
dim=max(size(x0)); % dimension of the problem

% set up adaptive parameters
alpha=1; 
beta=1+2/dim; 
gamma=0.75-0.5/dim; 
delta=1-1/dim;
% initial simplex
scalefactor = min(max(max(abs(x0)),1),10);
D0=eye(dim);
D0(dim+1,:)=(1-sqrt(dim+1))/dim*ones(1,dim);
for i=1:dim+1
    X(i,:)=x0+ scalefactor*D0(i,:);
    FX(i)=feval(myfunction,X(i,:));
end
ct=dim+1;
[FX,I]=sort(FX);
X=X(I,:);
% start iteration
 while max(max(abs(X(2:dim+1,:)-X(1:dim,:)))) >= scalefactor*tol 
   if ct>max_feval
        break;
   end   
    M=mean(X(1:dim,:));
    xref=(1+alpha)*M- alpha*X(dim+1,:);
    Fref=feval(myfunction,xref);
    ct=ct+1;
    if Fref<FX(1)
        % expansion
        xexp=(1+alpha*beta)*M-alpha*beta*X(dim+1,:);
        Fexp=feval(myfunction,xexp);
        ct=ct+1;
        if Fexp < Fref
         X(dim+1,:)=xexp;
         FX(dim+1)=Fexp;
        else
         X(dim+1,:)=xref;
         FX(dim+1)=Fref;
        end
    else
        if Fref<FX(dim)
            % accept reflection point
            X(dim+1,:)=xref;
            FX(dim+1)=Fref;
        else 
            if Fref<FX(dim+1)
             % Outside contraction
                xoc=(1+alpha*gamma)*M-alpha*gamma*X(dim+1,:);
                Foc=feval(myfunction,xoc);
                ct=ct+1;           
              if Foc<=Fref
                X(dim+1,:)=xoc;
                FX(dim+1)=Foc;
              else
                % shrink
                for i=2:dim+1
                    X(i,:)=X(1,:)+ delta*(X(i,:)-X(1,:));
                    FX(i)=feval(myfunction,X(i,:));
                end
                ct=ct+dim;
              end
            else
              %inside contraction
              xic=(1-gamma)*M+gamma*X(dim+1,:);
              Fic=feval(myfunction,xic);
              ct=ct+1;
            if Fic<FX(dim+1)
                X(dim+1,:)=xic;
                FX(dim+1)=Fic;
            else
                % shrink
                for i=2:dim+1
                    X(i,:)=X(1,:)+ delta*(X(i,:)-X(1,:));
                    FX(i)=feval(myfunction,X(i,:));
                end
                ct=ct+dim;
            end
            end
        end
    end
    [FX,I]=sort(FX);
    X=X(I,:);
end
d_fitting_list=X(1,:);
% fmin=FX(1);
end
