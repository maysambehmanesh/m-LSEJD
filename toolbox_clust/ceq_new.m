function [c,ceq,CG,CGeq] = ceq_new(x,Z)

    % CONSTR inputs a vector x (   = [A(:);B(:)]   )
    % and returns the constraint vector [A'A - I; B'B -I] = 0
    % and the Jacobian CGe of this vector of nonlinear equations
    % The nonlinear inequalities c and CG are [] (no inequalities)

        ZK2=size(x,1);
        ZK=ZK2/2;
        K=round(ZK/Z);

        a = x(1:ZK);
        b = x((ZK+1):end);

        % reshaping
        A = reshape( a, K, Z );
        B = reshape( b, K, Z );

        I = eye( Z );

        % constraints in matrix form
        res_vecA = A'*A - I;
        res_vecB = B'*B - I;

        % vectorization of constraints
        ceq = [ res_vecA(:); res_vecB(:) ];

        % there is NO inequality constraint
        c = [];

        % gradients

        if nargout > 2
            CG  = [];
            CGeq = grad_constr(A,B);
        end
        
        function CGeq=grad_constr(A,B)

        % GRAD_CONSTR computes the Jacobian of the constraints A'A-I=0, B'B-I=0
        % Input: matrices A, B
        % Output: matrix of partial derivatives (Jacobian) to the components of A and B
        %          
            [K,Z]=size(A);
            [K,Z]=size(B);
    
            CGeq = zeros(2*Z*Z,2*Z*K);
    
            CGeqA = (speye(Z^2)+vecperm(Z,Z))*kron(eye(Z),A');
            CGeqB = (speye(Z^2)+vecperm(Z,Z))*kron(eye(Z),B');
    
            CGeq(1:Z*Z,1:Z*K)=CGeqA;
            CGeq(Z*Z+1:2*Z*Z,Z*K+1:2*Z*K)=CGeqB;

            CGeq=CGeq';
        end     
        
end
