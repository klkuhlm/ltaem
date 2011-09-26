function f = dhksinv(F, t, M, alpha, tol, varargin)
%
% F: name of function that calculates complex vector values of Laplace transform
%    F(p,param1,param2,...) for complex vector argument p and optional parameters
%    param1,param2,... contained in varargin.
%
% t: real vector values of t at which inversion is to be calculated
%
% f: real vector values of inversion corresponding to input t vector
%
% alpha: real parameter
%
% tol: real parameter
%
% varagin: optional list of parameter values to be passed to function that
%          evaluates the Laplace transform
%
%example:
%
%
%
% de Hoog, F. R., Knight, J. H., and Stokes, A. N. (1982) An improved
%    method for numerical inversion of Laplace transforms. S.I.A.M. J. Sci.
%    and Stat. Comput., 3, 357-366.
%
% This is a May, 2011 vectorised MATLAB translation of the FORTRAN program
% written by  John Knight at CSIRO, Canberra, Australia in August 1978.
%
% contact:  john@johnknight.com.au      http://www.johnknight.com.au
%
st = size(t);
onet = complex(ones(st));
zerot = complex(zeros(st));
T = 2*max(max(t));
% M = 20;
M2 = 2*M;
gamma = alpha - log(tol)/(2*T);
kk = (0:2*M); % 2*M+1 terms
p = gamma + 1i*kk*pi/T;   % equally spaced points on Bromwich contour, 1i = sqrt(-1) in MATLAB
a = feval(F,p,varargin);  % evaluate complex vector Laplace transform
%                          for complex vector argument p and optional parameters
a(1) = 0.5*a(1); % 0 order coefficient, index of a is shifted to start at 1
%
% construct quotient difference table recursively one Q or E column
% at a time and extract continued fraction coefficients d
%
d = complex(ones(1,M2)); % preallocate
L1 = M2; L2 = L1 + 1;
Q = a(2:L2)./ a(1:L1) ;% initialise
L2 = L1; L1 = L2 - 1;
E = Q(2:L2) - Q(1:L1); % initialise
d(1) = -Q(1);    % continued fraction coefficient
d(2) = -E(1);    % continued fraction coefficient
I1 = 2;  % index of d starts at zero
for K = 2:M  %  indices of Q and E are shifted to start at 1 instead of 0
   L2 = L1; L1 = L2 - 1; I1 = I1 + 1;
   Q = Q(2:L2).* E(2:L2)./ E(1:L1); % next Q column in table
   d(I1) = -Q(1);                   % continued fraction coefficient
   L2 = L1; L1 = L2 - 1; I1 = I1 + 1;
   E = E(2:L2) + Q(2:L2) - Q(1:L1); % next E column in table
   d(I1) = -E(1);                   % continued fraction coefficient
end
%
% evaluate continued fraction form of diagonal Pade approximation to power series in z
%
z = cos(pi*t/T) + 1i*sin(pi*t/T); % complex argument of power series
AOLD = a(1)*onet; BOLD = onet; % initialise recurrence, d(0) = a(1)
AOLD2 = zerot; BOLD2 = onet;
for n = 1:M2
   if n == M2
      %H2M = 0.5*(onet+(d(M2-1)-d(M2))*z); % index of d starts at zero
      % H2MSQ = H2M.*H2M;
       term1 = d(M2)*z; %
       %term = -H2M.*(onet - sqrt(onet + term1./H2MSQ)); % remainder on last recurrence
       rem = term1.*(onet - term.*(onet - term1 - term)); % simplified remainder
       term = rem;
   else
       term = d(n)*z; % index of d starts at zero
   end
   A = AOLD + term.*AOLD2;
   AOLD2 = AOLD;
   AOLD = A;
   B = BOLD + term.*BOLD2;
   BOLD2 = BOLD;
   BOLD = B;
end
f = real(A./B).*exp(gamma*t)/T;

