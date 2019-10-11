function [U,T_TOTAL,T_ASSEMBLY,T_LINSOLVE] = compute_coeff_clean_case( L, phi, N )

% COMPUTE_COEFF_CLEAN_CASE:
% -------------------------
%
%   Computes the array of N approximate Fourier coefficients 
%   U = [E, d_1, d_2, ..., d_{N-1}] in the deviation streamfunction of the fluid flow 
%   problem described in:
% 
%       "A theory for the slip and drag of superhydrophobic surfaces with surfactant"
%
%                                    ---Authors:---
%
%                                   Julien R. Landel
%                                Francois J. Peaudecerf
%                               Fernando Temprano-Coleto
%                                    Frederic Gibou
%                                 Raymond E. Goldstein
%                                 Paolo Luzzatto-Fegiz
%
%   The function computes the coefficients for the clean case, in which the boundary
%   condition at the air gap is assumed to be purely no-shear and free of any Marangoni 
%   stress due to surfactants (gamma_Ma=0). Since the problem with surfactants assumes a 
%   spatially constant Marangoni shear and is linear in gamma_Ma, the coefficients for the
%   case with surfactant can easily be computed from the coefficients for the clean case 
%   multiplying every coefficient by (1-gamma_Ma).
%
%   The Fourier coefficients must be computed numerically by truncating the 
%   infinite series due to the mixed boundary conditions in the same segment of the
%   boundary. The linear system
%                   
%                                 A_{m,n}*U_{n} = B_{m}
%
%   is solved using the MATLAB built-in solver 'linsolve', which employs Gaussian 
%   elimination with partial pivoting. A is an NxN dense, non-singular matrix and B is a 
%   Nx1 right-hand-side (RHS) vector. The values of the entries for A and B, as well as 
%   more details about the problem, can be found in Section 4 of the manuscript.
%
%   We refer the interested reader to the manuscript for more information about the 
%   problem.
%
%   For the most updated version of the scripts, see: 
%
%                       https://github.com/feslab/shs-models-2d
%
%   WARNING: For a fast computation when a large number of terms is used, this code 
%            creates the matrix for the linear system via vectorization and matrix 
%            operations and without loops. This function have been tested in MATLAB 2018a,
%            but earlier versions might return errors due to different validity of the 
%            syntax for vector opeations used here.
%
%
% INPUTS:
% -------
%
%   L               : Nondimensional length of the problem domain, obtained by normalizing
%                     the dimensional length L^hat by the half-height of the domain h^hat.
%
%   phi             : Nondimensional gas fraction, which relates the length of the air gap
%                     g to the domain length by g = phi*L. The value of phi must therefore
%                     satisfy 0 <= phi <= 1.
%
%   N               : Number of computed Fourier coefficients in the truncated series. If
%                     N is not specified, it is set by default depending on the value of 
%                     the gas fraction as follows:
%                         a) If    0 <= phi <= 0.01, then N = 15000.
%                         b) If 0.01 <  phi <= 0.1 , then N = 2500.
%                         c) If 0.99 <= phi <= 1   , then N = 2500.
%                         d) If  0.1 <  phi <  0.99, then N = 500.
%
%
% OUTPUTS:
% --------
%
%   U               : Column vector of length N with the resulting computed coefficients,
%                     such that U = [E, d_1, d_2, ..., d_(N-1)].'.
%
%   T_TOTAL         : Total measured elapsed run time of the function.
%
%   T_ASSEMBLY      : Measured time of computation and assembly of the matrix and RHS
%                     coefficients.
%
%   T_LINSOLVE      : Measured time of solution of the linear system.
%
%
% LICENSE:
% ----------
%
%   This script is released under an MIT license. For more information, see:
%
%                       https://github.com/feslab/shs-models-2d
%
%   Copyright 2019 (c) Julien R. Landel, Francois J. Peaudecerf, Fernando Temprano-Coleto,
%   Frederic Gibou, Raymond E. Goldstein, Paolo Luzzatto-Fegiz.
% 
%   Permission is hereby granted, free of charge, to any person obtaining a copy of this 
%   software and associated documentation files (the "Software"), to deal in the Software
%   without restriction, including without limitation the rights to use, copy, modify, 
%   merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to the following 
%   conditions:
%
%   The above copyright notice and this permission notice shall be included in all copies 
%   or substantial portions of the Software.
%
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
%   INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
%   PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
%   CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%   
%   Contact email: ftempranocoleto@ucsb.edu
%                
%   October 2019

%----------------------------------------------------------------------------------------%
%------------------------------------CHECK OF INPUTS-------------------------------------%
%----------------------------------------------------------------------------------------%

% Check for omitted inputs
if nargin < 2

    error(['compute_coeff_clean_case: Not enough input arguments. Use "help '          ...
                                     'compute_coeff_clean_case" for more information.'])

elseif nargin < 3
    
    if( norm(phi)<=1e-2 )
        N = 15000;  
    elseif( norm(phi)<=1e-1 || norm(phi)>=0.99 )
        N = 2500;
    else
        N = 500;
    end

end

% Check for conflictive input of L
if ( ~isa(L, 'double') || size(L, 1) ~= 1 || size(L, 2) ~= 1 || isnan(L) ||            ...
     ~isfinite(L) || ~isreal(L) || L<=0 )
    error('compute_coeff_clean_case: L (1st input) must be a real and positive scalar')
end

if ( L <= 1e-3 )
    warning(['compute_coeff_clean_case: L (1st input) is very small. Please refer to ' ...
             'the asymptotic limits in the manuscript for an accurate and more '       ...
             'efficient evaluation of the coefficients in this limit.'])
elseif ( L >= 1e3 )
    warning(['compute_coeff_clean_case: L (1st input) is very large. Please refer to ' ...
             'the asymptotic limits in the manuscript for an accurate and more '       ...
             'efficient evaluation of the coefficients in this limit.'])
end

% Check for conflictive input of phi
if ( ~isa(phi, 'double') || size(phi, 1) ~= 1 || size(phi, 2) ~= 1 || isnan(phi) ||    ...
     ~isfinite(phi) || ~isreal(phi) || phi<0 || phi>1 )
    error(['compute_coeff_clean_case: phi (2nd input) must be a real scalar between 0 '...
           'and 1'])
end

% Check for conflictive input of N
if ( ~isa(N,'double') || size(N, 1) ~= 1 || size(N, 2) ~= 1 || isnan(N) ||             ...
     ~isfinite(N) || ~isreal(N) || fix(N)~=N || N<1 )
    error('compute_coeff_clean_case: N (3rd input) must be a real and positive integer')
end

if( (N<15000 && phi> 0    && phi<=0.01) ||                                             ...
    (N<2500  && phi> 0.01 && phi<=0.1 ) ||                                             ...
    (N<2500  && phi>=0.99 && phi<=1   ) ||                                             ...
    (N<500   && phi> 0.1  && phi< 0.99)    )
    warning(['compute_coeff_clean_case: N (3rd input) is small given the gas fraction '...
             'chosen. The computation might be under-resolved. A recommended '...
             'approximate number of terms to guarantee a fully converged truncated '   ...
             'series is: ' newline                                                     ...
                 'a) If    0 <= phi <= 0.01, then N = 15000.' newline                  ...
                 'b) If 0.01 <  phi <= 0.1 , then N = 2500.'  newline                  ...
                 'c) If 0.99 <= phi <= 1   , then N = 2500.'  newline                  ...
                 'd) If  0.1 <  phi <  0.99, then N = 500.'])
end

%----------------------------------------------------------------------------------------%
%--------------------DEFINITION OF WAVENUMBER AND COEFFICIENTS---------------------------%
%----------------------------------------------------------------------------------------%

k = @(n)(2*n*pi/L);

alpha = @(n)((-1+(4*k(n).*exp(-2*k(n)))+exp(-4*k(n)))./(1+exp(-2*k(n))));

beta  = @(n)( -2*k(n).*(1-(8*k(n).*exp(-4*k(n)))-exp(-8*k(n)))./...
                          ((1+exp(-2*k(n))).*(1+(4*k(n).*exp(-2*k(n)))-exp(-4*k(n)))) );

%----------------------------------------------------------------------------------------%
%---------------------------ASSEMBLY OF MATRIX AND RHS-----------------------------------%
%----------------------------------------------------------------------------------------%

% Start timer
tic;

% Define arrays for indices
n = (1:N-1) ;
m = (1:N-1)';

% Initialize matrix and RHS
A = zeros(N,N);
B = zeros(N,1);

% Compute and assemble matrix element for m = 0 and n = 0
A(1,1) = 1 - (phi/2);

% Compute and assemble matrix elements for m = 0 and n > 0
A(1,2:end) = (beta(n)-alpha(n)).*sin(pi*n*phi)./(2*pi*n);

% Compute and assemble matrix elements for m > 0 and n = 0
A(2:end,1) = - sin(pi*m*phi)./(2*pi*m);
            
% Compute and assemble the matrix elements with m ~= n, where m > 0 and n > 0           
A(2:end,2:end) = (1/4/pi)*(beta(n)-alpha(n))                                           ...
                        .*(sin(pi*(m+n)*phi)./(m+n) + sin(pi*(m-n)*phi)./(m-n));

% Compute and assemble the diagonal matrix elements with m = n > 0                       
A(2+N:N+1:N^2) = (alpha(n)/4) + (beta(n)-alpha(n)).*(phi/4 + sin(2*pi*n*phi)./(8*pi*n));
                 
% Compute and assemble the first RHS element
B(1) = phi/2;

% Compute and assemble ther remaining elements in RHS
B(2:end) = sin(pi*m*phi)./(2*pi*m);

% Stop timer
T_ASSEMBLY = toc;

%----------------------------------------------------------------------------------------%
%--------------------------SOLUTION OF THE LINEAR SYSTEM---------------------------------%
%----------------------------------------------------------------------------------------%

% Start timer
tic;
% Specify matrix characteristics for the linear solver
opts.LT      = false; 
opts.UT      = false;
opts.UHESS   = false; 
opts.SYM     = false; 
opts.RECT    = false;
opts.TRANSA  = false;
% Solve linear system
U = linsolve(A,B,opts);
% Stop timer
T_LINSOLVE = toc;

% Compute total time elapsed
T_TOTAL = T_ASSEMBLY + T_LINSOLVE;

end