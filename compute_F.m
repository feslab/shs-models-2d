function [F, coeff] = compute_F( x, g, phi, N)

% COMPUTE_F:
% ----------
%
%   Computes the value of the normalized slip velocity at the plastron in the fluid flow 
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
%   The function first computes the set of N approximate Fourier coefficients ( named E, 
%   d_1, d_2, ..., d_{N-1} ) of the deviation streamfunction in the surfactant-free 
%   (gamma_Ma = 0) case by calling 'compute_coeff_clean_case'. Since in a general case
%   with surfactant the coefficients are linear in (1-gamma_Ma), the function then 
%   computes the following normalized slip velocity:
%
%       F(x,g,phi) = u_I(x)/(2*(1-gamma_Ma)) =
%                  = E + SUM_{from n=1}^{to n=N-1}[ (d_n/2)*alpha_n*cos(k_n*x) ] ,
%
%   where alpha_n is a Fourier-type coefficient that is independent of the position (see 
%   the manuscript).
%
%   The computation of the coefficients is not necessary in the trivial cases {phi=0, g=0}
%   and {phi=1, g=L}, as well as in the asymptotic limits g/phi << 0  and g/phi >> 0. In
%   those cases, the normalized slip velocity follows a closed form (see the manuscript 
%   for details).
%
%   Note that F(x,g,phi) is independent of gamma_Ma, and depends on phi and g through the
%   Fourier coefficients E and d_n and through alpha_n. In order to obtain the actual
%   slip velocity for a general case, F should simply be multiplied by 2*(1-gamma_Ma).
%
%   We refer the interested reader to the manuscript for more information about the 
%   problem.
%
%   For the most updated version of the scripts, see: 
%
%                       https://github.com/feslab/shs-models-2d
%
%
% INPUTS:
% -------
%
%   x               : Nondimensional horizontal coordinate, obtained by normalizing the
%                     dimensional coordinate by the half-height of the domain h^hat. The 
%                     origin of the coordinate system is at the center of the domain, so x
%                     must satisfy -L/2 <= x <= L/2, where L is the length of the problem
%                     domain. This variable can be provided as a constant or as an array
%                     of values, in which case the function will return the value of F for
%                     each value of x in the array.
%
%   g               : Nondimensional gap length of the problem domain, obtained by 
%                     normalizing the dimensional gap length g^hat by the half-height of 
%                     the domain h^hat.
%
%   phi             : Nondimensional gas fraction, which relates the length of the air gap
%                     g to the domain length by g = phi*L. The plastron is therefore the 
%                     region of the bottom boundary (y=-1) of the domain given by
%                     |x| <= phi*L/2. The value of phi must satisfy 0 <= phi <= 1.
%
%   N               : Number of computed Fourier coefficients in the truncated series. If
%                     N is not specified, it is set by default by the function
%                     'compute_coeff_clean_case', which solves the linear system to obtain
%                     the coefficients. See "help compute_coeff_clean_case" for more
%                     information.
%
%
% OUTPUTS:
% --------
%
%   F               : Value of the normalized slip velocity for the values of x introduced
%                     as an input to the function, and for the given g and phi. F is an
%                     array with the same dimensions of x.
%
%   coeff           : Column array of N approximate Fourier coefficients 
%                     [E, d_1, d_2, ..., d_{N-1}] obtained in the calculation of F for
%                     given values of g and phi. The coefficients are obtained calling the
%                     function 'compute_coeff_clean_case'. In the trivial cases 
%                     {phi=0, g=0}, {phi=1, g=L} and in the asymptotic limits g/phi << 0 
%                     and g/phi >> 0 only the value of the first coefficient is returned
%                     coeff = E.
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
if nargin < 3

    error(['compute_F: Not enough input arguments. Use "help compute_F" for more '     ...
                      'information.'])

end

% Check for conflictive input of x
if ( ~isa(x, 'double') || ((size(x, 1) ~= 1) && (size(x, 2) ~= 1)) || any(isnan(x)) || ...
     any(~isfinite(x)) || any(~isreal(x)) ||  (phi~=0 && g~=0 && max(x)>g/(2*phi)) ||  ...
     (phi~=0 && g~=0 && min(x)<-g/(2*phi)) )
    error(['compute_F: x (1st input) must be an array of real values between '         ...
                      '-g/(2*phi) and g/(2*phi).'])
end

% Check for conflictive input of g
if ( ~isa(g, 'double') || size(g, 1) ~= 1 || size(g, 2) ~= 1 || isnan(g) ||            ...
     ~isfinite(g) || ~isreal(g) || g<0 )
    error('compute_F: g (2nd input) must be a real and non-negative scalar')
end

% Check for conflicts when g or phi are zero
if ( g==0 && phi~=0 ) 
    warning(['compute_F: The input g=0 requires phi=0 as well by definition. Changing '...
                        'phi to phi=0...'])
    phi = 0;
elseif (g~=0 && phi==0)
    warning(['compute_F: The input phi=0 requires g=0 as well by definition. Changing '...
                        'g to g=0...'])
    g = 0;     
end

% NOTE: The inputs phi and N are checked by the function 'compute_coeff_clean_case' called
%       below. In the cases {phi=0, g=0} and {phi=1, g=L} there is no need to obtain the 
%       coefficients or compute any sum, and the input N doesn't intervene in the 
%       calculations and therefore phi and N are not checked for conflicts.

%----------------------------------------------------------------------------------------%
%------------------------------COMPUTATION OF THE FUNCTION-------------------------------%
%----------------------------------------------------------------------------------------%

if(phi==0 || g ==0) % With no gap, the normalized slip velocity is zero everywhere
    
    F     = zeros(size(x,1),size(x,2));
    coeff = 0;
    
elseif(phi==1) % With no wall, the normalized slip velocity is equal to 1 everywhere
    
    F     = ones(size(x,1),size(x,2));
    coeff = 1;
    
elseif(g/phi > 1e3) % In the asymptotic limit L>>0, apply equation (C8) from the article,
                    % with the appropriate normalization
    
    F = zeros(size(x,1),size(x,2));
    
    F(logical(abs(x)<g/2)) = 1/(4-3*phi);

    coeff = phi/(4-3*phi);   
    
elseif(g/phi < 1e-3) % In the asymptotic limit L<<0, apply equation (C20) from the 
                     % article, with the appropriate normalization
   
    F = zeros(size(x,1),size(x,2));

    F(logical(abs(x)<g/2)) = ((g/4/pi/phi)/(1+(g/4/pi/phi)*log(sec(pi*phi/2))))*       ...
                             acosh(cos(pi*phi*x(logical(abs(x)<g/2))/g)/cos(pi*phi/2));
    
    coeff = (g/4/pi/phi)*log(sec(pi*phi/2))/(1+(g/4/pi/phi)*log(sec(pi*phi/2))); 
    
else % In the general case;
    
    % Calculate approximate Fourier coefficients
    k = @(n)(2*n*pi*phi/g);
    alpha_hat = @(n)((-1+(4*k(n)*exp(-2*k(n)))+exp(-4*k(n)))/(1+exp(-2*k(n))));
    
    if nargin < 4
        
        [coeff,~,~,~] = compute_coeff_clean_case(g/phi, phi);
        N = length(coeff);
        
    else
        
        [coeff,~,~,~] = compute_coeff_clean_case(g/phi, phi, N);
    
    end

    E = coeff(1);
    d_hat = @(n)(coeff(n+1));
    
    % Compute the sum for each value of x
    F = zeros(size(x,1),size(x,2));
    idx = logical(abs(x)<g/2); 
    F(idx) = E;

    for n=1:N-1        
        
        F(idx) = F(idx) +  (d_hat(n)/2)*alpha_hat(n)*cos(k(n)*x(idx));   
    
    end
    
end

end