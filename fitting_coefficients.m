function [a1,a2,delta_0,delta_1] = fitting_coefficients( SETTING )

% FITTING_COEFFICIENTS:
% ---------------------
%
%   Returns a set of fitting coefficients for the problem described in:
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
%   The function returns the value of the dimensionless fitting parameters. It either
%   returns the coefficients from the manuscript (see below), or allows the input of the  
%   coefficients manually (see below). Additionally, the user can uncomment the lines
%   at the bottom of the function script ('fitting_coefficients.m') in order to add new 
%   sets of fitting coefficients.
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
%   SETTING           : Vector of characters indicating the chosen option for the fitting
%                       coefficients. The available options are:
%
%                       -'FROM_MANUSCRIPT': Values from the manuscript of the article in 
%                                           the description above.
%
%                       -'CUSTOM_INPUT': The fitting coefficients will be inputted 
%                                        manually by the user via command window.
%                   
%                       If 'SETTING' is not specified, it is set to 'FROM_MANUSCRIPT' by 
%                       default.
%
%
% OUTPUTS:
% --------
%
%   a1                : Coefficient a1 in equation (4.29) of the article, dimensionless.
%
%   a2                : Coefficient a2 in equation (4.29) of the article, dimensionless.
%
%   delta_0           : Coefficient delta_0i in equations (3.20-3.22) of the article,
%                       dimensionless.
%
%   delta_1           : Coefficient delta_1i in equations (3.20-3.22) of the article,
%                       dimensionless.
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
if nargin < 1

    SETTING = 'FROM_MANUSCRIPT';

end

% Check for conflictive input of TYPE
if ( ~isa(SETTING, 'char') )
    
    error('fitting_coefficients: SETTING (1st input) must be a character vector')
    
end

%----------------------------------------------------------------------------------------%
%-------------------------------SELECT LIQUID PARAMETERS---------------------------------%
%----------------------------------------------------------------------------------------%

switch lower(SETTING)
    
    case 'from_manuscript' % Parameters for water at 296 K
        
        a1      = 2.30;
        a2      = 0.319;
        delta_0 = 1.68;
        delta_1 = 0.0528;
          
%     case 'ADD_YOUR_OWN_1' % Uncomment these lines to add a custom set of parameters for
%                           % a different set of fitting parameters to this database
%         
%         a1          = %%%%%%%;
%         a2          = %%%%%%%;
%         delta_0     = %%%%%%%;
%         delta_1     = %%%%%%%;

%     case 'ADD_YOUR_OWN_2' % Uncomment these lines to add a custom set of parameters for
%                           % a different set of fitting parameters to this database
%         
%         a1          = %%%%%%%;
%         a2          = %%%%%%%;
%         delta_0     = %%%%%%%;
%         delta_1     = %%%%%%%;

%     case 'ADD_YOUR_OWN_3' % Uncomment these lines to add a custom set of parameters for
%                           % a different set of fitting parameters to this database
%         
%         a1          = %%%%%%%;
%         a2          = %%%%%%%;
%         delta_0     = %%%%%%%;
%         delta_1     = %%%%%%%;

    case 'custom_input' % Custom input of the parameters via command window
        
        disp('%------------------------------------------------------------------------%')
        disp('Please introduce the following custom fitting coefficients:')
        a1            = input('1) a1: ');
        a2            = input('2) a2: ');
        delta_0       = input('3) delta_0: ');
        delta_1       = input('3) delta_1: ');
        disp('%------------------------------------------------------------------------%')                
        
    otherwise
        
        error(['fitting_coefficients: Unknown option. Type "help fitting_coefficients"'...
                                     ' for more information.'])
                                
end

%----------------------------------------------------------------------------------------%
%----------------------------------CHECK OF PARAMETERS-----------------------------------%
%----------------------------------------------------------------------------------------%

% Check for conflicitive parameters (in case of manual input or user addition)
if ( ~isa(a1, 'double') || size(a1, 1) ~= 1 || size(a1, 2) ~= 1 || isnan(a1) ||        ...
     ~isfinite(a1) || ~isreal(a1) || a1<=0 )
 
    error('surfactant_properties: Inputted a1 must be a real and positive scalar')
    
end

if ( ~isa(a2, 'double') || size(a2, 1) ~= 1 || size(a2, 2) ~= 1 || isnan(a2) ||        ...
     ~isfinite(a2) || ~isreal(a2) || a2<=0 )
 
    error('surfactant_properties: Inputted a2 must be a real and positive scalar')
    
end

if ( ~isa(delta_0, 'double') || size(delta_0, 1) ~= 1 || size(delta_0, 2) ~= 1 ||      ...
     isnan(delta_0) || ~isfinite(delta_0) || ~isreal(delta_0) || delta_0<=0 )
 
    error('surfactant_properties: Inputted delta_0 must be a real and positive scalar')
    
end

if ( ~isa(delta_1, 'double') || size(delta_1, 1) ~= 1 || size(delta_1, 2) ~= 1 ||      ...
     isnan(delta_1) || ~isfinite(delta_1) || ~isreal(delta_1) || delta_1<=0 )
 
    error('surfactant_properties: Inputted delta_1 must be a real and positive scalar')
    
end

end