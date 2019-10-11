function [rho_hat,mu_hat,T_hat] = liquid_properties( TYPE )

% LIQUID_PROPERTIES:
% ------------------
%
%   Returns a set of values of properties of a given liquid for the problem described in:
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
%   The function returns the value of the properties in SI units. It has only one option 
%   for the liquid type (see below), and another option to input the liquid type manually 
%   through the command window (see below). Additionally, the user can uncomment the lines
%   at the bottom of the function script ('liquid_properties.m') in order to add new 
%   liquids with custom properties.
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
%   TYPE              : Vector of characters indicating the liquid type. The available 
%                       options are:
%
%                       -'FROM_MANUSCRIPT': Properties used as reference case for the
%                                           numerical simulations in the manuscript. These 
%                                           are the values that the liquid properties take
%                                           in each simulation, unless specified otherwise
%                                           in the table provided in the supplementary
%                                           material. See the table for a complete
%                                           characterization of every parameter in each
%                                           simulation.
%
%                       -'WATER_296K': Properties for water at 296 K.
%
%                       -'CUSTOM_INPUT': The properties of the liquid will be manually 
%                                        inputted by the user via command window. 
%                   
%                       If 'TYPE' is not specified, it is set to 'FROM_MANUSCRIPT' by 
%                       default.
%
%
% OUTPUTS:
% --------
%
%   rho_hat           : Mass density, in [kg/m^3].
%
%   mu_hat            : Dynamic viscosity, in [kg/m/s].
%
%   T_hat             : Absolute temperature, in [K].
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

    TYPE = 'FROM_MANUSCRIPT';

end

% Check for conflictive input of TYPE
if ( ~isa(TYPE, 'char') )
    
    error('liquid_properties: TYPE (1st input) must be a character vector.')
    
end

%----------------------------------------------------------------------------------------%
%-------------------------------SELECT LIQUID PARAMETERS---------------------------------%
%----------------------------------------------------------------------------------------%

switch lower(TYPE)
    
    case 'water_296k' % Parameters for water at 296 K
        
        rho_hat     = 1e3;
        mu_hat      = 8.9e-4;
        T_hat       = 296;
        
    case 'from_manuscript' % Parameters used by default for the parametric sweep in the 
                           % manuscript.
        
        rho_hat     = 1e3;
        mu_hat      = 1e-3;
        T_hat       = 293;
          
%     case 'ADD_YOUR_OWN_1' % Uncomment these lines to add a custom set of parameters for
%                           % a different type of liquid to this database
%         
%         rho_hat     = %%%%%%%;
%         mu_hat      = %%%%%%%;
%         T_hat       = %%%%%%%;

%     case 'ADD_YOUR_OWN_2' % Uncomment these lines to add a custom set of parameters for
%                           % a different type of liquid to this database
%         
%         rho_hat     = %%%%%%%;
%         mu_hat      = %%%%%%%;
%         T_hat       = %%%%%%%;

%     case 'ADD_YOUR_OWN_3' % Uncomment these lines to add a custom set of parameters for
%                           % a different type of liquid to this database
%         
%         rho_hat     = %%%%%%%;
%         mu_hat      = %%%%%%%;
%         T_hat       = %%%%%%%;

    case 'custom_input' % Custom input of the parameters via command window
        
        disp('%------------------------------------------------------------------------%')
        disp('Please introduce the following liquid parameters:')
        rho_hat       = input('1) Mass density, in [kg/m^3]: ');
        mu_hat        = input('2) Dynamic viscosity, in [kg/m/s]: ');
        T_hat         = input('3) Absolute temperature, in [K]: ');
        disp('%------------------------------------------------------------------------%')                
        
    otherwise
        
        error(['liquid_properties: Unknown liquid type. Type "help liquid_properties" '...
                                  'for more information.'])
                                
end

%----------------------------------------------------------------------------------------%
%----------------------------------CHECK OF PARAMETERS-----------------------------------%
%----------------------------------------------------------------------------------------%

% Check for conflicitive parameters (in case of manual input or user addition)
if ( ~isa(rho_hat, 'double') || size(rho_hat, 1) ~= 1 || size(rho_hat, 2) ~= 1 ||      ...
     isnan(rho_hat) || ~isfinite(rho_hat) || ~isreal(rho_hat) || rho_hat<=0 )
 
    error('surfactant_properties: Inputted rho_hat must be a real and positive scalar')
    
end

if ( ~isa(mu_hat, 'double') || size(mu_hat, 1) ~= 1 || size(mu_hat, 2) ~= 1 ||         ...
     isnan(mu_hat) || ~isfinite(mu_hat) || ~isreal(mu_hat) || mu_hat<=0 )
 
    error('surfactant_properties: Inputted mu_hat must be a real and positive scalar')
    
end

if ( ~isa(T_hat, 'double') || size(T_hat, 1) ~= 1 || size(T_hat, 2) ~= 1 ||            ...
     isnan(T_hat) || ~isfinite(T_hat) || ~isreal(T_hat) || T_hat<=0 )
 
    error('surfactant_properties: Inputted T_hat must be a real and positive scalar')
    
end

end