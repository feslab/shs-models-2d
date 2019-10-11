function [D_hat,D_I_hat,k_a_hat,k_d_hat,Gamma_m_hat,n_sigma_hat] = ...
                                                             surfactant_properties( TYPE )

% SURFACTANT_PROPERTIES:
% ----------------------
%
%   Returns a set of values of properties of a given surfactant for the problem described 
%   in:
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
%   The function returns the value of the properties in SI units. It has three options
%   for the surfactant type (see below) and the option to input the surfactant type 
%   manually through the command window (see below). Additionally, the user can uncomment
%   the lines at the bottom of the function script ('surfactant_properties.m') in order to
%   add new surfactants with custom properties.
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
%   TYPE            : Vector of characters indicating the surfactant type. The available 
%                     options are:
%                       
%                     -'FROM_MANUSCRIPT': Properties used as reference case for the
%                                         numerical simulations in the manuscript. These 
%                                         are the values that the surfactant properties 
%                                         take in each simulation, unless specified 
%                                         otherwise in the table provided in the 
%                                         supplementary material. See the table for a
%                                         complete characterization of every parameter in
%                                         each simulation.
%
%                     -'SDS': Properties for sodium dodecyl sulfate, as reported in
%                             'CHANG & FRANSES, 1995' DOI:10.1016/0927-7757(94)03061-4 and
%                             'PROSSER & FRANSES, 2001' DOI:10.1016/S0927-7757(00)00706-8.
%
%                     -'SDS_INSOL': Properties for a hypothetical 'insoluble SDS', as 
%                                   described in the manuscript.
%
%                     -'STRONG': Properties for a hypothetical 'strong' surfactant, as
%                                described in 'PEAUDECERF ET AL., 2017' 
%                                DOI:10.1073/pnas.1702469114.
%
%                     -'WEAK': Properties for a hypothetical 'weak' surfactant, as
%                              described in 'PEAUDECERF ET AL., 2017' 
%                              DOI:10.1073/pnas.1702469114.
%
%                     -'CUSTOM_INPUT': The properties of the surfactant will be manually 
%                                      inputted by the user via command window. 
%                   
%                     If 'TYPE' is not specified, it is set to 'FROM_MANUSCRIPT' by 
%                     default.
%
%
% OUTPUTS:
% --------
%
%   D_hat           : Diffusion coefficient of the bulk surfactant, in [m^2/s].
%
%   D_I_hat         : Diffusion coefficient of the interface surfactant, in [m^2/s].
%
%   k_a_hat         : Adsorption rate constant, in [m^3/mol/s].
%
%   k_d_hat         : Desorption rate constant, in [s^(-1)].
%
%   Gamma_m_hat     : Maximum packing concentration of the interface surfactant, in    
%                     [mol/m^2].
%
%   n_sigma_hat     : Surfactant style constant, dimensionless.
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
    
    error('surfactant_properties: TYPE (1st input) must be a character vector.')
    
end

%----------------------------------------------------------------------------------------%
%-----------------------------SELECT SURFACTANT PARAMETERS-------------------------------%
%----------------------------------------------------------------------------------------%

switch lower(TYPE)
    
    case 'sds' % Parameters for sodium dodecyl sulfate, as reported in: 
               % 'CHANG & FRANSES (1995)' https://doi.org/10.1016/0927-7757(94)03061-4 
               % 'PROSSER & FRANSES (2001)' https://doi.org/10.1016/S0927-7757(00)00706-8
        
        D_hat       = 7e-10;
        D_I_hat     = 7e-10;
        k_a_hat     = 89.5;
        k_d_hat     = 500;
        Gamma_m_hat = 3.9e-6;
        n_sigma_hat = 2;    
        
    case 'strong' % Parameters for a hypothetical 'strong' surfactant, as described in
                  % 'PEAUDECERF ET AL. (2017)' https://doi.org/10.1073/pnas.1702469114
        
        D_hat       = 1e-11;
        D_I_hat     = 1e-11;
        k_a_hat     = 1e6;
        k_d_hat     = 1;
        Gamma_m_hat = 1e-5;
        n_sigma_hat = 2;
        
    case 'weak' % Parameters for a hypothetical 'weak' surfactant, as described in
                % 'PEAUDECERF ET AL. (2017)' https://doi.org/10.1073/pnas.1702469114
        
        D_hat       = 1e-9;
        D_I_hat     = 1e-9;
        k_a_hat     = 1e-1;
        k_d_hat     = 100;
        Gamma_m_hat = 1e-6;
        n_sigma_hat = 2;
        
    case 'from_manuscript' % Parameters used by default for the parametric sweep in the 
                           % manuscript.
        
        D_hat       = 1e-10;
        D_I_hat     = 1e-10;
        k_a_hat     = 10;
        k_d_hat     = 10;
        Gamma_m_hat = 5e-6;
        n_sigma_hat = 2; 
        
    case 'sds_insol' % Parameters for a hypothetical 'insoluble SDS', as described in the 
                     % manuscript.
        
        D_hat       = 7e-10;
        D_I_hat     = 7e-10;
        k_a_hat     = 89.5;
        k_d_hat     = 1;
        Gamma_m_hat = 3.9e-6;
        n_sigma_hat = 2;      
          
%     case 'ADD_YOUR_OWN_1' % Uncomment these lines to add a custom set of parameters for
%                           % a different type of surfactant to this database
%         
%         D_hat       = %%%%%%%;
%         D_I_hat     = %%%%%%%;
%         k_a_hat     = %%%%%%%;
%         k_d_hat     = %%%%%%%;
%         Gamma_m_hat = %%%%%%%;  
%         n_sigma_hat = %%%%%%%;

%     case 'ADD_YOUR_OWN_2' % Uncomment these lines to add a custom set of parameters for
%                           % a different type of surfactant to this database
%         
%         D_hat       = %%%%%%%;
%         D_I_hat     = %%%%%%%;
%         k_a_hat     = %%%%%%%;
%         k_d_hat     = %%%%%%%;
%         Gamma_m_hat = %%%%%%%;  
%         n_sigma_hat = %%%%%%%;

%     case 'ADD_YOUR_OWN_3' % Uncomment these lines to add a custom set of parameters for
%                           % a different type of surfactant to this database
%         
%         D_hat       = %%%%%%%;
%         D_I_hat     = %%%%%%%;
%         k_a_hat     = %%%%%%%;
%         k_d_hat     = %%%%%%%;
%         Gamma_m_hat = %%%%%%%;  
%         n_sigma_hat = %%%%%%%;

    case 'custom_input' % Custom input of the parameters via command window
        
        disp('%------------------------------------------------------------------------%')
        disp('Please introduce the following surfactant transport parameters:')
        D_hat       = input('1) Bulk surfactant diffusivity D_hat, in [m^2/s]: ');
        D_I_hat     = input('2) Interface surfactant diffusivity D_I_hat, in [m^2/s]: ');
        k_a_hat     = input('3) Adsorption rate constant k_a_hat, in [m^3/mol/s]: ');
        k_d_hat     = input('4) Desorption rate constant k_d_hat, in [s^(-1)]: ');
        Gamma_m_hat = input(['5) Maximum packing concentration Gamma_m_hat, in '       ...
                             '[mol/m^2]: ']);
        n_sigma_hat = input('6) Surfactant style constant (dimensionless): ');
        disp('%------------------------------------------------------------------------%')                
        
    otherwise
        
        error(['surfactant_properties: Unknown surfactant type. Type '                 ...
                                      '"help surfactant_properties" for more '         ...
                                      'information.'])
                                
end

%----------------------------------------------------------------------------------------%
%----------------------------------CHECK OF PARAMETERS-----------------------------------%
%----------------------------------------------------------------------------------------%

% Check for conflicitive parameters (in case of manual input or user addition)
if ( ~isa(D_hat, 'double') || size(D_hat, 1) ~= 1 || size(D_hat, 2) ~= 1 ||    ...
     isnan(D_hat) || ~isfinite(D_hat) || ~isreal(D_hat) || D_hat<=0 )
 
    error('surfactant_properties: Inputted D_hat must be a real and positive scalar')
    
end

if ( ~isa(D_I_hat, 'double') || size(D_I_hat, 1) ~= 1 || size(D_I_hat, 2) ~= 1 ||      ...
     isnan(D_I_hat) || ~isfinite(D_I_hat) || ~isreal(D_I_hat) || D_I_hat<=0 )
 
    error('surfactant_properties: Inputted D_I_hat must be a real and positive scalar')
    
end

if ( ~isa(k_a_hat, 'double') || size(k_a_hat, 1) ~= 1 || size(k_a_hat, 2) ~= 1 ||      ...
     isnan(k_a_hat) || ~isfinite(k_a_hat) || ~isreal(k_a_hat) || k_a_hat<=0 )
 
    error('surfactant_properties: Inputted k_a_hat must be a real and positive scalar')
    
end

if ( ~isa(k_d_hat, 'double') || size(k_d_hat, 1) ~= 1 || size(k_d_hat, 2) ~= 1 ||      ...
     isnan(k_d_hat) || ~isfinite(k_d_hat) || ~isreal(k_d_hat) || k_d_hat<=0 )
 
    error('surfactant_properties: Inputted k_d_hat must be a real and positive scalar')
    
end

if ( ~isa(Gamma_m_hat, 'double') || size(Gamma_m_hat, 1) ~= 1 ||                       ...
     size(Gamma_m_hat, 2) ~= 1 || isnan(Gamma_m_hat) || ~isfinite(Gamma_m_hat) ||      ...
     ~isreal(Gamma_m_hat) || Gamma_m_hat<=0 )
 
    error(['surfactant_properties: Inputted Gamma_m_hat must be a real and positive '  ...
           'scalar'])
       
end

if ( ~isa(n_sigma_hat, 'double') || size(n_sigma_hat, 1) ~= 1 ||                       ...
     size(n_sigma_hat, 2) ~= 1 || isnan(n_sigma_hat) || ~isfinite(n_sigma_hat) ||      ...
     ~isreal(n_sigma_hat) || n_sigma_hat<=0 )
 
    error(['surfactant_properties: Inputted n_sigma_hat must be a real and positive '  ...
           'scalar'])
       
end

end