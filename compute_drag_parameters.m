function [lambda_e, DR, gamma_Ma, lambda_e_hat, gamma_Ma_hat] = ...
                    compute_drag_parameters(h_hat, g_hat, L_hat, U_hat, c_0_hat, varargin)

% COMPUTE_DRAG_PARAMETERS:
% ------------------------
%   Computes the effective slip length, drag reduction and Marangoni shear of a 
%   flow over a surfactant-contaminated superhydrophobic surface using the model in: 
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
%   The function calculates the characteristic dimensionless groups of the problem from
%   dimensional inputs from the user. It then calculates the Marangoni shear using the
%   scaling law (4.29) in the manuscript, and the associated slip length and drag
%   reduction.
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
%   h_hat           : Half height of the channel, in [m].
%
%   g_hat           : Length of the gap (plastron), in [m].
%
%   L_hat           : Length of the domain (gap plus ridge), in [m].
%
%   U_hat           : Characteristic velocity of the flow, in [m/s]. It is twice the value
%                     of the maximum velocity for a plane Poiseuille flow with the same
%                     mean pressure gradient.
%
%   c_0_hat         : Background bulk surfactant concentration, in [Mol/m^3] (or, 
%                     equivalently, in [mMol/l]).
%
%   OPTIONS         : Following all the inputs above, it is possible to set different
%                     options for the calculation of the drag parameters. These options
%                     must be declared following the syntax:
%
%                     compute_drag_parameters(... , OPTION_NAME_1, OPTION_SETTING_1 , 
%                                                   OPTION_NAME_2, OPTION_SETTING_2 ,
%                                                         ...    ,       ...        ,
%                                                   OPTION_NAME_n, OPTION_SETTING_n    )
%
%                     OPTION_NAME_i is a vector of characters (text) that specifies the
%                     option name (for example, 'SURFACTANT_TYPE' enables the user to 
%                     switch properties between different surfactants). OPTION_SETTING_i
%                     denotes the specific setting for the option given by OPTION_NAME_i
%                     (see examples below).
%
%                     The possible options for this function are:
%
%                     - 'NUM_TERMS': Specifies the number of terms used in the truncated
%                     series for the calculation of the fluid flow (see the manuscript).
%                     It is followed by the actual value N of the number of terms.  For
%                     example, typing <'NUM_TERMS', 1200> among the function inputs will
%                     set the number of terms in the truncation to 1200. If NUM_TERMS is 
%                     not specified, it is set automatically as follows:
%                         a) If    0 <= phi <= 0.01, then N = 15000.
%                         b) If 0.01 <  phi <= 0.1 , then N = 2500.
%                         c) If 0.99 <= phi <= 1   , then N = 2500.
%                         d) If  0.1 <  phi <  0.99, then N = 500.
%                     Where phi is the gas fraction g_hat/L_hat. Type 
%                     "help compute_coeff_clean_case" for details about the number of 
%                     terms N.
%
%                     - 'SURFACTANT_TYPE': Specifies the type of surfactant present in the
%                     fluid flow by calling the external function 'surfactant_properties'.
%                     For example, typing <'SURFACTANT_TYPE', 'SDS'> among the inputs of
%                     this function will set the surfactant properties to those of sodium 
%                     dodecyl sulfate. If SURFACTANT_TYPE is not specified, it is set
%                     automatically to 'FROM_MANUSCRIPT', i.e. the values used as a 
%                     reference case for the numerical simulations in the manuscript. For 
%                     more information about available surfactant options and on how to 
%                     setup a custom set of surfactant properties, type "help 
%                     surfactant_properties".
%
%                     - 'LIQUID_TYPE': Specifies the type of liquid and its temperature
%                     by calling the external function 'liquid_properties'. For example, 
%                     typing <'LIQUID_TYPE', 'WATER_296K'> among the inputs of this 
%                     function will set the liquid properties to those of water at 296 K.
%                     If LIQUID_TYPE is not specified, it is set automatically to
%                     'FROM_MANUSCRIPT', i.e. the values used as a reference case for the 
%                     numerical simulations in the manuscript. For more information about 
%                     available liquid options and on how to input a custom set of liquid 
%                     properties, type "help liquid_properties".
%
%                     - 'FITTING PARAM': Specifies the values of the fitting parameters
%                     used in the scaling model (see the manuscript) by calling the 
%                     external function 'fitting_coefficients'. For example, typing
%                     <'FITTING_PARAM', 'FROM_MANUSCRIPT'> will set the values to those
%                     reported in the text. If FITTING_PARAM is not specified, it is set
%                     automatically to 'FROM_MANUSCRIPT'. For more information about 
%                     options for the fitting parameters and how to input a custom set of 
%                     coefficients, type "help fitting_coefficients".
%
%                     - 'BL_SCALING': Specifies the scaling law for the diffusive boundary
%                     layer thickness for the bulk surfactant, by calling the external
%                     function 'bl_scaling'. For example, typing <'BL_SCALING', 'LEVEQUE'>
%                     will set a boundary layer thickness following the Leveque regime
%                     (see the manuscript). Note that the dependence of the Marangoni 
%                     shear, slip length and drag reduction on the boundary layer 
%                     thickness is weak, so changing this functional form will likely not
%                     affect the results considerably. IF BL_SCALING is not specified, it
%                     is set automatically to 'LEVEQUE'. For more information about the
%                     boundary layer scaling settings, type "help bl_scaling".
%
%   A non-exhaustive list of examples of input lists for this function is the following:
%       
%    1) compute_drag_parameters(1e-4, 4e-4, 5e-4, 1e-5, 5e-3)
%
%    2) compute_drag_parameters(1e-4, 2e-4, 3e-4, 5e-5, 1e-3, 'BL_SCALING', 'CLASSIC')
%
%    3) compute_drag_parameters(2e-4, 1e-4, 2e-4, 1e-4, 1e-2, 'NUM_TERMS', 500,       ...
%                                                           'LIQUID_TYPE', 'CUSTOM_INPUT')
%
%    4) compute_drag_parameters(2e-4, 1e-4, 2e-4, 1e-4, 1e-2, 'SURFACTANT_TYPE', ...
%                                                               'WEAK', 'NUM_TERMS', 6000)
%
%
% OUTPUTS:
% --------
%
%   lambda_e        : Nondimensional effective slip length.
%
%   DR              : Nondimensional drag reduction.
%
%   gamma_Ma        : Nondimensional Marangoni shear rate.
%
%   lambda_e_hat    : Dimensional effective slip length in [m], with h_hat the constant of
%                     normalization.
%
%   gamma_Ma_hat    : Dimensional Marangoni shear rate in [kg*m/(s^2)], with 
%                     mu_hat*U_hat/h_hat the constant of normalization.
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
if nargin < 5

    error(['compute_drag_parameters: Not enough input arguments. Use "help '           ...
                                    'compute_drag_parameters" for more information.']) 
         
end

% Check for conflictive input of h_hat
if ( ~isa(h_hat, 'double') || size(h_hat, 1) ~= 1 || size(h_hat, 2) ~= 1 ||            ...
      isnan(h_hat) || ~isfinite(h_hat) || ~isreal(h_hat) || h_hat<=0 )
    error(['compute_drag_parameters: h_hat (1st input) must be a real and positive '   ...
                                     'scalar'])
end

% Check for conflictive input of g_hat
if ( ~isa(g_hat, 'double') || size(g_hat, 1) ~= 1 || size(g_hat, 2) ~= 1 ||            ...
      isnan(g_hat) || ~isfinite(g_hat) || ~isreal(g_hat) || g_hat<=0 )
    error(['compute_drag_parameters: g_hat (2nd input) must be a real and positive '   ...
                                     'scalar'])
end

% Check for conflictive input of L_hat
if ( ~isa(L_hat, 'double') || size(L_hat, 1) ~= 1 || size(L_hat, 2) ~= 1 ||            ...
      isnan(L_hat) || ~isfinite(L_hat) || ~isreal(L_hat) || L_hat<=0 )
    error(['compute_drag_parameters: L_hat (3rd input) must be a real and positive '   ...
                                     'scalar'])
end

% Check for conflictive input of U_hat
if ( ~isa(U_hat, 'double') || size(U_hat, 1) ~= 1 || size(U_hat, 2) ~= 1 ||            ...
      isnan(U_hat) || ~isfinite(U_hat) || ~isreal(U_hat) || U_hat<=0 )
    error(['compute_drag_parameters: U_hat (4th input) must be a real and positive '   ...
                                     'scalar'])
end

% Check for conflictive input of c_0_hat
if ( ~isa(c_0_hat, 'double') || size(c_0_hat, 1) ~= 1 || size(c_0_hat, 2) ~= 1 ||      ...
      isnan(c_0_hat) || ~isfinite(c_0_hat) || ~isreal(c_0_hat) || c_0_hat<=0 )
    error(['compute_drag_parameters: c_0_hat (5th input) must be a real and positive ' ...
                                     'scalar'])
end

%----------------------------------------------------------------------------------------%
%--------------------------------SELECTION OF OPTIONS------------------------------------%
%----------------------------------------------------------------------------------------%

% Select default options
[D_hat,D_I_hat,k_a_hat,k_d_hat,Gamma_m_hat,n_sigma_hat] = surfactant_properties('SDS');
[rho_hat,mu_hat,T_hat] = liquid_properties('WATER_296K');
[a1,a2,delta_0,delta_1] = fitting_coefficients('FROM_MANUSCRIPT');
BL_TYPE = 'LEVEQUE';
N = [];

flag_num_terms       = false;
flag_surfactant_type = false;
flag_liquid_type     = false;
flag_fitting_param   = false;
flag_bl_scaling      = false;

flag_warning_num_terms       = false;
flag_warning_surfactant_type = false;
flag_warning_liquid_type     = false;
flag_warning_fitting_param   = false;
flag_warning_bl_scaling      = false;

% Parse through the inputted options for custom options
while ~isempty(varargin)
    switch lower(varargin{1})
        
        case 'num_terms'
            
            if(flag_num_terms==false)
                N = varargin{2};
                flag_num_terms = true;
            elseif(flag_warning_num_terms==false)
                warning(['compute_drag_parameters: The option NUM_TERMS has been '     ...
                                                  'specified more than once. Only the '...
                                                  'first inputted value will be '      ...
                                                  'considered']) 
                flag_warning_num_terms = true;
            end
        
        case 'surfactant_type'
            
            if(flag_surfactant_type==false)                
                [D_hat,D_I_hat,k_a_hat,k_d_hat,Gamma_m_hat,n_sigma_hat] =              ...
                                                     surfactant_properties( varargin{2} );
                flag_surfactant_type = true;
            elseif(flag_warning_surfactant_type==false)
                warning(['compute_drag_parameters: The option SURFACTANT_TYPE has been'...
                                                  ' specified more than once. Only the'...
                                                  ' first inputted option will be '    ...
                                                  'considered'])              
                flag_warning_surfactant_type = true;                              
            end
                                                   
        case 'liquid_type'
            
            if(flag_liquid_type==false) 
                [rho_hat,mu_hat,T_hat] = liquid_properties( varargin{2} );
                flag_liquid_type = true;
            elseif(flag_warning_liquid_type==false)
                warning(['compute_drag_parameters: The option LIQUID_TYPE has been'    ...
                                                  ' specified more than once. Only the'...
                                                  ' first inputted option will be '    ...
                                                  'considered'])
                flag_warning_liquid_type = true;
            end
            
        case 'fitting_param'            
            
            if(flag_fitting_param==false)
                [a1,a2,delta_0,delta_1] = fitting_coefficients( varargin{2} );
                flag_fitting_param = true;
            elseif(flag_warning_fitting_param==false)
                warning(['compute_drag_parameters: The option FITTING_PARAM has been'  ...
                                                  ' specified more than once. Only the'...
                                                  ' first inputted option will be '    ...
                                                  'considered'])
                flag_warning_fitting_param = true;
            end
            
        case 'bl_scaling'
            
            if(flag_bl_scaling==false)
                BL_TYPE = varargin{2};
                flag_bl_scaling = true;
            elseif(flag_warning_bl_scaling==false)
                warning(['compute_drag_parameters: The option BL_TYPE has been'        ...
                                                  ' specified more than once. Only the'...
                                                  ' first inputted option will be '    ...
                                                  'considered'])
                flag_warning_bl_scaling = true;
            end
                
        otherwise
            
            error('compute_drag_parameters: Unknown input option')
        
    end
    
    varargin(1:2) = [];
end

%----------------------------------------------------------------------------------------%
%--------------CALCULATION OF THE DIMENSIONLESS NUMBERS AND PARAMETERS-------------------%
%----------------------------------------------------------------------------------------%

% Define the ideal gas constant in SI units
R_hat = 8.3144598;

% Calculation of the problem dimensionless groups
g     = g_hat/h_hat;
phi   = g_hat/L_hat;
k     = k_a_hat*c_0_hat/k_d_hat;
Pe    = h_hat*U_hat/D_hat;
Pe_I  = h_hat*U_hat/D_I_hat;
Bi    = k_d_hat*h_hat/U_hat;
chi   = k_d_hat*h_hat/k_a_hat/Gamma_m_hat;
Ma    = n_sigma_hat*R_hat*T_hat*Gamma_m_hat/mu_hat/U_hat; 

% Calculation of the boundary layer thickness
delta = bl_scaling(Pe,g,delta_0,delta_1,BL_TYPE);

% Computation of the approximate Fourier coefficients for the streamfunction
if(isempty(N))
    [F_0,coeff] = compute_F(0,g,phi);
else
    [F_0,coeff] = compute_F(0,g,phi,N);
end
E_0 = coeff(1);

% Calculation of other dimensionless parameters
Re     = rho_hat*h_hat*U_hat/mu_hat;
Pe_I_g = g*F_0*Pe_I;
K_I_g  = g*Bi*(1+k)/F_0;
D_I_g  = g*chi*(1+k)/sqrt(F_0*g*Pe);

%----------------------------------------------------------------------------------------%
%----------------------------WARNINGS OF MODEL BREAKDOWN---------------------------------%
%----------------------------------------------------------------------------------------%

if( Re > 1 && Re <= 1000)
    
    warning(['compute_drag_parameters: The Reynolds number is larger than 1, breaking '...
                                      'the Stokes flow assumption. Expect mild errors '... 
                                      'in all calculated quantities.'])
end

if( Re > 1000 ) 
    
    warning(['compute_drag_parameters: The Reynolds number is larger than 1000, '      ...
                                      'potentially breaking the laminar flow '         ...
                                      'assumption in a plane channel flow. Expect '    ... 
                                      'large errors in all calculated quantities.'])
end
                                  
if( k > 1 )
    
    warning(['compute_drag_parameters: The nondimensional concentration k is larger '  ...
                                      'than 1, breaking the dilute regime assumption.' ...
                                      'Expect a large underprediction of the slip '    ...
                                      'length and the drag reduction.'])
end
    
if( Pe_I_g >= 10 && ( K_I_g <= 50 || D_I_g <= 50 ) )

    warning(['compute_drag_parameters: The dimensionless groups suggest that the '     ...
                                      'surfactant distribution might be in a full '    ...
                                      'stagnant cap regime, potentially breaking the ' ...
                                      'uniform Marangoni shear approximation. Expect ' ...
                                      'a mild underprediction of the slip length and ' ...
                                      'the drag reduction.'])
                                  
elseif( Pe_I_g >= 1e3 && ( K_I_g <= 0.5 || D_I_g <= 0.5 ) )

    warning(['compute_drag_parameters: The dimensionless groups suggest that the '     ...
                                      'surfactant distribution might be in a partial ' ...
                                      'stagnant cap regime, breaking the uniform '     ...
                                      'Marangoni shear approximation. Expect a large ' ...
                                      'underprediction of the slip length and the '    ...
                                      'drag reduction.'])
    
end

%----------------------------------------------------------------------------------------%
%--------------------------CALCULATION OF DRAG PARAMETERS--------------------------------%
%----------------------------------------------------------------------------------------%

gamma_Ma = a1*k*Ma*F_0 / ((1/Pe_I) + (a2*(g^2)*Bi/(1+(Bi*Pe*delta/chi))) + (a1*k*Ma*F_0));

lambda_e = 2*(1-gamma_Ma)*E_0/(1-((1-gamma_Ma)*E_0));

DR       = 1 - 1/(1+(3*lambda_e)/(lambda_e+2))^2;

lambda_e_hat = lambda_e*h_hat;

gamma_Ma_hat = gamma_Ma*mu_hat*U_hat/h_hat;

end