function delta = bl_scaling( Pe, g, delta_0, delta_1, TYPE)

% BL_SCALING:
% -----------
%
%   Returns the scaling law for the bulk surfactant boundary layer thickness in the
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
%   The function returns the value of the boundary layer (BL) thickness based on scaling
%   arguments presented in Appendix B of the manuscript. The input 'TYPE' allows to switch
%   from a 'classic' (-1/2 power) scaling law to a 'Leveque regime' (-1/3 power) scaling 
%   law.
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
%   Pe                : Nondimesional Peclet number Pe = h^hat*U^hat/D^hat.
%
%   g                 : Nondimensional gap length of the problem domain, obtained by 
%                       normalizing the dimensional gap length g^hat by the half-height of 
%                       the domain h^hat.
%
%   delta_0           : First fitting parameter in the scaling law (see Appendix B of the 
%                       manuscript).
%
%   delta_1           : Second fitting parameter in the scaling law (see Appendix B of the 
%                       manuscript).
%
%   TYPE              : Vector of characters indicating the type of boundary layer 
%                       scaling. The possible options are:
%                       
%                       -'LEVEQUE': The scaling of the BL thickness is assumed to follow a
%                                   Leveque regime law of the form:
%                                   
%                                   delta ~ g*delta_0*(1+delta_1*g^2*Pe)^(-1/3)
%
%                       -'CLASSIC': The scaling of the BL thickness is assumed to follow a
%                                   'classic' scaling law of the form:
%
%                                   delta ~ g*delta_0*(1+delta_1*g^2*Pe)^(-1/2)  for g<~1
%
%                                   delta ~ g*delta_0*(1+delta_1*g*Pe)^(-1/2)    for g>~1
%
%                       If 'TYPE' is not specified, it is set to 'LEVEQUE' by default.
%
%
% OUTPUTS:
% --------
%
%   delta             : Nondimensional boundary layer thickness.
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
if nargin < 4

    error(['bl_scaling: Not enough input arguments. Use "help bl_scaling" for more '   ...
                       'information.'])

elseif nargin < 5

    TYPE = 'LEVEQUE';

end

% Check for conflictive input of Pe
if ( ~isa(Pe, 'double') || size(Pe, 1) ~= 1 || size(Pe, 2) ~= 1 || isnan(Pe) ||        ...
     ~isfinite(Pe) || ~isreal(Pe) || Pe<0 )
    error('bl_scaling: Pe (1st input) must be a real and non-negative scalar')
end

% Check for conflictive input of g
if ( ~isa(g, 'double') || size(g, 1) ~= 1 || size(g, 2) ~= 1 || isnan(g) ||            ...
     ~isfinite(g) || ~isreal(g) || g<0 )
    error('bl_scaling: g (2nd input) must be a real and non-negative scalar')
end

% Check for conflictive input of delta_0
if ( ~isa(delta_0, 'double') || size(delta_0, 1) ~= 1 || size(delta_0, 2) ~= 1 ||      ...
      isnan(delta_0) || ~isfinite(delta_0) || ~isreal(delta_0) || delta_0<0 )
    error('bl_scaling: delta_0 (3rd input) must be a real and non-negative scalar')
end

% Check for conflictive input of delta_1
if ( ~isa(delta_1, 'double') || size(delta_1, 1) ~= 1 || size(delta_1, 2) ~= 1 ||      ...
      isnan(delta_1) || ~isfinite(delta_1) || ~isreal(delta_1) || delta_1<0 )
    error('bl_scaling: delta_1 (4th input) must be a real and non-negative scalar')
end


%----------------------------------------------------------------------------------------%
%--------------------------SELECTION OF BOUNDARY LAYER SCALING---------------------------%
%----------------------------------------------------------------------------------------%

switch lower(TYPE)
    case 'leveque'
        
        delta = g*delta_0*(1+delta_1*(g^2)*Pe)^(-1/3);
        
    case 'classic'
        
        if(g>1)
            
            delta = g*delta_0*(1+delta_1*g*Pe)^(-1/2);
            
        else
            
            delta = g*delta_0*(1+delta_1*(g^2)*Pe)^(-1/2);
            
        end
        
    otherwise
        
        error(['bl_scaling: Unknown input option. Type "help bl_scaling" for more '    ...
                           'information'])
        
end

end