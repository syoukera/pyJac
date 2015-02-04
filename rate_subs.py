"""Module for writing species/reaction rate subroutines.

This is kept separate from Jacobian creation module in order
to create only the rate subroutines if desired.
"""

# Python 2 compatibility
from __future__ import division
from __future__ import print_function

# Standard libraries
import sys
import math

# Local imports
import chem_utilities as chem
import mech_interpret as mech
import utils

def rxn_rate_const(A, b, E):
    """Returns line with reaction rate calculation (after = sign).
    
    Notes
    -----
    Form of the reaction rate constant (from, e.g., Lu and Law [1]_):
    .. math::
        :nowrap:
        k_f = \begin{cases}
        A & \text{if } \beta = 0 \text{ and } T_a = 0 \\
        \exp \left( \log A + \beta \log T \right) & 
        \text{if } \beta \neq 0 \text{ and } \text{if } T_a = 0 \\
        \exp \left( \log A + \beta \log T - T_a / T \right)	& 
        \text{if } \beta \neq 0 \text{ and } T_a \neq 0 \\
        \exp \left( \log A - T_a / T \right) 
        & \text{if } \beta = 0 \text{ and } T_a \neq 0 \\
        A \prod^b T	& \text{if } T_a = 0 \text{ and } 
        b \in \mathbb{Z} \text{ (integers) }
        \end{cases}
    
    .. [1] TF Lu and CK Law, "Toward accommodating realistic fuel chemistry
       in large-scale computations," Progress in Energy and Combustion 
       Science, vol. 35, pp. 192-215, 2009. doi:10.1016/j.pecs.2008.10.002
    
    Parameters
    ----------
    A : float
        Arrhenius pre-exponential coefficient
    b : float
        Arrhenius temperature exponent
    E : float
        Arrhenius activation energy
    
    Returns
    -------
    line : str
        String with expression for reaction rate.
    
    """
    
    line = ''
    logA = math.log(A)
    
    if not E:
        # E = 0
        if not b:
            # b = 0
            line += str(A)
        else:
            # b != 0
            if isinstance(b, int):
                line += str(A)
                for i in range(b):
                    line += ' * T'
            else:
                line += 'exp({:.8e}'.format(logA)
                if b > 0:
                    line += ' + ' + str(b)
                else:
                    line += ' - ' + str(abs(b))
                line += ' * logT)'
    else:
        # E != 0
        if not b:
            # b = 0
            line += 'exp({:.8e}'.format(logA) + ' - ({:.8e} / T))'.format(E)
        else:
            # b!= 0
            line += 'exp({:.8e}'.format(logA)
            if b > 0:
                line += ' + ' + str(b)
            else:
                line += ' - ' + str(abs(b))
            line += ' * logT - ({:.8e} / T))'.format(E)
    
    return line


def write_rxn_rates(path, lang, specs, reacs):
    """Write reaction rate subroutine.
    
    Includes conditionals for reversible reactions.
    
    Parameters
    ----------
    path : str
        Path to build directory for file.
    lang : {'c', 'cuda', 'fortran', 'matlab'}
        Programming language.
    specs : list of SpecInfo
        List of species in the mechanism.
    reacs : list of ReacInfo
        ist of reactions in the mechanism.
    
    Returns
    _______
    None
    
    """
    
    num_s = len(specs)
    num_r = len(reacs)
    rev_reacs = [rxn for rxn in reacs if rxn.rev]
    num_rev = len(rev_reacs)
    
    pdep_reacs = []
    for reac in reacs:
        if reac.thd or reac.pdep:
            # add reaction index to list
            pdep_reacs.append(reacs.index(reac))
    
    # first write header file
    if lang == 'c':
        file = open(path + 'rates.h', 'w')
        file.write('#ifndef RATES_HEAD\n'
                   '#define RATES_HEAD\n'
                   '\n'
                   '#include "header.h"\n\n'
                   )
        if rev_reacs:
            file.write('void eval_rxn_rates (const Real, const Real*, '
                       'Real*, Real*);\n'
                       'void eval_spec_rates (const Real*, const Real*, '
                       'const Real*, Real*);\n'
                       )
        else:
            file.write('void eval_rxn_rates (const Real, const Real*, '
                       'Real*);\n'
                       'void eval_spec_rates (const Real*, const Real*, '
                       'Real*);\n'
                       )
        
        if pdep_reacs:
            file.write('void get_rxn_pres_mod (const Real, const Real, '
                       'const Real*, Real*);\n'
                       )
        
        file.write('\n'
                   '#endif\n'
                   )
        file.close()
    elif lang == 'cuda':
        file = open(path + 'rates.cuh', 'w')
        file.write('#ifndef RATES_HEAD\n'
                   '#define RATES_HEAD\n'
                   '\n'
                   '#include "header.h"\n'
                   '\n'
                   )
        if rev_reacs:
            file.write('__device__ void eval_rxn_rates (const Real, '
                       'const Real*, Real*, Real*);\n'
                       '__device__ void eval_spec_rates (const Real*, '
                       'const Real*, const Real*, Real*);\n'
                       )
        else:
            file.write('__device__ void eval_rxn_rates (const Real, const '
                       'Real*, Real*);\n'
                       '__device__ void eval_spec_rates (const Real*, const '
                       'Real*, Real*);\n'
                       )
        
        if pdep_reacs:
            file.write('__device__ void get_rxn_pres_mod (const Real, const '
                       'Real, const Real*, Real*);\n'
                       )
        
        file.write('\n')
        file.write('#endif\n')
        file.close()
    
    filename = 'rxn_rates' + utils.file_ext[lang]
    file = open(path + filename, 'w')
    
    if lang in ['c', 'cuda']:
        file.write('#include <math.h>\n'
                   '#include "header.h"\n'
                   '\n'
                   )
    
    line = ''
    if lang == 'cuda': line = '__device__ '
    
    if lang in ['c', 'cuda']:
        if rev_reacs:
            line += ('void eval_rxn_rates (const Real T, const Real * C, '
                     'Real * fwd_rxn_rates, Real * rev_rxn_rates) {\n'
                     )
        else:
            line += ('void eval_rxn_rates (const Real T, const Real * C, '
                     'Real * fwd_rxn_rates) {\n'
                     )
    elif lang == 'fortran':
        if rev_reacs:
            line += ('subroutine eval_rxn_rates (T, C, fwd_rxn_rates, '
                     'rev_rxn_rates)\n\n'
                     )
        else:
            line += 'subroutine eval_rxn_rates (T, C, fwd_rxn_rates)\n\n'
        
        # fortran needs type declarations
        line += ('  implicit none\n'
                 '  double precision, intent(in) :: T, C({})\n'.format(num_s)
                 )
        if rev_reacs:
            line += ('  double precision, intent(out) :: '
                     'fwd_rxn_rates({}), '.format(num_r) + 
                     'rev_rxn_rates({})\n'.format(num_rev)
                     )
        else:
            line += ('  double precision, intent(out) :: '
                     'fwd_rxn_rates({})\n'.format(num_r)
                     )
        line += ('  \n'
                 '  double precision :: logT\n'
                 )
        if rev_reacs and any(rxn.rev_par != [] for rxn in rev_reacs):
                line += '  double precision :: kf, Kc\n'
        line += '\n'
    elif lang == 'matlab':
        if rev_reacs:
            line += ('function [fwd_rxn_rates, rev_rxn_rates] = '
                     'eval_rxn_rates (T, C)\n\n'
                     '  fwd_rxn_rates = zeros({},1);\n'.format(num_r) + 
                     '  rev_rxn_rates = fwd_rxn_rates;\n'
                     )
        else:
            line += ('function fwd_rxn_rates = eval_rxn_rates (T, C)\n\n'
                     '  fwd_rxn_rates = zeros({},1);\n'.format(num_r)
                     )
    file.write(line)
    
    pre = '  '
    if lang == 'c':
        pre += 'Real '
    elif lang == 'cuda':
        pre += 'register Real '
    line = (pre + 'logT = log(T)' + 
            utils.line_end[lang]
            )
    file.write(line)
    file.write('\n')
    
    if rev_reacs and any(rxn.rev_par == [] for rxn in rev_reacs):
        if lang == 'c':
            file.write('  Real kf;\n'
                       '  Real Kc;\n'
                       )
        elif lang == 'cuda':
            file.write('  register Real kf;\n'
                       '  register Real Kc;\n'
                       )
    
    file.write('\n')
    
    for rxn in reacs:
        
        # if reversible, save forward rate constant for use
        if rxn.rev and not rxn.rev_par:
            line = ('  kf = ' + rxn_rate_const(rxn.A, rxn.b, rxn.E) + 
                    utils.line_end[lang]
                    )
            file.write(line)
        
        line = '  fwd_rxn_rates'
        if lang in ['c', 'cuda']:
            line += '[{}] = '.format(reacs.index(rxn))
        elif lang in ['fortran', 'matlab']:
            line += '({}) = '.format(reacs.index(rxn) + 1)
        
        # reactants
        for sp in rxn.reac:
            isp = next(i for i in xrange(len(specs)) if specs[i].name == sp)
            nu = rxn.reac_nu[rxn.reac.index(sp)]
            
            # check if stoichiometric coefficient is real or integer
            if isinstance(nu, float):
                if lang in ['c', 'cuda']:
                    line += 'pow(C[{}], {}) * '.format(isp, nu)
                elif lang in ['fortran', 'matlab']:
                    line += 'pow(C({}), {}) * '.format(isp + 1, nu)
            else:
                # integer, so just use multiplication
                for i in range(nu):
                    if lang in ['c', 'cuda']:
                        line += 'C[{}] * '.format(isp)
                    elif lang in ['fortran', 'matlab']:
                        line += 'C({}) * '.format(isp + 1)
        
        # Rate constant: print if not reversible, or reversible but 
        # with explicit reverse parameters.
        if not rxn.rev or rxn.rev_par:
            line += rxn_rate_const(rxn.A, rxn.b, rxn.E)
        else:
            line += 'kf'
        
        line += utils.line_end[lang]
        file.write(line)
        
        if rxn.rev:
            
            if not rxn.rev_par:
                
                line = '  Kc = 0.0' + utils.line_end[lang]
                file.write(line)
                
                # sum of stoichiometric coefficients
                sum_nu = 0
                
                # go through product species
                for prod_sp in rxn.prod:
                    isp = rxn.prod.index(prod_sp)
                    
                    # check if species also in reactants
                    if prod_sp in rxn.reac:
                        isp2 = rxn.reac.index(prod_sp)
                        nu = rxn.prod_nu[isp] - rxn.reac_nu[isp2]
                    else:
                        nu = rxn.prod_nu[isp]
                    
                    # Skip species with zero overall 
                    # stoichiometric coefficient.
                    if (nu == 0):
                        continue
                    
                    sum_nu += nu
                    
                    # get species object
                    sp = next((sp for sp in specs if 
                               sp.name == prod_sp), None)
                    if not sp:
                        print('Error: species ' + prod_sp + ' in reaction '
                              '{} not found.\n'.format(reacs.index(rxn))
                              )
                        sys.exit()
                    
                    # need temperature conditional for equilibrium constants
                    line = '  if (T <= {:})'.format(sp.Trange[1])
                    if lang in ['c', 'cuda']:
                        line += ' {\n'
                    elif lang == 'fortran':
                        line += ' then\n'
                    elif lang == 'matlab':
                        line += '\n'
                    file.write(line)
                    
                    line = '    Kc '
                    if lang in ['c', 'cuda']:
                        if nu < 0:
                            line += '-= '
                        elif nu > 0:
                            line += '+= '
                    elif lang in ['fortran', 'matlab']:
                        if nu < 0:
                            line += '= Kc - '
                        elif nu > 0:
                            line += '= Kc + '
                    line += '{:.2f} * '.format(abs(nu))
                    line += ('({:.8e} - '.format(sp.lo[6]) + 
                             '{:.8e} + '.format(sp.lo[0]) + 
                             '{:.8e} * '.format(sp.lo[0] - 1.0) + 
                             'logT + T * ('
                             '{:.8e} + T * ('.format(sp.lo[1] / 2.0) + 
                             '{:.8e} + T * ('.format(sp.lo[2] / 6.0) + 
                             '{:.8e} + '.format(sp.lo[3] / 12.0) + 
                             '{:.8e} * T))) - '.format(sp.lo[4] / 20.0) + 
                             '{:.8e} / T)'.format(sp.lo[5]) + 
                             utils.line_end[lang]
                             )
                    file.write(line)
                    
                    if lang in ['c', 'cuda']:
                        file.write('  } else {\n')
                    elif lang in ['fortran', 'matlab']:
                        file.write('  else\n')
                    
                    line = '    Kc '
                    if lang in ['c', 'cuda']:
                        if nu < 0:
                            line += '-= '
                        elif nu > 0:
                            line += '+= '
                    elif lang in ['fortran', 'matlab']:
                        if nu < 0:
                            line += '= Kc - '
                        elif nu > 0:
                            line += '= Kc + '
                    line += '{:.2f} * '.format(abs(nu))
                    line += ('({:.8e} - '.format(sp.hi[6]) + 
                             '{:.8e} + '.format(sp.hi[0]) + 
                             '{:.8e} * '.format(sp.hi[0] - 1.0) + 
                             'logT + T * ('
                             '{:.8e} + T * ('.format(sp.hi[1] / 2.0) + 
                             '{:.8e} + T * ('.format(sp.hi[2] / 6.0) + 
                             '{:.8e} + '.format(sp.hi[3] / 12.0) + 
                             '{:.8e} * T))) - '.format(sp.hi[4] / 20.0) + 
                             '{:.8e} / T)'.format(sp.hi[5]) + 
                             utils.line_end[lang]
                             )
                    file.write(line)
                    
                    if lang in ['c', 'cuda']:
                        file.write('  }\n\n')
                    elif lang == 'fortran':
                        file.write('  end if\n\n')
                    elif lang == 'matlab':
                        file.write('  end\n\n')
                
                # now loop through reactants
                for reac_sp in rxn.reac:
                    isp = rxn.reac.index(reac_sp)
                    
                    # Check if species also in products; 
                    # if so, already considered).
                    if reac_sp in rxn.prod: continue
                    
                    nu = rxn.reac_nu[isp]
                    sum_nu -= nu
                    
                    # get species object
                    sp = next((sp for sp in specs if sp.name == reac_sp), 
                              None)
                    if not sp:
                        print('Error: species ' + reac_sp + ' in reaction '
                              '{} not found.\n'.format(reacs.index(rxn))
                              )
                        sys.exit()
                        
                    # need temperature conditional 
                    # for equilibrium constants
                    line = '  if (T <= {:})'.format(sp.Trange[1])
                    if lang in ['c', 'cuda']:
                        line += ' {\n'
                    elif lang == 'fortran':
                        line += ' then\n'
                    elif lang == 'matlab':
                        line += '\n'
                    file.write(line)
                    
                    line = '    Kc '
                    if lang in ['c', 'cuda']:
                        line += '-= '
                    elif lang in ['fortran', 'matlab']:
                        line += '= Kc - '
                    line += '{:.2f} * '.format(nu)
                    line += ('({:.8e} - '.format(sp.lo[6]) + 
                             '{:.8e} + '.format(sp.lo[0]) + 
                             '{:.8e} * '.format(sp.lo[0] - 1.0) + 
                             'logT + T * ('
                             '{:.8e} + T * ('.format(sp.lo[1] / 2.0) + 
                             '{:.8e} + T * ('.format(sp.lo[2] / 6.0) + 
                             '{:.8e} + '.format(sp.lo[3] / 12.0) + 
                             '{:.8e} * T))) - '.format(sp.lo[4] / 20.0) + 
                             '{:.8e} / T)'.format(sp.lo[5]) + 
                             utils.line_end[lang]
                             )
                    file.write(line)
                    
                    if lang in ['c', 'cuda']:
                        file.write('  } else {\n')
                    elif lang in ['fortran', 'matlab']:
                        file.write('  else\n')
                    
                    line = '    Kc '
                    if lang in ['c', 'cuda']:
                        line += '-= '
                    elif lang in ['fortran', 'matlab']:
                        line += '= Kc - '
                    line += '{:.2f} * '.format(nu)
                    line += ('({:.8e} - '.format(sp.hi[6]) + 
                             '{:.8e} + '.format(sp.hi[0]) + 
                             '{:.8e} * '.format(sp.hi[0] - 1.0) + 
                             'logT + T * ('
                             '{:.8e} + T * ('.format(sp.hi[1] / 2.0) + 
                             '{:.8e} + T * ('.format(sp.hi[2] / 6.0) + 
                             '{:.8e} + '.format(sp.hi[3] / 12.0) + 
                             '{:.8e} * T))) - '.format(sp.hi[4] / 20.0) + 
                             '{:.8e} / T)'.format(sp.hi[5]) + 
                             utils.line_end[lang]
                             )
                    file.write(line)
                    
                    if lang in ['c', 'cuda']:
                        file.write('  }\n\n')
                    elif lang == 'fortran':
                        file.write('  end if\n\n')
                    elif lang == 'matlab':
                        file.write('  end\n\n')
                
                line = ('  Kc = '
                        '{:.8e}'.format((chem.PA / chem.RU)**sum_nu) + 
                        ' * exp(Kc)' + 
                        utils.line_end[lang]
                        )
                file.write(line)
            
            line = '  rev_rxn_rates'
            if lang in ['c', 'cuda']:
                line += '[{}] = '.format(rev_reacs.index(rxn))
            elif lang in ['fortran', 'matlab']:
                line += '({}) = '.format(rev_reacs.index(rxn) + 1)
            
            # reactants (products from forward reaction)
            for sp in rxn.prod:
                isp = next(i for i in xrange(len(specs)) 
                           if specs[i].name == sp)
                nu = rxn.prod_nu[rxn.prod.index(sp)]
            
                # check if stoichiometric coefficient is real or integer
                if isinstance(nu, float):
                    if lang in ['c', 'cuda']:
                        line += 'pow(C[{}], {}) * '.format(isp, nu)
                    elif lang in ['fortran', 'matlab']:
                        line += 'pow(C({}), {}) * '.format(isp + 1, nu)
                else:
                    # integer, so just use multiplication
                    for i in range(nu):
                        if lang in ['c', 'cuda']:
                            line += 'C[{}] * '.format(isp)
                        elif lang in ['fortran', 'matlab']:
                            line += 'C({}) * '.format(isp + 1)
        
            # rate constant
            if rxn.rev_par:
                # explicit reverse Arrhenius parameters
                line += rxn_rate_const(rxn.rev_par[0], 
                                       rxn.rev_par[1], 
                                       rxn.rev_par[2]
                                       )
            else:
                # use equilibrium constant
                line += 'kf / Kc'
            line += utils.line_end[lang]
            file.write(line)
            
    
    if lang in ['c', 'cuda']:
        file.write('} // end eval_rxn_rates\n\n')
    elif lang == 'fortran':
        file.write('end subroutine eval_rxn_rates\n\n')
    elif lang == 'matlab':
        file.write('end\n\n')
    
    file.close()
    
    return


def write_rxn_pressure_mod(path, lang, specs, reacs):
    """Write subroutine to for reaction pressure dependence modifications.
    
    Parameters
    ----------
    path : str
        Path to build directory for file.
    lang : {'c', 'cuda', 'fortran', 'matlab'}
        Language type.
    specs : list of SpecInfo
        List of species in mechanism.
    reacs : list of ReacInfo
        List of reactions in mechanism.
    
    Returns
    -------
    None
    
    """
    filename = 'rxn_rates_pres_mod' + utils.file_ext[lang]
    file = open(path + filename, 'w')
    
    # headers
    if lang in ['c', 'cuda']:
        file.write('#include <math.h>\n'
                   '#include "header.h"\n'
                   '\n'
                   )
    
    # list of reactions with third-body or pressure-dependence
    pdep_reacs = []
    thd_flag = False
    pdep_flag = False
    troe_flag = False
    sri_flag = False
    for reac in reacs:
        if reac.thd:
            # add reaction index to list
            thd_flag = True
            pdep_reacs.append(reacs.index(reac))
        if reac.pdep:
            # add reaction index to list
            pdep_flag = True
            if not reac.thd: pdep_reacs.append(reacs.index(reac))
            
            if reac.troe and not troe_flag: troe_flag = True
            if reac.sri and not sri_flag: sri_flag = True
    
    line = ''
    if lang == 'cuda': line = '__device__ '
    
    if lang in ['c', 'cuda']:
        line += ('void get_rxn_pres_mod (const Real T, const Real pres, '
                 'const Real * C, Real * pres_mod) {\n'
                 )
    elif lang == 'fortran':
        line += 'subroutine get_rxn_pres_mod ( T, pres, C, pres_mod )\n\n'
        
        # fortran needs type declarations
        line += ('  implicit none\n'
                 '  double precision, intent(in) :: T, pres, '
                 'C({})\n'.format(len(specs)) + 
                 '  double precision, intent(out) :: '
                 'pres_mod({})\n'.format(len(pdep_reacs)) + 
                 '  \n'
                 '  double precision :: logT, m\n')
    elif lang == 'matlab':
        line += ('function pres_mod = get_rxn_pres_mod (T, pres, C)\n\n'
                 '  pres_mod = zeros({},1);\n'.format(len(pdep_reacs))
                 )
    file.write(line)
    
    # declarations for third-body variables
    if thd_flag:
        if lang == 'c':
            file.write('  // third body variable declaration\n'
                       '  Real thd;\n'
                       '\n'
                       )
        elif lang == 'cuda':
            file.write('  // third body variable declaration\n'
                       '  register Real thd;\n'
                       '\n'
                       )
        elif lang == 'fortran':
            file.write('  ! third body variable declaration\n'
                       '  double precision :: thd\n'
                       )
    
    # declarations for pressure-dependence variables
    if pdep_flag:
        if lang == 'c':
            file.write('  // pressure dependence variable declarations\n'
                       '  Real k0;\n'
                       '  Real kinf;\n'
                       '  Real Pr;\n'
                       '\n'
                       )
            if troe_flag:
                # troe variables
                file.write('  // troe variable declarations\n'
                           '  Real logFcent;\n'
                           '  Real A;\n'
                           '  Real B;\n'
                           '\n'
                           )
            if sri_flag:
                # sri variables
                file.write('  // sri variable declarations\n')
                file.write('  Real x;\n'
                           '\n'
                           )
        elif lang == 'cuda':
            file.write('  // pressure dependence variable declarations\n')
#            if not thd_flag: file.write('  register Real thd;\n')
            file.write('  register Real k0;\n'
                       '  register Real kinf;\n'
                       '  register Real Pr;\n'
                       '\n'
                       )
            if troe_flag:
                # troe variables
                file.write('  // troe variable declarations\n'
                           '  register Real logFcent;\n'
                           '  register Real A;\n'
                           '  register Real B;\n'
                           '\n'
                           )
            if sri_flag:
                # sri variables
                file.write('  // sri variable declarations\n')
                file.write('  register Real x;\n'
                           '\n')
        elif lang == 'fortran':
            file.write('  ! pressure dependence variable declarations\n'
                       '  double precision :: k0, kinf, Pr\n'
                       '\n'
                       )
            if troe_flag:
                # troe variables
                file.write('  ! troe variable declarations\n'
                           '  double precision :: logFcent, A, B\n'
                           '\n'
                           )
            if sri_flag:
                # sri variables
                file.write('  ! sri variable declarations\n')
                file.write('  double precision :: X\n'
                           '\n')
    
    if lang == 'c':
        file.write('  Real logT = log(T);\n'
                   '  Real m = pres / ({:.8e} * T);\n'.format(chem.RU)
                   )
    elif lang == 'cuda':
        file.write('  register Real logT = log(T);\n'
                   '  register Real m = pres / ('
                   '{:.8e} * T);\n'.format(chem.RU)
                   )
    elif lang == 'fortran':
        file.write('  logT = log(T)\n'
                   '  m = pres / ({:.8e} * T)\n'.format(chem.RU)
                   )
    elif lang == 'matlab':
        file.write('  logT = log(T);\n'
                   '  m = pres / ({:.8e} * T);\n'.format(chem.RU)
                   )
    
    file.write('\n')
    
    # loop through third-body and pressure-dependent reactions
    for rind in pdep_reacs:
        reac = reacs[rind]              # index in reaction list
        pind = pdep_reacs.index(rind)   # index in list of third/pressure-dep reactions
        
        # print reaction index
        if lang in ['c', 'cuda']:
            line = '  // reaction ' + str(rind)
        elif lang == 'fortran':
            line = '  ! reaction ' + str(rind + 1)
        elif lang == 'matlab':
            line = '  % reaction ' + str(rind + 1)
        line += utils.line_end[lang]
        file.write(line)
        
        # third-body reaction
        if reac.thd:
            
            if reac.pdep and not reac.pdep_sp:
                line = '  thd = m'
            else:
                if lang in ['c', 'cuda']:
                    line = '  pres_mod[{}] = m'.format(pind)
                elif lang in ['fortran', 'matlab']:
                    line = '  pres_mod({}) = m'.format(pind + 1)
            
            for sp in reac.thd_body:
                isp = specs.index(next((s for s in specs 
                                  if s.name == sp[0]), None)
                                  )
                if sp[1] > 1.0:
                    line += ' + {}'.format(sp[1] - 1.0)
                elif sp[1] < 1.0:
                    line += ' - {}'.format(1.0 - sp[1])
                if lang in ['c', 'cuda']:
                    line += ' * C[{}]'.format(isp)
                elif lang in ['fortran', 'matlab']:
                    line += ' * C({})'.format(isp + 1)
            
            line += utils.line_end[lang]
            file.write(line)
        
        # pressure dependence
        if reac.pdep:
            
            # low-pressure limit rate
            line = '  k0 = '
            if reac.low:
                line += rxn_rate_const(reac.low[0], 
                                       reac.low[1], 
                                       reac.low[2]
                                       )
            else:
                line += rxn_rate_const(reac.A, reac.b, reac.E)
            
            line += utils.line_end[lang]
            file.write(line)
            
            # high-pressure limit rate
            line = '  kinf = '
            if reac.high:
                line += rxn_rate_const(reac.high[0], 
                                       reac.high[1], 
                                       reac.high[2]
                                       )
            else:
                line += rxn_rate_const(reac.A, reac.b, reac.E)
            
            line += utils.line_end[lang]
            file.write(line)
            
            # reduced pressure
            if reac.thd:
                line = '  Pr = k0 * thd / kinf'
            else:
                isp = next(i for i in xrange(len(specs))
                           if specs[i].name == reac.pdep_sp
                           )
                if lang in ['c', 'cuda']:
                    line = '  Pr = k0 * C[{}] / kinf'.format(isp)
                elif lang in ['fortran', 'matlab']:
                    line = '  Pr = k0 * C({}) / kinf'.format(isp + 1)
            line += utils.line_end[lang]
            file.write(line)
            
            simple = False
            if reac.troe:
                # Troe form
                line = ('  logFcent = log10( fmax('
                        '{:.8e} * '.format(1.0 - reac.troe_par[0])
                        )
                if reac.troe_par[1] > 0.0:
                    line += 'exp(-T / {:.8e})'.format(reac.troe_par[1])
                else:
                    line += 'exp(T / {:.8e})'.format(abs(reac.troe_par[1]))
                
                line += ' + {:.8e} * '.format(reac.troe_par[0])
                if reac.troe_par[2] > 0.0:
                    line += 'exp(-T / {:.8e})'.format(reac.troe_par[2])
                else:
                    line += 'exp(T / {:.8e})'.format(abs(reac.troe_par[2]))
                
                if len(reac.troe_par) == 4:
                    line += ' + '
                    if reac.troe_par[3] > 0.0:
                        val = reac.troe_par[3]
                        line += 'exp(-{:.8e} / T)'.format(val)
                    else:
                        val = abs(reac.troe_par[3])
                        line += 'exp({:.8e} / T)'.format(val)
                line += ', 1.0e-300))' + utils.line_end[lang]
                file.write(line)
                
                line = ('  A = log10(fmax(Pr, 1.0e-300)) - '
                        '0.67 * logFcent - 0.4' + 
                        utils.line_end[lang]
                        )
                file.write(line)
                
                line = ('  B = 0.806 - 1.1762 * logFcent - '
                        '0.14 * log10(fmax(Pr, 1.0e-300))' +
                        utils.line_end[lang]
                        )
                file.write(line)
                
                line = '  pres_mod'
                if lang in ['c', 'cuda']:
                    line += ('[{}]'.format(pind) + 
                             ' = ' + utils.exp_10_fun[lang]
                             )
                elif lang in ['fortran', 'matlab']:
                    # fortran & matlab don't have exp10
                    line += ('({})'.format(pind + 1) + 
                             ' = ' + utils.exp_10_fun[lang]
                             )
                line += 'logFcent / (1.0 + A * A / (B * B))) '
                
            elif reac.sri:
                # SRI form
                
                line = ('  X = 1.0 / (1.0 + log10(fmax(Pr, 1.0e-300)) * '
                        'log10(fmax(Pr, 1.0e-300)))' + 
                        utils.line_end[lang]
                        )
                file.write(line)
                
                line = '  pres_mod'
                if lang in ['c', 'cuda']:
                    line += '[{}]'.format(pind)
                elif lang in ['fortran', 'matlab']:
                    line += '({})'.format(pind + 1)
                line += ' = pow({:4} * '.format(reac.sri[0])
                # Need to check for negative parameters, and 
                # skip "-" sign if so.
                if reac.sri[1] > 0.0:
                    line += 'exp(-{:.4} / T)'.format(reac.sri[1])
                else:
                    line += 'exp({:.4} / T)'.format(abs(reac.sri[1]))
                
                if reac.sri[2] > 0.0:
                    line += ' + exp(-T / {:.4}), X) '.format(reac.sri[2])
                else:
                    line += ' + exp(T / {:.4}), X) '.format(abs(reac.sri[2]))
                    
                if len(reac.sri) == 5:
                    line += ('* {:.8e} * '.format(reac.sri[3]) + 
                             'pow(T, {:.4}) '.format(reac.sri[4])
                             )
            else:
                #simple falloff fn (i.e. F = 1)
                simple = True
                line = '  pres_mod'
                if lang in ['c', 'cuda']:
                    line += '[{}] ='.format(pind)
                elif lang in ['fortran', 'matlab']:
                    line += '({}) ='.format(pind + 1)
                 # regardless of F formulation
                if reac.low:
                    # unimolecular/recombination fall-off reaction
                    line += ' Pr / (1.0 + Pr)'
                elif reac.high:
                    # chemically-activated bimolecular reaction
                    line += '1.0 / (1.0 + Pr)'
            
            if not simple:
                # regardless of F formulation
                if reac.low:
                    # unimolecular/recombination fall-off reaction
                    line += '* Pr / (1.0 + Pr)'
                elif reac.high:
                    # chemically-activated bimolecular reaction
                    line += '/ (1.0 + Pr)'
            
            line += utils.line_end[lang]
            file.write(line)
        
        # space in between each reaction
        file.write('\n')
    
    if lang in ['c', 'cuda']:
        file.write('} // end get_rxn_pres_mod\n\n')
    elif lang == 'fortran':
        file.write('end subroutine get_rxn_pres_mod\n\n')
    elif lang == 'matlab':
        file.write('end\n\n')
    
    file.close()
    
    return


def write_spec_rates(path, lang, specs, reacs):
    """Write subroutine to evaluate species rates of production.
    
    Parameters
    ----------
    path : str
        Path to build directory for file.
    lang : {'c', 'cuda', 'fortran', 'matlab'}
        Programming language.
    specs : list of SpecInfo
        List of species in mechanism.
    reacs : list of ReacInfo
        List of reactions in mechanism.
    
    Returns
    -------
    None
    
    """
    
    filename = 'spec_rates' + utils.file_ext[lang]
    file = open(path + filename, 'w')
    
    if lang in ['c', 'cuda']:
        file.write('#include "header.h"\n'
                   '\n'
                   )
    
    num_s = len(specs)
    num_r = len(reacs)
    rev_reacs = [rxn for rxn in reacs if rxn.rev]
    num_rev = len(rev_reacs)
    
    # pressure dependent reactions
    pdep_reacs = []
    for reac in reacs:
        if reac.thd or reac.pdep:
            # add reaction index to list
            pdep_reacs.append(reacs.index(reac))
    num_pdep = len(pdep_reacs)
    
    line = ''
    if lang == 'cuda': line = '__device__ '
    
    if lang in ['c', 'cuda']:
        if rev_reacs:
            line += ('void eval_spec_rates (const Real * fwd_rates, '
                     'const Real * rev_rates, const Real * pres_mod, '
                     'Real * sp_rates) {\n'
                     )
        else:
            line += ('void eval_spec_rates (const Real * fwd_rates, '
                     'const Real * pres_mod, Real * sp_rates) {\n'
                     )
    elif lang == 'fortran':
        if rev_reacs:
            line += ('subroutine eval_spec_rates (fwd_rates, rev_rates, '
                     'pres_mod, sp_rates)\n\n'
                     )
        else:
            line += ('subroutine eval_spec_rates (fwd_rates, pres_mod, '
                     'sp_rates)\n\n'
                     )
        
        # fortran needs type declarations
        line += '  implicit none\n'
        if rev_reacs:
            line += ('  double precision, intent(in) :: '
                     'fwd_rates({0}), rev_rates({0}), '.format(num_r) + 
                     'pres_mod({})\n'.format(num_pdep)
                     )
        else:
            line += ('  double precision, intent(in) :: '
                     'fwd_rates({}), '.format(num_r) + 
                     'pres_mod({})\n'.format(num_pdep)
                     )
        line += ('  double precision, intent(out) :: '
                 'sp_rates({})\n'.format(num_s) + 
                 '\n'
                 )
    elif lang == 'matlab':
        if rev_reacs:
            line += ('function sp_rates = eval_spec_rates ( fwd_rates, '
                     'rev_rates, pres_mod )\n\n'
                     )
        else:
            line += ('function sp_rates = eval_spec_rates ( fwd_rates, '
                     'pres_mod )\n\n'
                     )
        line += '  sp_rates = zeros({},1);\n'.format(len(specs))
    file.write(line)
    
    # loop through species
    for sp in specs:
        line = '  sp_rates'
        if lang in ['c', 'cuda']:
            line += '[{}] = '.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line += '({}) = '.format(specs.index(sp) + 1)
        
        # continuation line
        cline = ' ' * ( len(line) - 3)
        
        isfirst = True
        
        inreac = False
        
        # loop through reactions
        for rxn in reacs:
            
            rind = reacs.index(rxn)
            
            pdep = False
            if rxn.thd or rxn.pdep: pdep = True
            
            # move to new line if current line is too long
            if len(line) > 85:
                line += '\n'
                # record position
                lastPos = file.tell()
                file.write(line)
                line = cline
            
            # first check to see if in both products and reactants
            if sp.name in rxn.prod and sp.name in rxn.reac:
                inreac = True
                pisp = rxn.prod.index(sp.name)
                risp = rxn.reac.index(sp.name)
                nu = rxn.prod_nu[pisp] - rxn.reac_nu[risp]
                
                if nu > 0.0:
                    if not isfirst: line += ' + '
                    if nu > 1:
                        if isinstance(nu, int):
                            line += '{} * '.format(float(nu))
                        else:
                            line += '{:3} * '.format(nu)
                    elif nu < 1.0:
                        line += '{} * '.format(nu)
                    
                    if rxn.rev:
                        line += '(fwd_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rind)
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rind + 1)
                        line += ' - rev_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rev_reacs.index(rxn))
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rev_reacs.index(rxn) + 1)
                        line += ')'
                    else:
                        line += 'fwd_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rind)
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rind + 1)
                elif nu < 0.0:
                    if isfirst:
                        line += '-'
                    else:
                        line += ' - '
                    
                    if nu < -1:
                        if isinstance(nu, int):
                            line += '{} * '.format(float(abs(nu)))
                        else:
                            line += '{:3} * '.format(abs(nu))
                    elif nu > -1:
                        line += '{} * '.format(abs(nu))
                    
                    if rxn.rev:
                        line += '(fwd_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rind)
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rind + 1)
                        line += ' - rev_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rev_reacs.index(rxn))
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rev_reacs.index(rxn) + 1)
                        line += ')'
                    else:
                        line += 'fwd_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rind)
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rind + 1)
                else:
                    inreac = False
                    continue
                
                if isfirst: isfirst = False
                
            # check products
            elif sp.name in rxn.prod:
                inreac = True
                isp = rxn.prod.index(sp.name)
                nu = rxn.prod_nu[isp]
                
                if not isfirst: line += ' + '
                
                if nu > 1:
                    if isinstance(nu, int):
                        line += '{} * '.format(float(nu))
                    else:
                        line += '{:3} * '.format(nu)
                elif nu < 1.0:
                    line += '{} * '.format(nu)
                
                if rxn.rev:
                        line += '(fwd_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rind)
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rind + 1)
                        line += ' - rev_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rev_reacs.index(rxn))
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rev_reacs.index(rxn) + 1)
                        line += ')'
                else:
                    line += 'fwd_rates'
                    if lang in ['c', 'cuda']:
                        line += '[{}]'.format(rind)
                    elif lang in ['fortran', 'matlab']:
                        line += '({})'.format(rind + 1)
                
                if isfirst: isfirst = False
                
            # check reactants
            elif sp.name in rxn.reac:
                inreac = True
                isp = rxn.reac.index(sp.name)
                nu = rxn.reac_nu[isp]
                
                if isfirst:
                    line += '-'
                else:
                    line += ' - '
                
                if nu > 1:
                    if isinstance(nu, int):
                        line += '{} * '.format(float(nu))
                    else:
                        line += '{:3} * '.format(nu)
                elif nu < 1.0:
                    line += '{} * '.format(nu)
                
                if rxn.rev:
                        line += '(fwd_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rind)
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rind + 1)
                        line += ' - rev_rates'
                        if lang in ['c', 'cuda']:
                            line += '[{}]'.format(rev_reacs.index(rxn))
                        elif lang in ['fortran', 'matlab']:
                            line += '({})'.format(rev_reacs.index(rxn) + 1)
                        line += ')'
                else:
                    line += 'fwd_rates'
                    if lang in ['c', 'cuda']:
                        line += '[{}]'.format(rind)
                    elif lang in ['fortran', 'matlab']:
                        line += '({})'.format(rind + 1)
                
                if isfirst: isfirst = False
            else:
                continue
            
            # pressure dependence modification
            if pdep:
                pind = pdep_reacs.index(rind)
                if lang in ['c', 'cuda']:
                    line += ' * pres_mod[{}]'.format(pind)
                elif lang in ['fortran', 'matlab']:
                    line += ' * pres_mod({})'.format(pind + 1)
        
        # species not participate in any reactions
        if not inreac: line += '0.0'
        
        # done with this species
        line += utils.line_end[lang] + '\n'
        file.write(line)
    
    if lang in ['c', 'cuda']:
        file.write('} // end eval_spec_rates\n\n')
    elif lang == 'fortran':
        file.write('end subroutine eval_spec_rates\n\n')
    elif lang == 'matlab':
        file.write('end\n\n')
    
    file.close()
    
    return

def write_chem_utils(path, lang, specs):
    """Write subroutine to evaluate species thermodynamic properties.
    
    Notes
    -----
    Thermodynamic properties include:  enthalpy, energy, specific heat
    (constant pressure and volume).
    
    Parameters
    ----------
    path : str
        Path to build directory for file.
    lang : {'c', 'cuda', 'fortran', 'matlab'}
        Programming language.
    specs : list of SpecInfo
        List of species in the mechanism.
    
    Returns
    -------
    None
    
    """
    
    num_s = len(specs)
    
    # first write header file
    if lang == 'c':
        file = open(path + 'chem_utils.h', 'w')
        file.write('#ifndef CHEM_UTILS_HEAD\n'
                   '#define CHEM_UTILS_HEAD\n'
                   '\n'
                   '#include "header.h"\n'
                   '\n'
                   'void eval_h (const Real, Real*);\n'
                   'void eval_u (const Real, Real*);\n'
                   'void eval_cv (const Real, Real*);\n'
                   'void eval_cp (const Real, Real*);\n'
                   '\n'
                   '#endif\n'
                   )
        file.close()
    elif lang == 'cuda':
        file = open(path + 'chem_utils.cuh', 'w')
        file.write('#ifndef CHEM_UTILS_HEAD\n'
                   '#define CHEM_UTILS_HEAD\n'
                   '\n'
                   '#include "header.h"\n'
                   '\n'
                   '__device__ void eval_h (const Real, Real*);\n'
                   '__device__ void eval_u (const Real, Real*);\n'
                   '__device__ void eval_cv (const Real, Real*);\n'
                   '__device__ void eval_cp (const Real, Real*);\n'
                   '\n'
                   '#endif\n'
                   )
        file.close()
    
    filename = 'chem_utils' + utils.file_ext[lang]
    file = open(path + filename, 'w')
    
    if lang in ['c', 'cuda']:
        file.write('#include "header.h"\n\n')
    
    pre = ''
    if lang == 'cuda': pre = '__device__ '
    
    ######################
    # enthalpy subroutine
    ######################
    line = pre
    if lang in ['c', 'cuda']:
        line += 'void eval_h (const Real T, Real * h) {\n\n'
    elif lang == 'fortran':
        line += ('subroutine eval_h (T, h)\n\n'
                 # fortran needs type declarations
                 '  implicit none\n'
                 '  double precision, intent(in) :: T\n'
                 '  double precision, intent(out) :: h({})\n'.format(num_s) + 
                 '\n'
                 )
    elif lang == 'matlab':
        line += 'function h = eval_h (T)\n\n'
    file.write(line)
    
    # loop through species
    for sp in specs:
        line = '  if (T <= {:})'.format(sp.Trange[1])
        if lang in ['c', 'cuda']:
            line += ' {\n'
        elif lang == 'fortran':
            line += ' then\n'
        elif lang == 'matlab':
            line += '\n'
        file.write(line)
        
        if lang in ['c', 'cuda']:
            line = '    h[{}]'.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line = '    h({})'.format(specs.index(sp) + 1)
        line += (' = {:.8e} * '.format(chem.RU / sp.mw) + 
                 '({:.8e} + T * ('.format(sp.lo[5]) + 
                 '{:.8e} + T * ('.format(sp.lo[0]) + 
                 '{:.8e} + T * ('.format(sp.lo[1] / 2.0) + 
                 '{:.8e} + T * ('.format(sp.lo[2] / 3.0) + 
                 '{:.8e} + '.format(sp.lo[3] / 4.0) + 
                 '{:.8e} * T)))))'.format(sp.lo[4] / 5.0) + 
                 utils.line_end[lang]
                 )
        file.write(line)
        
        if lang in ['c', 'cuda']:
            file.write('  } else {\n')
        elif lang in ['fortran', 'matlab']:
            file.write('  else\n')
        
        if lang in ['c', 'cuda']:
            line = '    h[{}]'.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line = '    h({})'.format(specs.index(sp) + 1)
        line += (' = {:.8e} * '.format(chem.RU / sp.mw) + 
                 '({:.8e} + T * ('.format(sp.hi[5]) + 
                 '{:.8e} + T * ('.format(sp.hi[0]) + 
                 '{:.8e} + T * ('.format(sp.hi[1] / 2.0) + 
                 '{:.8e} + T * ('.format(sp.hi[2] / 3.0) + 
                 '{:.8e} + '.format(sp.hi[3] / 4.0) + 
                 '{:.8e} * T)))))'.format(sp.hi[4] / 5.0) + 
                 utils.line_end[lang]
                 )
        file.write(line)
        
        if lang in ['c', 'cuda']:
            file.write('  }\n\n')
        elif lang == 'fortran':
            file.write('  end if\n\n')
        elif lang == 'matlab':
            file.write('  end\n\n')
    
    if lang in ['c', 'cuda']:
        file.write('} // end eval_h\n\n')
    elif lang == 'fortran':
        file.write('end subroutine eval_h\n\n')
    elif lang == 'matlab':
        file.write('end\n\n')
    
    #################################
    # internal energy subroutine
    #################################
    line = pre
    if lang in ['c', 'cuda']:
        line += 'void eval_u (const Real T, Real * u) {\n\n'
    elif lang == 'fortran':
        line += ('subroutine eval_u (T, u)\n\n'
                 # fortran needs type declarations
                 '  implicit none\n'
                 '  double precision, intent(in) :: T\n'
                 '  double precision, intent(out) :: u({})\n'.format(num_s) + 
                 '\n')
    elif lang == 'matlab':
        line += 'function u = eval_u (T)\n\n'
    file.write(line)
    
    # loop through species
    for sp in specs:
        line = '  if (T <= {:})'.format(sp.Trange[1])
        if lang in ['c', 'cuda']:
            line += ' {\n'
        elif lang == 'fortran':
            line += ' then\n'
        elif lang == 'matlab':
            line += '\n'
        file.write(line)
        
        if lang in ['c', 'cuda']:
            line = '    u[{}]'.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line = '    u({})'.format(specs.index(sp) + 1)
        line += (' = {:.8e} * '.format(chem.RU / sp.mw) + 
                 '({:.8e} + T * ('.format(sp.lo[5]) + 
                 '{:.8e} - 1.0 + T * ('.format(sp.lo[0]) + 
                 '{:.8e} + T * ('.format(sp.lo[1] / 2.0) + 
                 '{:.8e} + T * ('.format(sp.lo[2] / 3.0) + 
                 '{:.8e} + '.format(sp.lo[3] / 4.0) + 
                 '{:.8e} * T)))))'.format(sp.lo[4] / 5.0) + 
                 utils.line_end[lang]
                 )
        file.write(line)
        
        if lang in ['c', 'cuda']:
            file.write('  } else {\n')
        elif lang in ['fortran', 'matlab']:
            file.write('  else\n')
        
        if lang in ['c', 'cuda']:
            line = '    u[{}]'.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line = '    u({})'.format(specs.index(sp) + 1)
        line += (' = {:.8e} * '.format(chem.RU / sp.mw) + 
                 '({:.8e} + T * ('.format(sp.hi[5]) + 
                 '{:.8e} - 1.0 + T * ('.format(sp.hi[0]) + 
                 '{:.8e} + T * ('.format(sp.hi[1] / 2.0) + 
                 '{:.8e} + T * ('.format(sp.hi[2] / 3.0) + 
                 '{:.8e} + '.format(sp.hi[3] / 4.0) + 
                 '{:.8e} * T)))))'.format(sp.hi[4] / 5.0) + 
                 utils.line_end[lang]
                 )
        file.write(line)
        
        if lang in ['c', 'cuda']:
            file.write('  }\n\n')
        elif lang == 'fortran':
            file.write('  end if\n\n')
        elif lang == 'matlab':
            file.write('  end\n\n')
    
    if lang in ['c', 'cuda']:
        file.write('} // end eval_u\n\n')
    elif lang == 'fortran':
        file.write('end subroutine eval_u\n\n')
    elif lang == 'matlab':
        file.write('end\n\n')
    
    ##################################
    # cv subroutine
    ##################################
    if lang in ['c', 'cuda']:
        line = pre + 'void eval_cv (const Real T, Real * cv) {\n\n'
    elif lang == 'fortran':
        line = ('subroutine eval_cv (T, cv)\n\n'
                # fortran needs type declarations
                '  implicit none\n'
                '  double precision, intent(in) :: T\n'
                '  double precision, intent(out) :: cv({})\n'.format(num_s) + 
                '\n'
                )
    elif lang == 'matlab':
        line = 'function cv = eval_cv (T)\n\n'
    file.write(line)
    
    # loop through species
    for sp in specs:
        line = '  if (T <= {:})'.format(sp.Trange[1])
        if lang in ['c', 'cuda']:
            line += ' {\n'
        elif lang == 'fortran':
            line += ' then\n'
        elif lang == 'matlab':
            line += '\n'
        file.write(line)
        
        if lang in ['c', 'cuda']:
            line = '    cv[{}]'.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line = '    cv({})'.format(specs.index(sp) + 1)
        line += (' = {:.8e} * '.format(chem.RU / sp.mw) + 
                 '({:.8e} - 1.0 + T * ('.format(sp.lo[0]) + 
                 '{:.8e} + T * ('.format(sp.lo[1]) + 
                 '{:.8e} + T * ('.format(sp.lo[2]) + 
                 '{:.8e} + '.format(sp.lo[3]) + 
                 '{:.8e} * T))))'.format(sp.lo[4]) + 
                 utils.line_end[lang]
                 )
        file.write(line)
        
        if lang in ['c', 'cuda']:
            file.write('  } else {\n')
        elif lang in ['fortran', 'matlab']:
            file.write('  else\n')
        
        if lang in ['c', 'cuda']:
            line = '    cv[{}]'.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line = '    cv({})'.format(specs.index(sp) + 1)
        line += (' = {:.8e} * '.format(chem.RU / sp.mw) + 
                 '({:.8e} - 1.0 + T * ('.format(sp.hi[0]) + 
                 '{:.8e} + T * ('.format(sp.hi[1]) + 
                 '{:.8e} + T * ('.format(sp.hi[2]) + 
                 '{:.8e} + '.format(sp.hi[3]) + 
                 '{:.8e} * T))))'.format(sp.hi[4]) +
                 utils.line_end[lang]
                 )
        file.write(line)
        
        if lang in ['c', 'cuda']:
            file.write('  }\n\n')
        elif lang == 'fortran':
            file.write('  end if\n\n')
        elif lang == 'matlab':
            file.write('  end\n\n')
    
    if lang in ['c', 'cuda']:
        file.write('} // end eval_cv\n\n')
    elif lang == 'fortran':
        file.write('end subroutine eval_cv\n\n')
    elif lang == 'matlab':
        file.write('end\n\n')
    
    ###############################
    # cp subroutine 
    ###############################
    if lang in ['c', 'cuda']:
        line = pre + 'void eval_cp (const Real T, Real * cp) {\n\n'
    elif lang == 'fortran':
        line = ('subroutine eval_cp (T, cp)\n\n'
                # fortran needs type declarations
                '  implicit none\n'
                '  double precision, intent(in) :: T\n'
                '  double precision, intent(out) :: cp({})\n'.format(num_s) + 
                '\n'
                )
    elif lang == 'matlab':
        line = 'function cp = eval_cp (T)\n\n'
    file.write(line)
    
    # loop through species
    for sp in specs:
        line = '  if (T <= {:})'.format(sp.Trange[1])
        if lang in ['c', 'cuda']:
            line += ' {\n'
        elif lang == 'fortran':
            line += ' then\n'
        elif lang == 'matlab':
            line += '\n'
        file.write(line)
        
        if lang in ['c', 'cuda']:
            line = '    cp[{}]'.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line = '    cp({})'.format(specs.index(sp) + 1)
        line += (' = {:.8e} * '.format(chem.RU / sp.mw) + 
                 '({:.8e} + T * ('.format(sp.lo[0]) + 
                 '{:.8e} + T * ('.format(sp.lo[1]) + 
                 '{:.8e} + T * ('.format(sp.lo[2]) + 
                 '{:.8e} + '.format(sp.lo[3]) + 
                 '{:.8e} * T))))'.format(sp.lo[4]) + 
                 utils.line_end[lang]
                 )
        file.write(line)
        
        if lang in ['c', 'cuda']:
            file.write('  } else {\n')
        elif lang in ['fortran', 'matlab']:
            file.write('  else\n')
        
        if lang in ['c', 'cuda']:
            line = '    cp[{}]'.format(specs.index(sp))
        elif lang in ['fortran', 'matlab']:
            line = '    cp({})'.format(specs.index(sp) + 1)
        line += (' = {:.8e} * '.format(chem.RU / sp.mw) + 
                 '({:.8e} + T * ('.format(sp.hi[0]) + 
                 '{:.8e} + T * ('.format(sp.hi[1]) + 
                 '{:.8e} + T * ('.format(sp.hi[2]) + 
                 '{:.8e} + '.format(sp.hi[3]) + 
                 '{:.8e} * T))))'.format(sp.hi[4]) + 
                 utils.line_end[lang]
                 )
        file.write(line)
        
        if lang in ['c', 'cuda']:
            file.write('  }\n\n')
        elif lang == 'fortran':
            file.write('  end if\n\n')
        elif lang == 'matlab':
            file.write('  end\n\n')
    
    if lang in ['c', 'cuda']:
        file.write('} // end eval_cp\n\n')
    elif lang == 'fortran':
        file.write('end subroutine eval_cp\n\n')
    elif lang == 'matlab':
        file.write('end\n\n')
    
    file.close()
    
    return

def write_derivs(path, lang, specs, reacs):
    """Writes derivative function file and header.
    
    Parameters
    ----------
    path : str
        Path to build directory for file.
    lang : {'c', 'cuda', 'fortran', 'matlab'}
        Programming language.
    specs : list of SpecInfo
        List of species in the mechanism.
    reacs : list of ReacInfo
        List of reactions in the mechanism.
    
    Returns
    -------
    None
    
    """
    
    # first write header file
    if lang == 'c':
        file = open(path + 'dydt.h', 'w')
        file.write('#ifndef DYDT_HEAD\n'
                   '#define DYDT_HEAD\n'
                   '\n'
                   '#include "header.h"\n'
                   '\n'
                   'void dydt (const Real, const Real, '
                   'const Real*, Real*);\n'
                   '\n'
                   '#endif\n'
                   )
        file.close()
    elif lang == 'cuda':
        file = open(path + 'dydt.cuh', 'w')
        file.write('#ifndef DYDT_HEAD\n'
                   '#define DYDT_HEAD\n'
                   '\n'
                   '#include "header.h"\n'
                   '\n'
                   '__device__ void dydt (const Real, const Real, '
                   'const Real*, Real*);\n'
                   '\n'
                   '#endif\n'
                   )
        file.close()
    
    filename = 'dydt' + utils.file_ext[lang]
    file = open(path + filename, 'w')

    pre = ''
    if lang == 'cuda': pre = '__device__ '
    
    file.write('#include "header.h"\n')
    if lang == 'c':
        file.write('#include "chem_utils.h"\n'
                   '#include "rates.h"\n'
                   )
    elif lang == 'cuda':
        file.write('#include "chem_utils.cuh"\n'
                   '#include "rates.cuh"\n'
                   )
    file.write('\n')
    
    # constant pressure
    file.write('#if defined(CONP)\n\n')
    
    line = (pre + 'void dydt (const Real t, const Real pres, '
            'const Real * y, Real * dy) {\n\n'
            )
    file.write(line)
    
    # calculation of density
    file.write('  // mass-averaged density\n'
               '  Real rho;\n'
               )
    line = '  rho = '
    isfirst = True
    for sp in specs:
        if len(line) > 70:
            line += '\n'
            file.write(line)
            line = '     '
        
        if not isfirst: line += ' + '
        if lang in ['c', 'cuda']:
            line += '(y[{}] / {})'.format(specs.index(sp) + 1, sp.mw)
        elif lang in ['fortran', 'matlab']:
            line += '(y[{}] / {})'.format(specs.index(sp) + 2, sp.mw)
        
        isfirst = False
    
    line += ';\n'
    file.write(line)
    line = '  rho = pres / ({:.8e} * y[0] * rho);\n\n'.format(chem.RU)
    file.write(line)
    
    # calculation of species molar concentrations
    file.write('  // species molar concentrations\n'
               '  Real conc[{}];\n'.format(len(specs))
               )
    # loop through species
    for sp in specs:
        isp = specs.index(sp)
        line = '  conc'
        if lang in ['c', 'cuda']:
            line += '[{}] = rho * y[{}] / '.format(isp, isp + 1)
        elif lang in ['fortran', 'matlab']:
            line += '({}) = rho * y({}) / '.format(isp + 1, isp + 2)
        line += '{}'.format(sp.mw) + utils.line_end[lang]
        file.write(line)
    
    file.write('\n')
    
    # evaluate reaction rates
    rev_reacs = [rxn for rxn in reacs if rxn.rev]
    if rev_reacs:
        file.write('  // local arrays holding reaction rates\n'
                   '  Real fwd_rates[{}];\n'.format(len(reacs)) + 
                   '  Real rev_rates[{}];\n'.format(len(rev_reacs)) + 
                   '  eval_rxn_rates (y[0], conc, fwd_rates, rev_rates);\n'
                   '\n'
                   )
    else:
        file.write('  // local array holding reaction rates\n'
                   '  Real rates[{}];\n'.format(len(reacs)) + 
                   '  eval_rxn_rates (y[0], conc, rates);\n'
                   '\n'
                   )
    
    # reaction pressure dependence
    pdep_reacs = []
    for reac in reacs:
        if reac.thd or reac.pdep:
            # add reaction index to list
            pdep_reacs.append(reacs.index(reac))
    num_pdep = len(pdep_reacs)
    if pdep_reacs:
        if lang in ['c', 'cuda']:
            file.write('  // get pressure modifications to reaction rates\n'
                       '  Real pres_mod[{}];\n'.format(num_pdep) + 
                       '  get_rxn_pres_mod (y[0], pres, conc, pres_mod);\n'
                       )
        elif lang == 'fortran':
            file.write('  ! get and evaluate pressure modifications to '
                       'reaction rates\n'
                       '  get_rxn_pres_mod (y[0], pres, conc, pres_mod)\n'
                       )
        elif lang == 'matlab':
            file.write('  % get and evaluate pressure modifications to '
                       'reaction rates\n'
                       '  pres_mod = get_rxn_pres_mod (y[0], pres, conc, '
                       'pres_mod);\n'
                       )
        file.write('\n')
    
    # species rate of change of molar concentration
    file.write('  // evaluate rate of change of species molar '
               'concentration\n'
               )
    if rev_reacs and pdep_reacs:
        file.write('  eval_spec_rates (fwd_rates, rev_rates, pres_mod, '
                   '&dy[1]);\n'
                   '\n'
                   )
    elif rev_reacs:
        file.write('  eval_spec_rates (fwd_rates, rev_rates, &dy[1]);\n\n')
    else:
        file.write('  eval_spec_rates (rates, &dy[1] );\n\n')
    
    # evaluate specific heat
    file.write('  // local array holding constant pressure specific heat\n'
               '  Real cp[{}];\n'.format(len(specs)) + 
               '  eval_cp (y[0], cp);\n'
               '\n'
               )
    
    file.write('  // constant pressure mass-average specific heat\n')
    line = '  Real cp_avg = '
    isfirst = True
    for sp in specs:
        if len(line) > 70:
            line += '\n'
            file.write(line)
            line = '             '
        
        if not isfirst: line += ' + '
        
        isp = specs.index(sp)
        line += '(cp[{}] * y[{}])'.format(isp, isp + 1)
        
        isfirst = False
    
    line += ';\n\n'
    file.write(line)
    
    # evaluate enthalpy
    file.write('  // local array for species enthalpies\n'
               '  Real h[{}];\n'.format(len(specs)) + 
               '  eval_h (y[0], h);\n'
               '\n'
               )
    
    # energy equation
    file.write('  // rate of change of temperature\n')
    line = '  dy[0] = (-1.0 / (rho * cp_avg)) * ( '
    isfirst = True
    for sp in specs:
        if len(line) > 70:
            line += '\n'
            file.write(line)
            line = '       '
        
        if not isfirst: line += ' + '
        
        isp = specs.index(sp)
        line += '(dy[{}] * h[{}] * {})'.format(isp + 1, isp, sp.mw)
        
        isfirst = False
    
    line += ' );\n\n'
    file.write(line)
    
    # rate of change of species mass fractions
    file.write('  // calculate rate of change of species mass fractions\n')
    for sp in specs:
        line = '  dy[{}] *= ({} / rho);\n'.format(specs.index(sp) + 1, sp.mw)
        file.write(line)
    
    file.write('\n')
    file.write('} // end dydt\n\n')
    
    # constant volume
    file.write('#elif defined(CONV)\n\n')
    
    line = (pre + 'void dydt (const Real t, const Real rho, '
            'const Real * y, Real * dy) {\n'
            '\n'
            )
    file.write(line)
    
    # just use y[0] for temperature
    #file.write('  Real T = y[0];\n\n')
    
    # calculation of pressure
    file.write('  // pressure\n'
               '  Real pres;\n'
               )
    line = '  pres = '
    isfirst = True
    for sp in specs:
        if len(line) > 70:
            line += '\n'
            file.write(line)
            line = '      '
        
        if not isfirst: line += ' + '
        line += '(y[{}] / {})'.format(specs.index(sp) + 1, sp.mw)
        
        isfirst = False
    
    line += ';\n'
    file.write(line)
    line = '  pres = rho * {:.8e} * y[0] * pres;\n\n'.format(chem.RU)
    file.write(line)
    
    # calculation of species molar concentrations
    file.write('  // species molar concentrations\n'
               '  Real conc[{}];\n'.format(len(specs))
               )
    # loop through species
    for sp in specs:
        isp = specs.index(sp)
        line = '  conc[{}] = rho * y[{}] / {};\n'.format(isp, isp + 1, sp.mw)
        file.write(line)
    
    file.write('\n')
    
    # evaluate reaction rates
    file.write('  // local array holding reaction rates\n'
               '  Real rates[{}];\n'.format(len(reacs)) + 
               '  eval_rxn_rates (y[0], pres, conc, rates);\n'
               '\n'
               )
    
    # species rate of change of molar concentration
    file.write('  // evaluate rate of change of species molar '
               'concentration\n'
               '  eval_spec_rates (rates, &dy[1]);\n'
               '\n'
               )
    
    # evaluate specific heat
    file.write('  // local array holding constant volume specific heat\n'
               '  Real cv[{}];\n'.format(len(specs)) + 
               '  eval_cv (y[0], cv);\n'
               '\n'
               )
    
    file.write('  // constant volume mass-average specific heat\n')
    line = '  Real cv_avg = '
    isfirst = True
    for sp in specs:
        if len(line) > 70:
            line += '\n'
            file.write(line)
            line = '             '
        
        if not isfirst: line += ' + '
        
        isp = specs.index(sp)
        line += '(cv[{}] * y[{}])'.format(isp, isp + 1)
        
        isfirst = False
    
    line += ';\n\n'
    file.write(line)
    
    # evaluate internal energy
    file.write('  // local array for species internal energies\n'
               '  Real u[{}];\n'.format(len(specs)) + 
               '  eval_u (y[0], u);\n'
               '\n'
               )
    
    # energy equation
    file.write('  // rate of change of temperature\n')
    line = '  dy[0] = (-1.0 / (rho * cv_avg)) * ( '
    isfirst = True
    for sp in specs:
        if len(line) > 70:
            line += '\n'
            file.write(line)
            line = '       '
        
        if not isfirst: line += ' + '
        
        isp = specs.index(sp)
        line += '(dy[{}] * u[{}] * {})'.format(isp + 1, isp, sp.mw)
        
        isfirst = False
    
    line += ' );\n\n'
    file.write(line)
    
    # rate of change of species mass fractions
    file.write('  // calculate rate of change of species mass fractions\n')
    for sp in specs:
        isp = specs.index(sp)
        line = '  dy[{}] *= ({} / rho);\n'.format(isp + 1, sp.mw)
        file.write(line)
    
    file.write('\n')
    file.write('} // end dydt\n\n')
    
    file.write('#endif\n')
    
    file.close()
    return

def write_mass_mole(path, lang, specs):
    """Writes files for mass/molar concentration and density conversion utility.
    
    Parameters
    ----------
    path : str
        Path to build directory for file.
    lang : {'c', 'cuda', 'fortran', 'matlab'}
        Programming language.
    specs : list of SpecInfo
        List of species in mechanism.
    
    Returns
    -------
    None
    
    """
    
    # Create header file
    if lang in ['c', 'cuda']:
        file = open(path + 'mass_mole.h', 'w')
    
        file.write('#ifndef MASS_MOLE_H\n'
                   '#define MASS_MOLE_H\n'
                   '\n'
                   '#include "header.h"\n'
                   '\n'
                   '#ifdef __cplusplus\n'
                   '  extern "C" {\n'
                   '#endif\n'
                   '\n'
                   'void mole2mass (const Real*, Real*);\n'
                   'void mass2mole (const Real*, Real*);\n'
                   'Real getDensity (const Real, const Real, const Real*);\n'
                   '\n'
                   '#ifdef __cplusplus\n'
                   '  }\n'
                   '#endif\n'
                   '#endif\n'
                   )
        file.close()
    
    # Open file; both C and CUDA programs use C file (only used on host)
    if lang in ['c', 'cuda']:
        filename = 'mass_mole.c'
    elif lang == 'fortran':
        filename = 'mass_mole.f90'
    elif lang == 'matlab':
        filename = 'mass_mole.m'
    file = open(path + filename, 'w')
    
    if lang in ['c', 'cuda']:
        file.write('#include "mass_mole.h"\n\n')
    
    ###################################################
    # Documentation and function/subroutine initializaton for mole2mass
    if lang in ['c', 'cuda']:
        file.write('/** Function converting species mole fractions to '
                   'mass fractions.\n'
                   ' *\n'
                   ' * \param[in]  X  array of species mole fractions\n'
                   ' * \param[out] Y  array of species mass fractions\n'
                   ' */\n'
                   'void mole2mass (const Real * X, Real * Y) {\n'
                   '\n'
                   )
    elif lang == 'fortran':
        file.write('!-----------------------------------------------------------------\n'
                   '!> Subroutine converting species mole fractions to mass fractions.\n'
                   '!! @param[in]  X  array of species mole fractions\n'
                   '!! @param[out] Y  array of species mass fractions\n'
                   '!-----------------------------------------------------------------\n'
                   'subroutine mole2mass (X, Y)\n'
                   '  implicit none\n'
                   '  double, dimension(:), intent(in) :: X\n'
                   '  double, dimension(:), intent(out) :: X\n'
                   '  double :: mw_avg\n'
                   '\n'
                   )
    
    # calculate molecular weight
    if lang in ['c', 'cuda']:
        file.write('  // average molecular weight\n'
                   '  Real mw_avg = 0.0;\n'
                   )
        for sp in specs:
            file.write('  mw_avg += X[{}] * '.format(specs.index(sp)) + 
                       '{};\n'.format(sp.mw)
                       )
    elif lang == 'fortran':
        file.write('  ! average molecular weight\n'
                   '  mw_avg = 0.0\n'
                   )
        for sp in specs:
            file.write('  mw_avg = mw_avg + '
                       'X({}) * '.format(specs.index(sp) + 1) + 
                       '{}\n'.format(sp.mw)
                       )
    file.write('\n')
    
    # calculate mass fractions
    if lang in ['c', 'cuda']:
        file.write('  // calculate mass fractions\n')
        for sp in specs:
            file.write('  Y[{0}] = X[{0}] * '.format(specs.index(sp)) + 
                       '{} / mw_avg;\n'.format(sp.mw)
                       )
        file.write('\n'
                   '} // end mole2mass\n'
                   '\n'
                   )
    elif lang == 'fortran':
        file.write('  ! calculate mass fractions\n')
        for sp in specs:
            file.write('  Y({0}) = X({0}) * '.format(specs.index(sp) + 1) + 
                       '{} / mw_avg\n'.format(sp.mw)
                       )
        file.write('\n'
                   'end subroutine mole2mass\n'
                   '\n'
                   )
    
    ################################
    # Documentation and function/subroutine initialization for mass2mole
    
    if lang in ['c', 'cuda']:
        file.write('/** Function converting species mass fractions to mole '
                   'fractions.\n'
                   ' *\n'
                   ' * \param[in]  Y  array of species mass fractions\n'
                   ' * \param[out] X  array of species mole fractions\n'
                   ' */\n'
                   'void mass2mole (const Real * Y, Real * X) {\n'
                   '\n'
                   )
    elif lang == 'fortran':
        file.write('!-------------------------------------------------------'
                   '----------\n'
                   '!> Subroutine converting species mass fractions to mole '
                   'fractions.\n'
                   '!! @param[in]  Y  array of species mass fractions\n'
                   '!! @param[out] X  array of species mole fractions\n'
                   '!-------------------------------------------------------'
                   '----------\n'
                   'subroutine mass2mole (Y, X)\n'
                   '  implicit none\n'
                   '  double, dimension(:), intent(in) :: Y\n'
                   '  double, dimension(:), intent(out) :: X\n'
                   '  double :: mw_avg\n'
                   '\n'
                   )
    
    # calculate average molecular weight
    if lang in ['c', 'cuda']:
        file.write('  // average molecular weight\n')
        file.write('  Real mw_avg = 0.0;\n')
        for sp in specs:
            file.write('  mw_avg += Y[{}] / '.format(specs.index(sp)) + 
                       '{};\n'.format(sp.mw)
                       )
        file.write('  mw_avg = 1.0 / mw_avg;\n')
    elif lang == 'fortran':
        file.write('  ! average molecular weight\n')
        file.write('  mw_avg = 0.0\n')
        for sp in specs:
            file.write('  mw_avg = mw_avg + '
                       'Y({}) / '.format(specs.index(sp) + 1) + 
                       '{}\n'.format(sp.mw)
                       )
    file.write('\n')
    
    # calculate mass fractions
    if lang in ['c', 'cuda']:
        file.write('  // calculate mass fractions\n')
        for sp in specs:
            file.write('  X[{0}] = Y[{0}] * '.format(specs.index(sp)) + 
                       'mw_avg / {};\n'.format(sp.mw)
                       )
        file.write('\n'
                   '} // end mass2mole\n'
                   '\n'
                   )
    elif lang == 'fortran':
        file.write('  ! calculate mass fractions\n')
        for sp in specs:
            file.write('  X({0}) = Y({0}) * '.format(specs.index(sp) + 1) + 
                       'mw_avg / {}\n'.format(sp.mw)
                       )
        file.write('\n'
                   'end subroutine mass2mole\n'
                   '\n'
                   )
    
    ###############################
    # Documentation and subroutine/function initialization for getDensity
    
    if lang in ['c', 'cuda']:
        file.write('/** Function calculating density from mole fractions.\n'
                   ' *\n'
                   ' * \param[in]  temp  temperature\n'
                   ' * \param[in]  pres  pressure\n'
                   ' * \param[in]  X     array of species mole fractions\n'
                   r' * \return     rho  mixture mass density' + '\n'
                   ' */\n'
                   'Real getDensity (const Real temp, const Real pres, '
                   'const Real * X) {\n'
                   '\n'
                   )
    elif lang == 'fortran':
        file.write('!-------------------------------------------------------'
                   '----------\n'
                   '!> Function calculating density from mole fractions.\n'
                   '!! @param[in]  temp  temperature\n'
                   '!! @param[in]  pres  pressure\n'
                   '!! @param[in]  X     array of species mole fractions\n'
                   '!! @return     rho   mixture mass density' + '\n'
                   '!-------------------------------------------------------'
                   '----------\n'
                   'function mass2mole (temp, pres, X) result(rho)\n'
                   '  implicit none\n'
                   '  double, intent(in) :: temp, pres\n'
                   '  double, dimension(:), intent(in) :: X\n'
                   '  double :: mw_avg, rho\n'
                   '\n'
                   )
    
    # get molecular weight
    if lang in ['c', 'cuda']:
        file.write('  // average molecular weight\n'
                   '  Real mw_avg = 0.0;\n'
                   )
        for sp in specs:
            file.write('  mw_avg += X[{}] * '.format(specs.index(sp)) + 
                      '{};\n'.format(sp.mw)
                      )
        file.write('\n')
    elif lang == 'fortran':
        file.write('  ! average molecular weight\n'
                   '  mw_avg = 0.0\n'
                   )
        for sp in specs:
            file.write('  mw_avg = mw_avg + '
                       'X({}) * '.format(specs.index(sp) + 1) + 
                       '{}\n'.format(sp.mw)
                       )
        file.write('\n')
    
    # calculate density
    if lang in ['c', 'cuda']:
        file.write('  return pres * mw_avg / ({:.8e} * temp);'.format(chem.RU))
        file.write('\n')
    else:
        line = '  rho = pres * mw_avg / ({:.8e} * temp)'.format(chem.RU)
        line += utils.line_end[lang]
        file.write(line)
    
    if lang in ['c', 'cuda']:
        file.write('} // end getDensity\n\n')
    elif lang == 'fortran':
        file.write('end function getDensity\n\n')
    
    file.close()
    return


def create_rate_subs(lang, mech_name, therm_name = None):
    """Create rate subroutines from mechanism.
    
    Parameters
    ----------
    lang : {'c', 'cuda', 'fortran', 'matlab'}
        Language type.
    mech_name : str
        Reaction mechanism filename (e.g. 'mech.dat').
    therm_name : str, optional
        Thermodynamic database filename (e.g. 'therm.dat') 
        or nothing if info in mechanism file.
    
    Returns
    -------
    None
    
    """
        
    lang = lang.lower()
    if lang not in utils.langs:
        print('Error: language needs to be one of: ')
        for l in utils.langs:
            print(l)
        sys.exit()
    
    # create output directory if none exists
    build_path = './out/'
    utils.create_dir(build_path)
    
    # interpret reaction mechanism file
    [elems, specs, reacs, units] = mech.read_mech(mech_name)
    
    # interpret thermodynamic database file (if it exists & needed)
    therm_flag = True
    if therm_name:
        # check for any species missing molecular weight
        for sp in specs:
            if not sp.mw:
                therm_flag = False
                break
        if not therm_flag:
            # need to read thermo file
            file = open(therm_name, 'r')
            mech.read_thermo(file, elems, specs)
            file.close()
    
    # convert activation energy units to K (if needed)
    if 'kelvin' not in units:
        efac = 1.0
        
        if 'kcal/mole' in units:
            efac = 4184.0 / chem.RU_JOUL
        elif 'cal/mole' in units:
            efac = 4.184 / chem.RU_JOUL
        elif 'kjoule' in units:
            efac = 1000.0 / chem.RU_JOUL
        elif 'joules' in units:
            efac = 1.00 / chem.RU_JOUL
        elif 'evolt' in units:
            efac = 11595.0
        else:
            # default is cal/mole
            efac = 4.184 / chem.RU_JOUL
        
        for rxn in reacs:
            rxn.E *= efac
        
        for rxn in [rxn for rxn in reacs if rxn.low]:
            rxn.low[2] *= efac
        
        for rxn in [rxn for rxn in reacs if rxn.high]:
            rxn.high[2] *= efac
    
    # now begin writing subroutines
    
    # print reaction rate subroutine
    write_rxn_rates(build_path, lang, specs, reacs)
    
    # if third-body/pressure-dependent reactions, 
    # print modification subroutine
    if next((r for r in reacs if (r.thd or r.pdep)), None):
        write_rxn_pressure_mod(build_path, lang, specs, reacs)
    
    # write species rates subroutine
    write_spec_rates(build_path, lang, specs, reacs)
    
    # write chem_utils subroutines
    write_chem_utils(build_path, lang, specs)
    
    # write derivative subroutines
    write_derivs(build_path, lang, specs, reacs)
    
    # write mass-mole fraction conversion subroutine
    write_mass_mole(build_path, lang, specs)
    
    return


if __name__ == "__main__":
    import argparse
    
    # command line arguments
    parser = argparse.ArgumentParser(description = 'Generates source code '
                                     'for species and reaction rates.'
                                     )
    parser.add_argument('-l', '--lang',
                        type = str,
                        choices = utils.langs, 
                        required = True, 
                        help = 'Programming language for output '
                        'source files.'
                        )
    parser.add_argument('-i', '--input',
                        type = str,
                        required = True, 
                        help = 'Input mechanism filename (e.g., mech.dat).'
                        )
    parser.add_argument('-t', '--thermo',
                        type = str,
                        default = None, 
                        help = 'Thermodynamic database filename (e.g., '
                        'therm.dat), or nothing if in mechanism.'
                        )
    
    args = parser.parse_args()
    
    create_rate_subs(args.lang, args.input, args.thermo)

