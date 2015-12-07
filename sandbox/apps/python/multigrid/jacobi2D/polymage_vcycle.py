import sys
from polymage_smoother    import wJacobi
from polymage_defect      import defect
from polymage_restrict    import restrict
from polymage_interpolate import interpolate

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def vCycle(impipeData, appData):
    N = impipeData['N']
    L = appData['L']

    nu1 = appData['nu1']
    nu2 = appData['nu2']
    nuc = appData['nuc']

    # initial guess
    V = Image(Double, "V_", [N[L]+2, N[L]+2])
    # rhs
    F = Image(Double, "F_", [N[L]+2, N[L]+2])

    jacobi_c = impipeData['jacobi_c']

    # pre-smoothing and coarse-smoothing outputs
    smoothP1 = {}
    # post-smoothing outputs
    smoothP2 = {}

    # defect outputs
    r_h = {}
    # restrict outputs
    r_2h = {}
    # error outputs
    e_h = {}
    # interpolated error outputs
    e_2h = {}
    # corrected error outputs
    ec = {}

    #######################################################
    def recVCycle(v, f, l):

        # coarsest level
        if l == 0:
            ''' COARSE-SMOOTHING '''
            if nuc == 0:
                return v

            smoothP1[l] = {}
            for t in range(0, nuc):
                if l == L and t == nuc-1:
                    fname = appData['cycle_name']
                else:
                    fname = "T"+str(t)+"_coarse"

                if t == 0:
                    inFunc = v
                else:
                    inFunc = smoothP1[l][t-1]

                smoothP1[l][t] = wJacobi(inFunc, f, l, fname,
                                         impipeData, appData)

            return smoothP1[l][nuc-1]
        ###################################################
        # all other finer levels
        else:
            ''' PRE-SMOOTHING '''
            smoothP1[l] = {}
            for t in range(0, nu1):
                fname = "T"+str(t)+"_pre_L"+str(l)
                if t == 0:
                    inFunc = v
                else:
                    inFunc = smoothP1[l][t-1]

                smoothP1[l][t] = wJacobi(inFunc, f, l, fname,
                                         impipeData, appData)

            if nu1 <= 0:
                smoothOut = v
            else:
                smoothOut = smoothP1[l][nu1-1]

            ###############################################
            ''' RESIDUAL '''

            r_h[l] = defect(smoothOut, f, l, "defect_L"+str(l),
                            impipeData)

            ###############################################
  
            ''' RESTRICTION '''
            r_2h[l] = restrict(r_h[l], l, "restrict_L"+str(l-1),
                               impipeData)

            ###############################################
 
            ''''''''''''''''''
            ''' NEXT LEVEL '''
            ''''''''''''''''''
            # e_2h <- 0
            e_2h[l] = recVCycle(None, r_2h[l], l-1)

            ###############################################

            ''' INTERPOLATION & CORRECTION '''
            if l == L and nu2 <= 0:
                fname = appData['cycle_name']
            else:
                fname = "interp_correct_L"+str(l)

            if nu1 <= 0:
                correctIn = v
            else:
                correctIn = smoothP1[l][nu1-1]

            ec[l] = interpolate(e_2h[l], correctIn, l, fname,
                                impipeData)

            if nu2 <= 0:
                return ec[l]
 
            ###############################################

            ''' POST-SMOOTHING '''
            smoothP2[l] = {}
            for t in range(0, nu2):
                fname = "T"+str(t)+"_post_L"+str(l)
                if l == L and t == nu2-1:
                    fname = appData['cycle_name']

                if t == 0:
                    inFunc = ec[l]
                else:
                    inFunc = smoothP2[l][t-1]

                smoothP2[l][t] = wJacobi(inFunc, f, l, fname,
                                         impipeData, appData)
 
            return smoothP2[l][nu2-1]
    #######################################################

    # one whole v-cycle beginning at the finest level
    u = recVCycle(V, F, L)

    return u
