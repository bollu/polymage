import sys
from polymage_smoother    import wJacobi
from polymage_defect      import defect
from polymage_restrict    import restrict
from polymage_interpolate import interpolate

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def wCycle(pipeData, appData):
    N = pipeData['N']
    L = appData['L']

    nu1 = appData['nu1']
    nu2 = appData['nu2']
    nuc = appData['nuc']

    # initial guess
    V = Image(Double, "V_", [N[L]+2, N[L]+2, N[L]+2])
    # rhs
    F = Image(Double, "F_", [N[L]+2, N[L]+2, N[L]+2])

    jacobi_c = pipeData['jacobi_c']

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
    def recWCycle(v, f, l, visit):

        visit[l] += 1

        # coarsest level
        if l == 0:
            ''' COARSE-SMOOTHING '''
            if nuc == 0:
                return v

            if visit[l] == 1:
                smoothP1[l] = {}


            smoothP1[l][visit[l]] = {}

            for t in range(0, nuc):
                if l == L and t == nuc-1:
                    fname = "Wcycle"
                else:
                    fname = "T"+str(t)+"_coarse"+"__"+str(visit[l])

                if t == 0:
                    inFunc = v
                else:
                    inFunc = smoothP1[l][visit[l]][t-1]

                smoothP1[l][visit[l]][t] = wJacobi(inFunc, f, l, fname,
                                                   pipeData, appData)

            return smoothP1[l][visit[l]][nuc-1]
        ###################################################
        # all other finer levels
        else:
            ''' PRE-SMOOTHING '''
            if visit[l] == 1:
                smoothP1[l] = {}

            smoothP1[l][visit[l]] = {}

            for t in range(0, nu1):
                fname = "T"+str(t)+"_pre_L"+str(l)+"__"+str(visit[l])
                if t == 0:
                    inFunc = v
                else:
                    inFunc = smoothP1[l][visit[l]][t-1]

                smoothP1[l][visit[l]][t] = wJacobi(inFunc, f, l, fname,
                                                   pipeData, appData)

            if nu1 <= 0:
                smoothOut = v
            else:
                smoothOut = smoothP1[l][visit[l]][nu1-1]

            ###############################################
            ''' RESIDUAL '''

            if visit[l] == 1:
                r_h[l] = {}

            r_h[l][visit[l]] = defect(smoothOut, f, l, 
                                      "defect_L"+str(l)+"__"+str(visit[l]),
                                      pipeData)

            ###############################################
  
            ''' RESTRICTION '''
            if visit[l] == 1:
                r_2h[l] = {}

            r_2h[l][visit[l]] = restrict(r_h[l][visit[l]], l, 
                                         "restrict_L"+str(l-1)+"__"+str(visit[l]),
                                         pipeData)

            ###############################################
 
            ''''''''''''''''''
            ''' NEXT LEVEL '''
            ''''''''''''''''''

            if visit[l] == 1:
                e_2h[l] = {}

            # e_2h <- 0
            e_2h[l][visit[l]] = recWCycle(None, 
                                          r_2h[l][visit[l]],
                                          l-1, visit)

            e_2h[l][visit[l]] = recWCycle(e_2h[l][visit[l]],
                                          r_2h[l][visit[l]],
                                          l-1, visit)

            ###############################################

            ''' INTERPOLATION & CORRECTION '''
            if l == L and nu2 <= 0:
                fname = "Wcycle"
            else:
                fname = "interp_correct_L"+str(l)+"__"+str(visit[l])

            if nu1 <= 0:
                correctIn = v
            else:
                correctIn = smoothP1[l][visit[l]][nu1-1]

            if visit[l] == 1:
                ec[l] = {}

            ec[l][visit[l]] = interpolate(e_2h[l][visit[l]],
                                          correctIn, l, fname,
                                          pipeData)

            if nu2 <= 0:
                return ec[l][visit[l]]
 
            ###############################################

            ''' POST-SMOOTHING '''

            if visit[l] == 1:
                smoothP2[l] = {}

            smoothP2[l][visit[l]] = {}
            for t in range(0, nu2):
                fname = "T"+str(t)+"_post_L"+str(l)+"__"+str(visit[l])
                if l == L and t == nu2-1:
                    fname = "Wcycle"

                if t == 0:
                    inFunc = ec[l][visit[l]]
                else:
                    inFunc = smoothP2[l][visit[l]][t-1]

                smoothP2[l][visit[l]][t] = wJacobi(inFunc, f, l, fname,
                                                   pipeData, appData)
 
            return smoothP2[l][visit[l]][nu2-1]
    #######################################################

    visit = {}
    visit[L] = 0
    for l in range(0, L):
        visit[l] = 0

    # one whole v-cycle beginning at the finest level
    u = recWCycle(V, F, L, visit)

    return u
