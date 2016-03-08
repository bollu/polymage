import sys
from polymage_smoother import w_jacobi
from polymage_defect import defect
from polymage_restrict import restrict
from polymage_interpolate import interpolate

sys.path.insert(0, '../../../../')

from compiler   import *
from constructs import *

def w_cycle(app_data):
    pipe_data = app_data['pipe_data']
    N = pipe_data['N']
    L = app_data['L']

    nu1 = app_data['nu1']
    nu2 = app_data['nu2']
    nuc = app_data['nuc']

    # initial guess
    V = Image(Double, "V_", [N[L]+2, N[L]+2])
    # rhs
    F = Image(Double, "F_", [N[L]+2, N[L]+2])

    jacobi_c = pipe_data['jacobi_c']

    # pre-smoothing and coarse-smoothing outputs
    smooth_p1 = {}
    # post-smoothing outputs
    smooth_p2 = {}

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
    def rec_w_cycle(v, f, l, visit):

        visit[l] += 1

        # coarsest level
        if l == 0:
            ''' COARSE-SMOOTHING '''
            if nuc == 0:
                return v

            if visit[l] == 1:
                smooth_p1[l] = {}

            smooth_p1[l][visit[l]] = {}

            for t in range(0, nuc):
                if l == L and t == nuc-1:
                    fname = "Wcycle"
                else:
                    fname = "T"+str(t)+"_coarse"+"__"+str(visit[l])

                if t == 0:
                    in_func = v
                else:
                    in_func = smooth_p1[l][visit[l]][t-1]

                smooth_p1[l][visit[l]][t] = \
                    w_jacobi(in_func, f, l, fname, app_data)

            return smooth_p1[l][visit[l]][nuc-1]
        ###################################################
        # all other finer levels
        else:
            ''' PRE-SMOOTHING '''
            if visit[l] == 1:
                smooth_p1[l] = {}

            smooth_p1[l][visit[l]] = {}

            for t in range(0, nu1):
                fname = "T"+str(t)+"_pre_L"+str(l)+"__"+str(visit[l])
                if t == 0:
                    in_func = v
                else:
                    in_func = smooth_p1[l][visit[l]][t-1]

                smooth_p1[l][visit[l]][t] = \
                    w_jacobi(in_func, f, l, fname, app_data)

            if nu1 <= 0:
                smooth_out = v
            else:
                smooth_out = smooth_p1[l][visit[l]][nu1-1]

            ###############################################
            ''' RESIDUAL '''

            if visit[l] == 1:
                r_h[l] = {}

            name = "defect_L"+str(l)+"__"+str(visit[l])
            r_h[l][visit[l]] = defect(smooth_out, f, l, name, pipe_data)

            ###############################################
  
            ''' RESTRICTION '''
            if visit[l] == 1:
                r_2h[l] = {}

            name = "restrict_L"+str(l-1)+"__"+str(visit[l])
            r_2h[l][visit[l]] = restrict(r_h[l][visit[l]], l, name, pipe_data)

            ###############################################
 
            ''''''''''''''''''
            ''' NEXT LEVEL '''
            ''''''''''''''''''

            if visit[l] == 1:
                e_2h[l] = {}

            # e_2h <- 0
            e_2h[l][visit[l]] = \
                rec_w_cycle(None, r_2h[l][visit[l]], l-1, visit)

            e_2h[l][visit[l]] = \
                rec_w_cycle(e_2h[l][visit[l]], r_2h[l][visit[l]], l-1, visit)

            ###############################################

            ''' INTERPOLATION & CORRECTION '''
            if l == L and nu2 <= 0:
                fname = "Wcycle"
            else:
                fname = "interp_correct_L"+str(l)+"__"+str(visit[l])

            if nu1 <= 0:
                correct_in = v
            else:
                correct_in = smooth_p1[l][visit[l]][nu1-1]

            if visit[l] == 1:
                ec[l] = {}

            ec[l][visit[l]] = \
                interpolate(e_2h[l][visit[l]], correct_in, l, fname, pipe_data)

            if nu2 <= 0:
                return ec[l][visit[l]]
 
            ###############################################

            ''' POST-SMOOTHING '''

            if visit[l] == 1:
                smooth_p2[l] = {}

            smooth_p2[l][visit[l]] = {}
            for t in range(0, nu2):
                fname = "T"+str(t)+"_post_L"+str(l)+"__"+str(visit[l])
                if l == L and t == nu2-1:
                    fname = "Wcycle"

                if t == 0:
                    in_func = ec[l][visit[l]]
                else:
                    in_func = smooth_p2[l][visit[l]][t-1]

                smooth_p2[l][visit[l]][t] = \
                    w_jacobi(in_func, f, l, fname, app_data)
 
            return smooth_p2[l][visit[l]][nu2-1]
    #######################################################

    visit = {}
    visit[L] = 0
    for l in range(0, L):
        visit[l] = 0

    # one whole v-cycle beginning at the finest level
    u = rec_w_cycle(V, F, L, visit)

    return u
