# import numpy as np
# from management.enumerator import TimeIntegrator

# def SetTimeIntegratorParameters(ud):
#     if ud.advec_time_integrator == TimeIntegrator.EXPL_MIDPT:
#         ud.tips.dt_frac = 0.5
#         ud.tips.flux_frac[0][0] = 0.0
#         ud.tips.flux_frac[0][1] = 1.0
#         ud.tips.flux_frac[1][0] = -0.5
#         ud.tips.flux_frac[1][1] = 1.0

#         for k in range(2,ud.tips.NO_OF_RK_STAGES):
#             ud.tips.flux_frac[k][1] = 1e10

#         ud.tips.update_frac[0] = 0.5
#         ud.tips.update_frac[1] = 1.0

#         for k in range(2,ud.tips.NO_OF_RK_STAGES):
#             ud.tips.update_frac[k] = 1e10
        
#         ud.tips.multiD_updt = True
    
#     elif ud.advec_time_integrator == TimeIntegrator.HUEN:
#         ud.tips.dt_frac = 0.5
#         ud.tips.flux_frac[0][0] = 0.0
#         ud.tips.flux_frac[0][1] = 1.0
#         ud.tips.flux_frac[1][0] = -0.5
#         ud.tips.flux_frac[1][1] = 0.5

#         for k in range(ud.tips.NO_OF_RK_STAGES):
#             ud.tips.flux_frac[k,:] = 1e10

#         ud.tips.update_frac[0] = 1.0
#         ud.tips.update_frac[1] = 1.0

#         for k in range(ud.tips.NO_OF_RK_STAGES):
#             ud.tips.update_frac[k] = 1e10
        
#         ud.tips.multiD_updt = True
        
#     else:
#         ud.tips.dt_frac = 0.5
#         ud.tips.flux_frac[0][0] = 0.0
#         ud.tips.flux_frac[0][1] = 1.0
#         ud.tips.flux_frac[1][0] = 0.5
#         ud.tips.flux_frac[1][1] = 0.5

#         for k in range(ud.tips.NO_OF_RK_STAGES):
#             ud.tips.flux_frac[k,:] = 1e10

#         ud.tips.update_frac[0] = 0.5
#         ud.tips.update_frac[1] = 1.0

#         for k in range(ud.tips.NO_OF_RK_STAGES):
#             ud.tips.update_frac[k] = 1e10
        
#         ud.tips.multiD_updt = False
