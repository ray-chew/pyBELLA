import logging


######################################################
# Nonhydrostatic - Hydrostatic blending
######################################################
def do_nonhydro_to_hydro_conv(
    sol, flux, npf, bld, elem, node, th, ud, label, writer, step, window_step, t, dt
):
    logging.info("nonhydrostatic to hydrostatic conversion...")
    # bld.convert_p2n(npf.p2_nodes)
    # bld.update_sol(sol,elem,node,th,ud,npf,'bef',label=label,writer=writer)
    # sol.rhov = sol.rhov_half

    # sol_tmp = deepcopy(sol)
    # flux_tmp = deepcopy(flux)
    # npf_tmp = deepcopy(npf)

    # nonhydro to hydro blending incomplete.
    # ret = data.time_update(sol,flux,npf, t, t+1*dt, ud, elem, node, [0,0], th, bld=None, writer=None, debug=False)

    # sol = sol_tmp
    # flux = flux_tmp
    # npf = npf_tmp
    # sol = ret[0]
    # flux = ret[1]
    # npf = ret[2]
    # sol = deepcopy(ret[0])
    # npf = deepcopy(ret[2])
    # sol.rhov[...] = sol.rhov_half
    # t += 0.5*dt
    # t += 1*dt
    return sol, npf, t


def do_hydro_to_nonhydro_conv(
    sol, flux, npf, bld, elem, node, th, ud, label, writer, step, window_step, t, dt
):
    logging.info("hydrostatic to nonhydrostatic conversion...")
    logging.info(f"Blending... step = {step}")

    # sol_tmp = deepcopy(sol)
    # flux_tmp = deepcopy(flux)
    # npf_tmp = deepcopy(npf)

    # ret = data.time_update(sol,flux,npf, t, t+dt, ud, elem, node, [0,step-1], th, bld=None, writer=None, debug=False)

    # sol = sol_tmp
    # flux = flux_tmp
    # npf = npf_tmp

    # retv_half = ret[0].rhov_half / ret[0].rho_half
    # retv_full = ret[0].rhov / ret[0].rho

    # solv_half = sol.rhov_half / sol.rho_half
    # solv_full = sol.rhov / sol.rho

    # fac_full = 0.5
    # fac_half = 1.0 - fac_full

    # # logging.info(np.sum(solv_full))
    # # logging.info(np.sum(retv_half))
    # # logging.info(np.sum((fac_full * solv_full + fac_half * retv_half)))
    # # logging.info(np.sum(fac_half * retv_half))

    # fac_full = 0.5
    # fac_half = 0.5

    # # logging.info(np.sum(solv_full))
    # # logging.info(np.sum(retv_half))
    # # logging.info(np.sum((fac_full * solv_full + fac_half * retv_half)))
    # # logging.info(np.sum(fac_half * retv_half))

    # if writer != None: writer.populate(str(label)+'_after_full_step', 'ret_half', ret[0].rhov_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'ret_full', ret[0].rhov)

    # if writer != None: writer.populate(str(label)+'_after_full_step', 'solv_half', sol.rhov_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'solv_full', sol.rhov)

    # sol.rhov = sol.rho * (fac_full * solv_full + fac_half * retv_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_end', ret[2].p2_nodes)

    # fac_npf_half = 0.5
    # fac_npf_full = 1.0 - fac_npf_half
    # npf.p2_nodes = fac_npf_half * npf.p2_nodes + fac_npf_full * ret[2].p2_nodes
    # dp2n = ret[2].p2_nodes_half
    # bld.convert_p2n(dp2n)
    # bld.update_sol(sol,elem,node,th,ud,npf,'aft',label=label,writer=writer)
    # bld.update_p2n(sol,npf,node,th,ud)
    #

    ###############################
    # alternative version
    ###############################

    # if c1 or c2:
    #     logging.info(
    #         termcolor.colored("hydrostatic to nonhydrostatic conversion...", "blue")
    #     )

    # writer.write_all(mem, str(label) + "_half_full")
    # writer.populate(str(label) + "_ic", "pwchi", sol.pwchi)

    # if test_hydrob == False:
    #     sol = copy.deepcopy(sol_half_old)
    #     # npf = copy.deepcopy(npf_half_old)

    #     logging.info(termcolor.colored("test_hydrob == False", "red"))
    #     writer.write_all(mem, str(label) + "_quarter")

    #     writer.populate(str(label) + "_quarter", "pwchi", sol.pwchi)

    #     logging.info("quarter dt = %.8f" % (dt * 0.5))

    #     ret = do(
    #         sol_half_old,
    #         flux_half_old,
    #         npf_half_old,
    #         dt - 0.5 * dt,
    #         dt + 0.5 * dt,
    #         ud,
    #         elem,
    #         node,
    #         [0, 0],
    #         th,
    #         bld=None,
    #         writer=None,
    #         debug=False,
    #     )

    #     sol_tu = copy.deepcopy(ret[0])
    #     # npf_tu = copy.deepcopy(ret[2])
    #     sol.rho[...] = sol_tu.rho_half
    #     sol.rhou[...] = sol_tu.rhou_half
    #     sol.rhov[...] = sol_tu.rhov_half
    #     sol.rhow[...] = sol_tu.rhow_half
    #     sol.rhoX[...] = sol_tu.rhoX_half
    #     sol.rhoY[...] = sol_tu.rhoY_half
    #     sol.pwchi[...] = sol_tu.pwchi

    #     # npf.p2_nodes[...] = npf_tu.p2_nodes_half

    #     writer.write_all(mem, str(label) + "_half")

    #     writer.populate(str(label) + "_half", "pwchi", sol.pwchi)

    #     ret = do(
    #         sol,
    #         flux,
    #         npf,
    #         dt,
    #         2.0 * dt,
    #         ud,
    #         elem,
    #         node,
    #         [0, 0],
    #         th,
    #         bld=None,
    #         writer=None,
    #         debug=False,
    #     )

    #     sol = copy.deepcopy(ret[0])
    #     flux = copy.deepcopy(ret[1])
    #     npf = copy.deepcopy(ret[2])

    # if test_hydrob == True:
    #     sol = copy.deepcopy(sol_half_old)
    #     # npf = copy.deepcopy(npf_half_old)

    #     logging.info(termcolor.colored("test_hydrob == False", "red"))
    #     writer.write_all(mem, str(label) + "_quarter")

    #     # writer.populate(str(label)+'_quarter', 'pwchi', sol.pwchi)

    #     logging.info("quarter dt = %.8f" % (dt * 0.5))

    #     ret = do(
    #         sol_half_old,
    #         flux_half_old,
    #         npf_half_old,
    #         dt - 0.5 * dt,
    #         dt + 0.5 * dt,
    #         ud,
    #         elem,
    #         node,
    #         [0, 0],
    #         th,
    #         bld=None,
    #         writer=None,
    #         debug=False,
    #     )

    #     sol_tu = copy.deepcopy(ret[0])
    #     # npf_tu = copy.deepcopy(ret[2])
    #     sol.rho[...] = sol_tu.rho_half
    #     sol.rhou[...] = sol_tu.rhou_half
    #     sol.rhov[...] = sol_tu.rhov_half
    #     sol.rhow[...] = sol_tu.rhow_half
    #     sol.rhoX[...] = sol_tu.rhoX_half
    #     sol.rhoY[...] = sol_tu.rhoY_half
    #     sol.pwchi[...] = sol_tu.pwchi

    #     # npf.p2_nodes[...] = npf_tu.p2_nodes_half

    #     # writer.write_all(sol,npf,elem,node,th,str(label)+'_half')

    #     # writer.populate(str(label)+'_half', 'pwchi', sol.pwchi)

    #     ret = do(
    #         sol,
    #         flux,
    #         npf,
    #         dt,
    #         2.0 * dt,
    #         ud,
    #         elem,
    #         node,
    #         [0, 0],
    #         th,
    #         bld=None,
    #         writer=None,
    #         debug=False,
    #     )

    #     sol = copy.deepcopy(ret[0])
    #     flux = copy.deepcopy(ret[1])
    #     npf = copy.deepcopy(ret[2])
    #     # writer.write_all(sol,npf,elem,node,th,str(label)+'_half')
    #     # writer.populate(str(label)+'_half', 'pwchi', sol.pwchi)

    #     logging.info(termcolor.colored("test_hydrob == True", "red"))

    # if test_hydrob == False:
    #     dt *= 2.0
    # if c2:
    # ud.is_nonhydrostatic = 1

    return sol, npf
