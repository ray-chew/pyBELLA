import numpy as np
import logging

from . import ensemble_access


def da_interface(results, obs, obs_attributes, delta, times, tout, N):
    igx, igy = results[0].elem.igx, results[0].elem.igy
    inner = (slice(igx, -igx), slice(igy, -igy))

    attributes = ["rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"]
    # attributes = ['rho', 'rhou', 'rhov']
    tmp = ensemble_access.stack_fields(results, obs_attributes[0], slc=inner)
    tmp = tmp[:, np.newaxis, ...]

    ensemble = ensemble_access.stack_fields(results, attributes[0], slc=inner)
    ensemble = ensemble[:, np.newaxis, ...]

    obs_current = np.array(
        obs[np.where(np.isclose(times, tout))[0][0]][obs_attributes[0]]
    )[inner]

    # obs_current = bin_func(obs_current,(Nx,Ny))
    obs_current = obs_current[np.newaxis, ...]

    for attr in obs_attributes[1:]:
        tmp = np.hstack(
            (
                tmp,
                ensemble_access.stack_fields(results, attr, slc=inner)[
                    :, np.newaxis, ...
                ],
            )
        )
        tmp01 = np.array(obs[np.where(np.isclose(times, tout))[0][0]][attr])[inner]
        # tmp01 = bin_func(tmp01,(Nx,Ny))
        # print(tmp01.shape)
        tmp01 = tmp01[np.newaxis, ...]
        obs_current = np.vstack((obs_current, tmp01))

    for attr in attributes[1:]:
        ensemble = np.hstack(
            (
                ensemble,
                ensemble_access.stack_fields(results, attr, slc=inner)[
                    :, np.newaxis, ...
                ],
            )
        )

    obs_current = obs_current[np.newaxis, ...]
    # r = tmp - obs_current

    Hx = tmp.reshape(N, -1)
    obs_current = obs_current.reshape(-1)
    # ensemble = ensemble.reshape(N,-1)
    # print(ensemble.shape)

    etpf = analysis(ensemble, delta)
    etpf.analyse(obs_current, 1.0, Hx, N)

    analysis_ens = etpf.get_ensemble_from_X()

    for n in range(N):
        cnt = 0
        current = analysis_ens[n]
        for attr in attributes:
            data = current[cnt]
            data = np.pad(data, ((igx, igx), (igy, igy)), mode="wrap")

            ensemble_access.set_field(results[n], attr, data)
            cnt += 1
    return results


class analysis(object):
    def __init__(self, ensemble, delta, identifier=None):
        self.ensemble = np.array(ensemble)
        self.ensemble_shape = self.ensemble.shape

        self.X = self.state_vector(ensemble)
        self.identifier = identifier

        # rejuvenation factor
        self.delta = delta

    def analyse(self, obs_current, obs_covar, Hx, N):
        try:
            import ot
        except ImportError as e:
            raise ImportError(
                "the ETPF needs the POT package: pip install 'pybella[da]'"
            ) from e

        logging.info("starting ETPF analysis...")

        r = (Hx - obs_current) ** 2
        r = np.sum(r, axis=1)
        # print(r.shape)

        ww = np.exp(-r / (2.0 * obs_covar))
        ww /= np.sum(ww)

        # print(ww)

        Co = self.X @ self.X.T
        diag = np.diag(Co)
        Co = diag * np.ones((1, N)) - 2.0 * Co + np.ones((N, 1)) * diag.T

        # print(Co)

        # exact optimal-transport plan; equivalent to the reference
        # pyemd.emd_with_flow call (both histograms sum to one, so the
        # extra-mass penalty there was irrelevant)
        T = ot.emd(ww, np.ones(N) / N, Co)
        T = T * N

        self.X = np.dot(self.X.T, T).T + self.delta * np.random.randn(
            self.X.shape[0], self.X.shape[1]
        )
        # print(self.X.shape)

    @staticmethod
    def state_vector(ensemble):
        return ensemble.reshape(ensemble.shape[0], -1)

    def get_ensemble_from_X(self):
        return self.X.reshape(self.ensemble_shape)
