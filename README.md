<p align="center">
  <a href="">
  <img alt="pyBELLA Logo" src="/docs/source/_static/logo.png">
  </a>
</p>

<h2 align="center"><b><font color="#417b95">B</font></b>lended s<b><font color="#417b95">E</font></b>am<b><font color="#417b95">L</font></b>ess so<b><font color="#417b95">L</font></b>ver for <b><font color="#417b95">A</font></b>tmospheric dynamics</h2>


<p align="center">
<a href="https://github.com/ray-chew/pyBELLA/actions/workflows/documentation.yml">
<img alt="GitHub Actions: docs" src=https://github.com/ray-chew/pyBELLA/actions/workflows/documentation.yml/badge.svg>
</a>
<a href="https://github.com/ray-chew/pyBELLA/issues">
<img alt="open issues" src=https://img.shields.io/github/issues/ray-chew/pyBELLA>
</a>
<a href="https://opensource.org/licenses/BSD-3-Clause">
<img alt="License: BSD-3" src=https://img.shields.io/badge/License-BSD_3--Clause-blue.svg>
</a>
<!-- <a href="https://github.com/psf/black">
<img alt="Code style: black" src=https://img.shields.io/badge/code%20style-black-000000.svg>
</a> -->
</p>


The Blended sEamLess soLver for Atmospheric dynamics (pyBELLA) is a Python-based numerical flow solver for atmospheric dynamics. The current version features PyBELLA+ as it is coupled to an ensemble data assimillation engine based on the Local Ensemble Transform Kalman Filter.

The numerical scheme for pyBELLA was introduced by [Bennachio and Klein (2019)](https://journals.ametsoc.org/view/journals/mwre/147/11/mwr-d-19-0073.1.xml) and the seamless blending between model regimes within a simulation run was extended in [Chew et al. (2022)](https://journals.ametsoc.org/view/journals/mwre/150/9/MWR-D-21-0175.1.xml).

This code was used to produce the results in [Chew (2022)](https://refubium.fu-berlin.de/bitstream/handle/fub188/37313/thesis_final.pdf?sequence=1&isAllowed=y) and [Chew et al. (2023)](https://tinyurl.com/2dc7hjqa).


---

**[Read the documentation here](https://ray-chew.github.io/pyBELLA/index.html)**

---

## Requirements

See [`requirements.txt`](https://github.com/ray-chew/pyBELLA/blob/main/requirements.txt)

**Note**:  The development dependencies can be found in [`dev-requirements.py`](https://github.com/ray-chew/pyBELLA/blob/main/dev-requirements.txt).


## Usage

### Installation

Fork this repository and clone your remote fork. `cd` into your local forked repository and execute:

```console
pip install -e . 
```

**Note**: depending on your IDE, you may need to add the `--config-settings editable_mode` at the end of the above command for the pyBELLA package to be visible to the linter. However, this comes with limitations to the development mode, see [here for more details](https://setuptools.pypa.io/en/latest/userguide/development_mode.html).

### Configuration

The user-defined input parameters are in the [`inputs`](https://github.com/ray-chew/pyBELLA/tree/main/inputs) subpackage. These parameters are imported into the run scripts in [`run_scripts`](https://github.com/ray-chew/pyBELLA/tree/main/run_scripts). 

### Execution

A simple test can be found in [`run_scripts.test_suite`](https://github.com/ray-chew/pyBELLA/blob/main/ru_scripts/test_suite.py). To execute this run script from the pyBELLA parent directory:

```console
python3 ./run_scripts/test_suite.py
```

However, the codebase is structured such that the user can easily assemble a run script to define their own experiments. Refer to the documentation for the [available APIs](https://ray-chew.github.io/pyBELLA/apis.html).

## License

[BSD License v3.0](https://fossa.com/blog/open-source-software-licenses-101-bsd-3-clause-license/)

## Contributions

Refer to the [open issues](https://github.com/ray-chew/pyBELLA/issues) that require attention.

Any changes, improvements, or bug fixes can be submitted to from your remote to upstream via a pull request.

