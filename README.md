# The Kilonova Code (KNC)
_Simulates relativistic jets, winds, and explosions in 2D axisymmety._

## Quick Start

Copy + paste these lines at your terminal to clone and build the code:
```bash
git clone https://github.com/clemson-cal/app-kilonova
cd app-kilonova
cargo build --release
```
From within the `app-kilonova`, copy + paste these lines to build and install the Python loader extension:
```bash
pip3 install maturin
cd knc_loader
maturin build
pip3 install target/wheels/*.whl
```

To remove the Python module at a later time, you can do `pip3 uninstall knc_loader`. You can also install the code locally by running `cargo install --path .` from the project root directory. This will make an executable called `kilonova` available in your path (it will be in the `~/.cargo/bin` directory). To run the code and generate a plot, you can use one of the preset configurations. For example, to run the `jet_in_cloud` problem for 0.1 seconds, type the following:
```bash
kilonova jet_in_cloud control.final_time=1.1
```
This will write two "products" files to the current working directory, `prods.0000.cbor` and `prods.0001.cbor`, and then exit. The first one is the solution at the simulation start time (1 second for this setup) and at the end time at 1.1 seconds. To plot either of these files, you can use the inclded plotting script:
```bash
python3 knc_tools/plot.py prods.0001.cbor --field=ur
```
This will show a relief plot of the gas radial four-velocity. To see more plotting options, run `python3 knc_tools/plot.py --help`.

## Developers
KNC is written and maintained by the [Computational Astrophysics Lab](https://jzrake.people.clemson.edu) at the [Clemson University Department of Physics and Astronomy](http://www.clemson.edu/science/departments/physics-astro). The core developer/maintainer is presently Jonathan Zrake.
