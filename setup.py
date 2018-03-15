from cx_Freeze import setup, Executable

base = None


executables = [Executable("run_meas.py", base=base)]

packages = ["idna"]
options = {
    'build_exe': {

        'packages':packages,
    },

}

setup(
    name = "TheMeasurator",
    options = options,
    version = "0.3",
    description = 'InDevelopment',
    executables = executables
)
