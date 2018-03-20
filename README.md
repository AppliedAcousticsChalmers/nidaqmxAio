# nidaqmxAio 

Input/output acquisition automation using the [nidaqmx](https://github.com/ni/nidaqmx-python) python API and analog input/output modules.

### Prerequisites

You need to be able to run python scripts on your system with version >3.6 installed. nidaqmx supports only Windows OS.
Python prerequisites found in `requirements.txt`. To install them, run `pip install -r requirements.txt`.

### Usage

Run `python run_meas.py <simulation time> <sample rate>` with additional arguments for channel input/outputs. 

The script is built without external triggering/clocking devices in mind. Due to system delays that cannot be deterministically calculated, the test signals are zero padded to avoid truncated ends. Type `python run_meas.py --help` for information on the additional arguments. Measurements are saved in the directory `acquired_data\measurements_<YYMMDD>\` in the format `meas_<YYMMDD>.py` by specifying the desired filename, using `-sv filename`. During measurements, a live post-processing preview is displayed, presenting the spectra of reference and measured channel, transfer function (using the H1 estimator), impulse response and coherence. The reference channel can be selected from a list after running the script (leave blank to select the first channel on the list).

#### Input arguments

If you need to run the script with other than the default values, specify the corresponding arguments explicitly.

`-tp <string>` or `--signalType <string>`:
The type of the output signal to be generated. Currently, `white_noise` and `pink_noise` are supported (default: `pink_noise`)

`-pad <int>` or `--pad_samples <int>`:
Pad N samples at the end of the signal to avoid input truncation due to system delays (default `5000`)

`-ai <int> <int> ...` or `--channelsIn <int> <int> ...`:
Set the number of analog input channels per module. Begin from the lowest numbered slot and continue. The number of input channels per module are chosen serially. For example, if the argument is set as `-ai 2 4 0`, the first 2 input channels are used from the module connected at the lowest numbered slot on the chassis, 4 from the second and none from the third. The number of channels and modules set here should match your physical setup (default: `0`)

`-ao <int> <int> ...` or `--channelsOut <int> <int> ...`:
Set the number of analog output channels per module. Begin from the lowest numbered slot and continue. The number of output channels per module are chosen serially. For example, if the argument is set as `-ao 2 4 0`, the first 2 output channels are used from the module connected at the lowest numbered slot on the chassis, 4 from the second and none from the third. The number of channels and modules set here should match your physical setup. Note that the test signal is copied to each analog output channel. This is due to change and be configurable by the user in future updates (default: `0`)

`-aip <float>` or `--aiPeak <float>`:
The analog input absolute peak expected in volts (default: `5.0`)

`-aop <float>` or `--aoPeak <float>`:
The analog output absolute peak expected in volts (default: `1.0`)

`-sv <string>` or `--save_file <string>`:
If specified the collected data are saved, in the given filename.  (default: `no`)

`-bf <int>` or `--bufferSize <int>`:
Buffer size in samples of the input channels. This is used also as the time window for the analysis done at both the live post-processing preview and the Transfer function post-processing script. (default: `8192`)

`-cal <string>`or `--calibration <string>`:
If set as `new`, the calibration measurement script is run and the calibration coefficient and microphone sensitivity are saved. You will be asked to specify the name of the saved file (leave blank for the default name `calibration_cal`) and the microphone sensitivity. If not defined, the sensitivity used is for the microphone B&K type 4190.

`-pp <string>` or `--postProcess <string>`:
Run the chosen post processing script. You will be asked to select the directory and a the file from a list (not selecting a file will cause the script to run for all the files in the selected directory with the file extension`.np[yz]`, except the calibration files). Results will be saved as `<measurement_name>_TFs`. Current choices:
`no` - don't run post process scripts
`TF` - run transfer function process script using the H1 estimator (default: `no`)


### TODO
The acquisition algorithm works as is, but further post processing methods are needed to get results other than a transfer function.

* reverberation time (offline post processing)
* reduction index (on- and off-line)
* FIX: live preview displays up to 2 input channels. More channels should be able to be displayed meaningfully.
* move CLI input arguments to a settings file
* fix plotly compatibility with LaTeX interpreter 
