# nidaqmxAio 

Input/output acquisition automation using the [nidaqmx](https://github.com/ni/nidaqmx-python) python API and analog input/output modules.

### Prerequisites

You need to be able to run python scripts on your system with version >3.6 installed. nidaqmx supports only Windows OS.
Python prerequisites found in `requirements.txt`.

### Usage
Run `python run_meas.py -c <configuration_file_name.yml` to run using the specified configuration file. An example configuration file is provided.

or

Run `python run_meas.py <simulation time> <sample rate>` `<new meassurement>` with additional arguments for channel input/outputs. If `--newMeasurement` is set to `1` a new measurement is going to take place (default `0`).

The script is built without external triggering/clocking devices in mind and therefore only specific sample rates can be used(refer to the NI cDAQ modules manuals for more information). Due to system delays that cannot be deterministically calculated, the test signals are zero padded to avoid truncated ends. Type `python run_meas.py --help` for information on the additional arguments. Measurements are saved in the directory `acquired_data\measurements_<YYMMDD>\` in the format `meas_<YYMMDD>.py` by specifying the desired filename, using `-sv filename`. During measurements, a live post-processing preview is displayed, presenting the spectra of reference and measured channel, transfer function (using the H1 estimator), impulse response and coherence. In case the H1 estimator calculation is not possible, the corresponding time signals are displayed. The reference channel can be selected from a list after running the script (leave blank to select the first channel on the list). The option for no reference channel is also available. 

#### Input arguments

If you need to run the script with other than the default values, specify the corresponding arguments explicitly.

`-sig <string> <float> <float>` or `--signalType <string> <float> <float>`:
The type of the output signal to be generated. Currently, `[noise_white]`, `[noise_pink]`, `[sweep_linear, f0, f1]`, `[sweep_logarithmic, f0, f1]`, `[tone, f0]` and `[matLab, filename]`  are supported (default: `pink_noise`). The values `f0`, `f1` correspond to the starting and stopping frequency of the sweep respectively. In the case of the tone signal `f0` corresponds to the frequency of the signal. The `matLab` reads a `.mat` file and outputs it. The file should contain a structure with a variable named `audio` in it. 

`-pad <int>` or `--pad_samples <int>`:
Pad N samples at the end of the signal to avoid input truncation due to system delays (default `5000`)

`-ai <int> <int> ...` or `--channelsIn <int> <int> ...`:
Set the number of analog input channels per module. Begin from the lowest numbered slot and continue. The number of input channels per module are chosen serially. For example, if the argument is set as `-ai 2 4 0`, the first 2 input channels are used from the module connected at the lowest numbered slot on the chassis, 4 from the second and none from the third. The number of channels and modules set here should match your physical setup (default: `0`)

`-ao <int> <int> ...` or `--channelsOut <int> <int> ...`:
Set the number of analog output channels per module. Begin from the lowest numbered slot and continue. The number of output channels per module are chosen serially. For example, if the argument is set as `-ao 2 4 0`, the first 2 output channels are used from the module connected at the lowest numbered slot on the chassis, 4 from the second and none from the third. The number of channels and modules set here should match your physical setup. Note that the test signal is copied to each analog output channel. This is due to change and be configurable by the user in future updates (default: `0`)

`-aip <float>` or `--aiPeak <float>`:
The analog input absolute peak expected in volts (default: `5.0`)

`-aop <float>` or `--aoPeak <float>`:
The analog output absolute peak expected in volts (default: `0.001`)

`-sv <string>` or `--save_file <string>`:
 Name of the saved output file.  (default: `measurement`)

`-bf <int>` or `--bufferSize <int>`:
Buffer size in samples of the input channels. This is used also as the time window for the analysis done at both the live post-processing preview and the Transfer function post-processing script. In the deconvolution analysis, corresponds to the block size of the spectrogram and the window that the IR is expected to end in samples (default: `8192`).

`-cal <string>` or `--calibration <string>`:
Specifies the filename of the calibration data file. If it's not present ask if a new calibration measurement needs to be performed. If the case a new measurement is needed, the calibration measurement script is run and the calibration coefficient and microphone sensitivity are saved. If not defined, the sensitivity used is for the microphone B&K type 4190.

`-micA <float>` or `--micAmp <float>`:
 Specifies any external microphone amplification used during the measurement. (default: `1`)
 
`-sens <float>` or `--sensitivity <float>`:
 Specifies the sensitivity of the used microphone default `47.1` mV/Pa.

`-pp <string>` or `--postProcess <string>`:
Run the chosen post processing script. You will be asked to select the directory and a the file from a list (typing `all` will cause the script to run for all the files in the selected directory with the file extension`.np[yz]`, except the calibration files). Results will be saved as `<measurement_name>_TFs`. Current choices:
`TF` - Estimates the transfer function using the H1 estimator, or the deconvolution method. The correct method is chosen automatically depending on the signal type used in the measurement.
`RAC` - Estimates the reverberation time using the back-wards integration method. This script always causes the `TF` post-porssesing script to run in order to calculate the impulse response using the desired calibration data and block size.

`-cT <int>` or `--cutoffTime <int>`:
 Causes the algorithm to keep measuring for the specified amount of seconds after the output has stopped. (default: `0`)
 
`-plt <string> <string> ...` or `--plotting <string> <string> ...`:
  Causes the algorithm to output the specified plots, if the corresponding script is run. Current choices:
  `live` - Plots live, the time signals or processed data during the data acquisition.
  `TF` - Plots the transfer function, coherence/spectrogram (depending on the method), impulse response and phase, when the `TF` post-prossesing script is run.
  `timeSig` - Plots the raw time signals, when the measurement script is run.
  `T60_one_band` - Plots each step of the back-wards integration method for each octave band calculated, when the `RAC` script is run.
  `T60_3rd` - Plots the T60 in 3rd octave bands, when the `RAC` script is run.
  
`-fRange <float> <float> <string>` or `--frequencyRange <float> <float> <string>`:
 Sets the frequency range for the `RAC` script. The format is a list as - [lowest frequency band, highest frequency band, bandwidth]. Bandwidth can be set to `'third'` or `'one'` for third and one octave bands respectively (default: `[20, 10000, 'third']`).

`-refCH <float>` or `--refferenceChannel <float>`: 
 Presets the channel that is going to be used as reference, bypassing the command line dialog, when the measurement script is run (default=`'"""'`).
 
 `-mt <string>` or `--note <string>`:
 Saves a text note on the output file (default=`""`). 
