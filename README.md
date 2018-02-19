# nidaqmxAio 

Input/output acquisition automation using the [nidaqmx](https://github.com/ni/nidaqmx-python) python API and analog input/output modules.

### Prerequisites

You need to be able to run python scripts on your system with version >3.6 installed. nidaqmx supports only Windows OS.
Python prerequisites: numpy, matplotlib, collections, colorednoise, time, nidaqmx.

### Usage

Run `python run_meas.py <simulation time> <sample rate>` with additional arguments for channel input/outputs. 

The script is built without external triggering/clocking devices in mind. Due to system delays that cannot be deterministically calculated, the test signals are zero padded to avoid truncated ends. Type `python run_meas.py --help` for information on the additional arguments. 

#### Input arguments

If you need to run the script with other than the default values, specify the corresponding arguments explicitly.

`-tp <string>` or `--signalType <string>`:
The type of the output signal to be generated. Currently, `white_noise` and `pink_noise` are supported (default: `pink_noise`)

`-pad N` or `--pad_samples N`:
Pad N samples at the end of the signal to avoid input truncation due to system delays (default `5000`)

`-ai X Y ...` or `--channelsIn X Y ...`:
Set the number of analog input channels per module. Begin from the lowest numbered slot and continue. The number of input channels per module are chosen serially. For example, if the argument is set as `-ai 2 4 0`, the first 2 input channels are used from the module connected at the lowest numbered slot on the chassis, the 4 from the second and none from the third. The number of channels and modules set here should match your physical setup (default: `0`)


`-ao X Y ...` or `--channelsOut X Y ...`:
Set the number of analog output channels per module. Begin from the lowest numbered slot and continue. The number of output channels per module are chosen serially. For example, if the argument is set as `-ao 2 4 0`, the first 2 output channels are used from the module connected at the lowest numbered slot on the chassis, the 4 from the second and none from the third. The number of channels and modules set here should match your physical setup. Note that the test signal is copied to each analog output channel. This is due to change and be configurable by the user in future updates (default: `0`)

`-airms N` or `--aiRange N`:
The analog input expected peak to peak range in rms volts (default: `3`)

`-aorms N` or `--aoRange N`:
The analog output expected peak to peak range in rms volts (default: `3`)

`-sv True/False` or `--save_file True/False`:
Choose to save the collected data. You will be asked to provide the desired filename after the measurement (default: `False`)

`-bf N` or `--bufferSize N`:
Buffer size in samples of the input channels (default: `2048`)

### TODO
The acquisition algorithm works as is, but further post processing methods are needed to get meaningful results.

* Include a basic live post-processing procedure for preview and feedback
* Add pros-processing algorithms for transfer functions and impulse responses
