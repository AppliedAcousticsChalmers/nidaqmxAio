# nidaqAio 

Input/output acquisition automation using the [nidaqmx](https://github.com/ni/nidaqmx-python) python API and analog input/output modules.

### Prerequisites

You need to be able to run python scripts on your system with version >3.6 installed. nidaqmx supports only Windows OS.
Python prerequisites: numpy, matplotlib, collections, colorednoise, nidaqmx.

### Usage

Run with `python nidaqAio.py <simulation time> <sample rate>` to perform the transfer function measurements. This is built without external triggering/clocking devices in mind. Due to system delays that cannot be deterministically calculated, the test signals are zero padded to avoid truncated ends. Type `python nidaqAio.py --help` for more information for more script inputs. Inputs and outputs are used serially. To define them correctly, define the ins/outs for each module in a list, starting from the
first on the chassis module. For example, if 2 inputs are needed from the module in slot 2 and 1 output from the module in slot 4, run the script as `python nidaqAio.py <simulation time in s> <sample rate> -ai 2 -ao 1`. In the case more modules are connected but not needed, include them using 0 ins/outs. For example if an additional input module is connected at slot 5, run  `python nidaqAio.py <simulation time in s> <sample rate> -ai 2 0 -ao 1`.

### TODO
The acquisition algorithm works as is, but further post processing methods are needed to get meaningful results.

* Include a basic live post-processing procedure for preview and feedback
* Add pros-processing algorithms for transfer functions and impulse responses
