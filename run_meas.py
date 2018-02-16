import argparse
import numpy as np
import ni_tools as ni

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("time", help="Measurement time in s", type=float)
    parser.add_argument("sampleRate", help="Measurement sample rate in Hz", type=float)

    # Optional arguments
    parser.add_argument("-tp", "--signalType", help="Type of test signal. Options: \"pink_noise\", \"white_noise\"", type=str, choices={"white_noise", "pink_noise", "debug"} )
    parser.add_argument("-pad", "--pad_samples", help="Number of zeros added to the end of the signal.", type=int)
    parser.add_argument("-ai", "--channelsIn", help="Number of channels per input module", nargs='+', type=int)
    parser.add_argument("-ao", "--channelsOut", help="Number of channels per output module", nargs='+', type=int)
    parser.add_argument("-airms", "--aiRange", help="Analog input in Vrms", type=float)
    parser.add_argument("-aorms", "--aoRange", help="Analog output in Vrms", type=float)
    parser.add_argument("-sv", "--save_file", help="Option for exporting data to file.", type=bool)
    parser.add_argument("-bf", "--bufferSize", help="Input buffer size in samples", type=int)
    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 0.0')

    # Parse arguments
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parseArguments()
    if not args.signalType: args.signalType = "pink_noise"
    if not args.pad_samples: args.pad_samples = 5000
    if not args.channelsIn: args.channelsIn = [0]
    if not args.channelsOut: args.channelsOut = [0]
    if not args.aiRange: args.aiRange = 3
    if not args.aoRange: args.aoRange = 3
    if not args.save_file: args.save_file = False
    if not args.bufferSize: args.bufferSize = 2048

    print("You are running the script with arguments: ")
    for a in args.__dict__:
        print(str(a) + ": " + str(args.__dict__[a]))
    meas = ni.ni_io_tf(args)
    if args.save_file==True:
        filename = input("Please enter the desired name, for the output file:")
        np.save(filename, meas)


