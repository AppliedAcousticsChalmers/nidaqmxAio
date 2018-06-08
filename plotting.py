import numpy as np
import plotly.offline as py
import plotly.graph_objs as go
from plotly import tools
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

import utils as gz


def singlePlot(plots):
    """
    Outputs a signle plot with arbitary many lines of data.
    The data should be given in a dictionary formated as:

    plots = {'0':[tVec,IR.real[1,:],'Impulse Responce'],
              'xAxisTitle': 'frequency in Hz',
              'yAxisTitle': 'H',
              'plotTitle': 'Plot title'}

    The each line of data should be registered to a key,
    named '0', '1', etc. The data assigned to the key should
    be a list in the form:
    [x axis data, y axis data, legend name]

    The last three keys refer to the tittles of the x,y axis
    and the title of the plot respectivelly.
    """

    trace = []
    for idx in range(0, len(plots.keys())-4):
        trace.append(go.Scatter(
            x = plots[str(idx)][0],
            y = plots[str(idx)][1],
            name = plots[str(idx)][2]
        ))
    layout = go.Layout(
        title=plots['plotTitle'],
        xaxis=dict(
        title=plots['xAxisTitle'],
        type=plots['scale'],
        autorange=True
        ),
        yaxis=dict(
        title=plots['yAxisTitle'],
        )
    )
    fig = go.Figure(data=trace, layout=layout)
    py.plot(fig, filename=plots['plotTitle'] + '.html')

    return

def irSubPlot(plots, filename, title):
    """
    Outputs a plot with 4 subplots for an IR measurement conducted
    using a Noise or a Sweep signal.

    The plots pressented are:

    Noise signal:
                 [1,1] - Transfer function,
                 [1,2] - Coherence,
                 [2,1] - Impulse responce,
                 [2,2] - Phase

    Sweep signal:
                 [1,1] - Transfer function,
                 [1,2] - Spectrogram,
                 [2,1] - Impulse responce,
                 [2,2] - Phase

    The type of the [1,2] subplot is determined by the name of the key.
    If 'gamma2' is used, a list with [frequency vector, coherence] is expected.
    If 'spectrogram' is used, a list with [frequency vector, time vector, magnitude] is expected.

    The input is a ditionary in the form:

    plots = {'TF':[frequency vector, Transfer function],
             'gamma2':[frequency vector, coherence], or
             'spectrogram':[frequency vector, time vector, magnitude],
             'IR':[time vectror, Impulse responce],
             'phase':[frequency vectror, Phase]}

    filename - the name of the output file
    title - The title of the plot
    """

    trace1 = go.Scatter(
            x=plots['TF'][0],
            y=plots['TF'][1]
        )

    if ("gamma2" in plots.keys()):
        trace2 = go.Scatter(
            x=plots['gamma2'][0],
            y=plots['gamma2'][1]
           )
    elif ("spectrogram" in plots.keys()):
        trace2 = go.Heatmap(
                x=plots['spectrogram'][1],
                y=plots['spectrogram'][0],
                z=plots['spectrogram'][2],
                colorbar=dict(y=0.815, len=0.41)
            )

    trace3 = go.Scatter(
            x=plots['IR'][0],
            y=plots['IR'][1]
        )

    trace4 = go.Scatter(
            x=plots['phase'][0],
            y=plots['phase'][1],
        )

    if ("gamma2" in plots.keys()):
        fig = tools.make_subplots(rows=2, cols=2, subplot_titles=('HdB', 'gamma2', 'IR', 'H_phase'))
    elif ("spectrogram" in plots.keys()):
        fig = tools.make_subplots(rows=2, cols=2, subplot_titles=('HdB', 'Spectrogram', 'IR', 'H_phase'))

    fig.append_trace(trace1, 1, 1)
    fig.append_trace(trace2, 1, 2)
    fig.append_trace(trace3, 2, 1)
    fig.append_trace(trace4, 2, 2)

    if ("gamma2" in plots.keys()):
        fig['layout'].update(legend=dict(y=0.5),
                             title=title,
                             xaxis1=dict(title='frequency in Hz', type='log', autorange=True),
                             xaxis2=dict(title='frequency in Hz', type='log', autorange=True),
                             xaxis3=dict(title='time in s', type='lin', autorange=True),
                             xaxis4=dict(title='frequency in Hz', type='log', autorange=True)
        )
    elif ('spectrogram' in plots.keys()):
        fig['layout'].update(legend=dict(y=0.5),
                         title=title,
                         xaxis1=dict(title='frequency in Hz', type='log', autorange=True),
                         xaxis2=dict(title='time in s', type='lin', autorange=True),
                         yaxis2=dict(title='frequency in Hz', type='log', autorange=True),
                         xaxis3=dict(title='time in s', type='lin', autorange=True),
                         xaxis4=dict(title='frequency in Hz', type='log', autorange=True)
    )
    py.plot(fig, filename=filename + '.html')


class livePlot(object):
    '''
    Class that handles the preallocation, update, and clipping of the live plotting of nidaqio.

    Methods:
    __init__: Initializes the plotting object.
    livePlotClipp: Detects and displays the clipping of the time signals.
    livePlot_plotPreallocation: creates the figure, prior to data acquisition, according to the system settings.
    livePlot_2CH_Noise: used by livePlot_plotPreallocation if the signal is noise and the number of channels is 2.
    livePlotCreation: used by the livePlot_plotPreallocation in all the other cases.
    livePlotUpdate: updates the values of the plots in all cases except the 2CH_Noise.
    livePlotUpdate_2CH_Noise: updates the values of the plots in the 2CH_Noise case.


    '''

    def __init__(self, args, bufferSize, sample_rate, selection, channelNames):
        '''
        Method initialazes the plotting object.

        Inputs:
        - number_of_channels: the list containing the number of channels used.
        Ex [2,1] for 2 modules using 2 and 1 channels respectivelly.

        - args: The parser arguments input of the CLI.
        - bufferSize: The bufferSize used by the program. This is not inputed by the args because the program
        ensures that bufferSize is a power of 2 after the CLI input.
        - sample_rate: the sample rate used. This is not inputed by the args for the same reason as above.
        - selection: the selected by the user reference channel.
        - channelNames: the list with the channelNames provided by the niDAQmx API.
        '''

        #Initial arguements
        self.args = args
        self.number_of_channels = sum(self.args.channelsIn)
        self.bufferSize = bufferSize
        self.sample_rate = sample_rate
        self.selection = selection
        self.channelNames = channelNames
        self.signalType = args.signalType[0]
        self.signalType = self.signalType.split("_")[0]

        self.chidx = list(np.arange(self.number_of_channels))
        if self.selection != 'none':
            self.channelNames.remove(self.channelNames[self.selection])
            self.chidx.remove(self.chidx[self.selection])


        #Time and frequency vectors
        self.tVec = np.linspace(0, self.bufferSize / self.sample_rate, self.bufferSize)
        self.fftfreq = np.fft.rfftfreq(self.bufferSize, 1 / self.sample_rate)
        self.fftfreq = self.fftfreq[:-1]

        #Initial clipping parameters
        self.pen_colours = [(173, 255, 47, 130), #Reference channel
                            (200, 200, 200, 130), #Regular channel
                            (255, 0, 0, 130), #Clipped block (red)
                            (255, 127, 80, 130) #Previous block clipped (orange)
        ]
        self.pen = []
        self.clipped = []
        for i in range(0, self.number_of_channels):
            if i == self.selection and self.selection != 'none':
                self.pen.append(self.pen_colours[0])
                self.clipped.append(False)
            else:
                self.pen.append(self.pen_colours[1])
                self.clipped.append(False)



        #Window initialazation
        global app
        app = QtGui.QApplication([])
        global winGraph
        winGraph = pg.GraphicsLayoutWidget()
        winGraph.setWindowTitle('And awaaaaay we go!')
        winGraph.resize(1000, 600)
        winGraph.show()
        pg.setConfigOptions(antialias=True)
        self.p = []
        self.curve = []
        self.downsample = self.number_of_channels*0+1


        self.livePlot_plotPreallocation()

        return

    def livePlotClipp(self,data):
        '''
        The method checks if the time signals clipp and changes the pen colours acordingly.

        Red - If the signal clips durring the current block.
        Orange - If the signal has clipped in one of the previous blocks.
        '''
        current_max_val = np.amax(abs(data),1)

        for i in range(0, len(self.pen)):
            if current_max_val[i] >= self.args.aiRange:
                self.pen[i] = (255, 0, 0, 130)
                self.clipped[i] = True
            elif not self.clipped[i]:
                continue
            else:
                self.pen[i] = (255, 127, 80, 130)

        return


    def livePlot_plotPreallocation(self):
        """
        The method creates the graphs according to the input provided by the user in the __init__ part

        The cases are:

            - Measurement without reference: Only the time signals are displayed in one subplot each
        (for obvious reasons this is always the case for single channel measurements)

            - Measurement wth reference: Each subplot pressents the time signals for one channel
        and the reference channel for easier comparison.

            - 2 Channel with reference using noise signal: The time signals are displayed as above, as well as
        the instantaneus spectra, the H1 transfer function, the IR and coherence in different subplots.

        """

        #Single channel measurement
        if self.number_of_channels == 1:
            self.livePlot_plotCreation(str("Time Signals: " + self.channelNames[0]), xAxisMode="time",
                                  yRange=[-self.args.aiRange, self.args.aiRange], colspan=1, nPlots=1, plt_idx=0)
            return

        #With reference channel
        if self.selection != 'none':
            #2CH Noise measurements
            if self.number_of_channels == 2 and "noise" in self.signalType:
                self.livePlot_2CH_Noise(chName=self.channelNames[0])
            #All other measurements using reference
            else:
                for i in range(0, self.number_of_channels - 1):
                    self.livePlot_plotCreation(title='Time Signals: ' + self.channelNames[i], xAxisMode="time",
                                               yRange=[-self.args.aiRange, self.args.aiRange], colspan=1, nPlots=2, plt_idx=i)
                    winGraph.nextRow()
        #No reference channel
        else:
            for i in range(0, self.number_of_channels):
                self.livePlot_plotCreation(title='Time Signals: ' + self.channelNames[i], xAxisMode="time",
                                           yRange=[-self.args.aiRange, self.args.aiRange], colspan=1, nPlots=1, plt_idx=i)

                winGraph.nextRow()

        return



    def livePlot_2CH_Noise(self, chName=""):
        '''
        This method preallocates the plot for the TF measurement, using noise signal,
        2 channels one of which is reference

        Input: chName: The input channel name (not the reference channel).

        '''
        #Time signals
        self.livePlot_plotCreation(title="Time Signals: " + chName,
                              xAxisMode='time',
                              yRange=[-self.args.aiRange, self.args.aiRange],
                              colspan=2,
                              nPlots=2,
                              plt_idx=0)

        winGraph.nextRow()

        #Spectrum
        self.livePlot_plotCreation(title="Spectrum: " + chName,
                             xAxisMode="freq",
                             yRange=[],
                             colspan=1,
                             nPlots=2,
                             plt_idx=1)

        #Transfer function
        self.livePlot_plotCreation(title="Transfer function: " + chName,
                             xAxisMode="freq",
                             yRange=[],
                             colspan=1,
                             nPlots=1,
                             plt_idx=2)

        winGraph.nextRow()

        #Impulse Responce
        self.livePlot_plotCreation(title="Impulse Responce: " + chName,
                             xAxisMode="time",
                             yRange=[],
                             colspan=1,
                             nPlots=1,
                             plt_idx=3)


        #Coherence
        self.livePlot_plotCreation(title="Coherence: " + chName,
                             xAxisMode="freq",
                             yRange=[0, 1.1],
                             colspan=1,
                             nPlots=1,
                             plt_idx=4)


        return

    def livePlot_plotCreation(self, title, xAxisMode="time", yRange=[], colspan=1, nPlots=1, plt_idx=0):
        '''
        Creates a plot with a set number of curves in it.

        Inputs:
        - title: The plot title
        - xAxisMode: sets the x axis to time or frequency and the scale to lin or log respectivelly.
        The values of the axis are calculated by __init__.
        - yRange: sets the limits for the y axis
        - colspan: sets the number of colums the subplot will occupy.
        - nPlots: the number of curves created in the same plot.
        - plt_idx: indexing for the current plot. Used in the case many subplots with different amount of curves are needed.
        Example: If a figure with 2 subplots containg 2 cureves each is needed, the method should be called 2 times.
        First using nPlots=2 and plt_idx=0 and the second time with nPlots=2 and plt_idx=1
        '''

        #Setting up the x Axis
        if xAxisMode == "time":
            xVec = self.tVec
            logMode = False
        elif xAxisMode == "freq":
            xVec = self.fftfreq
            logMode = True

        self.p.append(winGraph.addPlot(title=title, colspan=colspan))
        self.p[plt_idx].showGrid(True, True)
        self.p[plt_idx].setLogMode(logMode, False)
        if yRange: self.p[plt_idx].setRange(yRange=[yRange[0], yRange[1]])
        for i in range(0,nPlots):
            self.curve.append(self.p[plt_idx].plot(xVec, np.zeros(len(xVec))))


        return

    def livePlotUpdate(self, data):
        '''
        The method updates the data for the plots, displaying the time signals, when the number of channels
        is other than 2 and/or the signal used is not noise.

        Input: The vector of the current block of data.
        '''
        if self.selection == 'none':
            for i in range(0, self.number_of_channels):
                self.curve[i].setData(self.tVec, data[i,:], pen=self.pen[i], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
        else:
            ii = 0
            for i in range(0, 2 * (self.number_of_channels - 1) - 1, 2):
                self.curve[i].setData(self.tVec, data[self.selection,:], pen=self.pen[self.selection], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
                self.curve[i+1].setData(self.tVec, data[self.chidx[ii],:], pen=self.pen[self.chidx[ii]], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
                ii += 1

        pg.QtGui.QApplication.processEvents()

        return

    def livePlotUpdate_H1(self, current, spectra, HdB, IR, gamma2):
        '''
        The method updates the plots created for the TF estimator, using Noise signals and a reference channel.

        Inputs:
        - current: The vector with the current block of time signals.
        - spectra: The instantaneus spectra of the current block.
        - HdB: The estimated TF in dB.
        - IR: The corresponding IR.
        - gamma2: The coherence of the measurement.
        '''

        self.curve[0].setData(self.tVec, current[self.selection,:], pen=self.pen[self.selection], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
        self.curve[1].setData(self.tVec, current[self.chidx[0],:], pen=self.pen[self.chidx[0]], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
        self.curve[2].setData(self.fftfreq, gz.amp2db(spectra[self.selection, 0:int(self.bufferSize//2)]), pen=self.pen_colours[0], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
        self.curve[3].setData(self.fftfreq, gz.amp2db(spectra[self.chidx[0], 0:int(self.bufferSize//2)]), pen=self.pen_colours[1], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
        self.curve[4].setData(self.fftfreq, HdB[self.chidx[0], ...], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
        self.curve[5].setData(self.tVec, IR.real[self.chidx[0], ...], antialias=True, downsample=self.downsample, downsampleMethod='subsample')
        self.curve[6].setData(self.fftfreq, gamma2[self.chidx[0], ...], antialias=True, downsample=self.downsample, downsampleMethod='subsample')

        pg.QtGui.QApplication.processEvents()

        return

