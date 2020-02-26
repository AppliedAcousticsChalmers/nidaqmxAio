# Python libraries
import numpy as np
from scipy import polyfit, polyval
from tqdm import tqdm
import matplotlib.pyplot as plt

# Program libraries
import filtersAndMathtools as flm


def T60(IR, sr, tVec, f_min=20, f_max=10e3, bandWidth="third"):

    """
    Calculates the RT values in the bandwidth specified, using Schroeder's
    backward integration method. Evaluates RT at T_20, T_30 and T_60. the beginning
    of the background noise level needs to be set by clicking in an interactive plot.
    The curvature C Is also estimated.
    """
    db = lambda x: flm.NormAmp2dB(x)


    print("Converting IR into 3rd octave bands...", end="")
    IR_3rd, fvec = flm.band_filter(IR, sr, f_min, f_max, 3, bandWidth)
    IR_3rd = np.array(IR_3rd)
    print("completed")

    T_60 = []
    T_30 = []
    T_20 = []
    IRs = []
    IR2s = []
    Rev_int_curves = []
    line_fits_T20 = []
    line_fits_T30 = []
    line_fits_T60 = []

    print("Calculating T60 for each band")
    pbar1 = tqdm(total=len(IR_3rd-1))
    for idx in range(0,len(IR_3rd-1)):

        IR_current = IR_3rd.real[idx,:]

        #Selecting the start of the noise floor
        mng = plt.get_current_fig_manager()
        mng.window.state('zoomed')
        plt.title(str(fvec[idx]) + " Hz " + str(bandWidth) + "-octave band")
        plt.figtext(.5,.85,"Please select the beginning of the nose floor.", ha="center")
        plt.plot(tVec, db(IR_current**2))
        plt.xlabel('time in sec')
        plt.ylabel('Amplitude in dB')
        plt.ylim(-120, 10)
        x = ""
        while not x:
            x = plt.ginput(1, timeout = -1)
        tau = flm.find_nearest(tVec,x[0][0])[0]
        plt.show(block = False)
        plt.close()

        while True:
            #Calculating the Schroeder curve
            IR2 = IR_current ** 2
            Rev_int = np.zeros((len(IR2[:tau])))
            Rev_int[tau::-1] =(np.cumsum(IR2[tau:0:-1]) / np.sum(IR2[0:tau]))
            Rev_int_dB = db(Rev_int)

            #Fitting a line in the Schroeder curve
            fit_start = flm.find_nearest(Rev_int_dB, -5)[0]
            fit_T20 = flm.find_nearest(Rev_int_dB, -25)[0]
            fit_T30 = flm.find_nearest(Rev_int_dB, -35)[0]
            fit_T60 = flm.find_nearest(Rev_int_dB, -65)[0]

            #T20
            x_T20 = tVec[fit_start:fit_T20]
            Rev_int_fit_T20 = Rev_int_dB[fit_start:fit_T20]
            (b_T20, a_T20) = polyfit(x_T20, Rev_int_fit_T20, 1)
            xr_T20=polyval([b_T20,a_T20],tVec)

            #T30
            x_T30 = tVec[fit_start:fit_T30]
            Rev_int_fit_T30 = Rev_int_dB[fit_start:fit_T30]
            (b_T30, a_T30) = polyfit(x_T30, Rev_int_fit_T30, 1)
            xr_T30=polyval([b_T30,a_T30],tVec)

            #T60
            x_T60 = tVec[fit_start:fit_T60]
            Rev_int_fit_T60 = Rev_int_dB[fit_start:fit_T60]
            (b_T60, a_T60) = polyfit(x_T60, Rev_int_fit_T60, 1)
            xr_T60=polyval([b_T60,a_T60],tVec)

            #Fitting the guiding lines
            x_5dB = polyval([0,max(Rev_int_dB)-5],tVec)
            x_25dB = polyval([0,max(Rev_int_dB)-25],tVec)
            x_35dB = polyval([0,max(Rev_int_dB)-35],tVec)
            x_65dB = polyval([0,max(Rev_int_dB)-65],tVec)

            #Calculating T_20
            t20_5dB = tVec[np.argmin(abs(xr_T20-x_5dB))]
            t20_25dB = tVec[np.argmin(abs(xr_T20-x_25dB))]
            T20 = -60 / b_T20

            #Calculating T_30
            t30_5dB = tVec[np.argmin(abs(xr_T30-x_5dB))]
            t30_35dB = tVec[np.argmin(abs(xr_T30-x_35dB))]
            T30 = -60 / b_T30

            #Calculating T_60
            t60_5dB = tVec[np.argmin(abs(xr_T60-x_5dB))]
            t60_65dB = tVec[np.argmin(abs(xr_T60-x_65dB))]
            T60 = -60 / b_T60

            #Calculation of nonlinearity
            xi_20 = degree_of_nonlinearity(Rev_int_dB, xr_T20[:tau])
            xi_30 = degree_of_nonlinearity(Rev_int_dB, xr_T30[:tau])
            xi_60 = degree_of_nonlinearity(Rev_int_dB, xr_T60[:tau])

            #Calculation of curvature
            C = degree_of_curvature(T20, T30)


            #Plotting the results
            mng = plt.get_current_fig_manager()
            mng.window.state('zoomed')
            plt.plot(tVec, db(IR2), label=r"$IR^2$", linewidth=1.2, color='blue' )
            plt.plot(tVec[:tau], Rev_int_dB, label='Schroeder curve, C=' + str(round(C, 2)) + "%", linewidth=1.2, color='red' )
            plt.plot(tVec, x_5dB, label='-5 dB', linestyle='--', linewidth=0.5, color='black' )
            plt.plot(tVec, x_25dB, label='-25 dB', linestyle='--', linewidth=0.5, color='black')
            plt.plot(tVec, x_35dB, label='-35 dB', linestyle='--', linewidth=0.5, color='black')
            plt.plot(tVec, x_65dB, label='-65 dB', linestyle='--', linewidth=0.5, color='black')
            plt.plot(tVec, xr_T20, label='T20:' + str(round(T20, 2)) + r"s,$\xi_{20}$=" + str(round(xi_20, 2)) + r"$\perthousand$" , linestyle=':', linewidth=1.2, color='dimgray')
            plt.plot(tVec, xr_T30, label='T30:' + str(round(T30, 2)) + r"s,$\xi_{30}$=" + str(round(xi_30, 2)) + r"$\perthousand$", linestyle='--', linewidth=1.2, color='dimgray')
            plt.plot(tVec, xr_T60, label='T60:' + str(round(T60, 2)) + r"s,$\xi_{60}$=" + str(round(xi_60, 2)) + r"$\perthousand$", linewidth=1.2, color='dimgray')
            try:
                plt.xlim(0, tVec[int(tau + sr)])
            except:
                plt.xlim(0, tVec[-1])
            plt.ylim(-120, 10)
            plt.title(str(fvec[idx]) + " Hz " + str(bandWidth) + "-octave band")
            plt.figtext(.5,.85,r"Click to readjust $\tau$ or press CR to continue", ha="center")
            plt.legend(loc=1)
            plt.xlabel('time in sec')
            plt.ylabel('Amplitude in dB')
            x = plt.ginput(1, timeout = -1)
            if x:tau = flm.find_nearest(tVec,x[0][0])[0]
            tVec_cut = tVec[:tau]
            plt.legend()
            plt.show(block=False)
            plt.close()
            if not x:
                T_60.append(T60)
                T_30.append(T30)
                T_20.append(T20)
                IRs.append(IR_current)
                IR2s.append(IR2)
                Rev_int_curves.append([Rev_int, {"C": C}])
                line_fits_T20.append([[xr_T20, x_5dB, x_25dB], {'xi': xi_20}])
                line_fits_T30.append([[xr_T30, x_5dB, x_35dB], {'xi': xi_30}])
                line_fits_T60.append([[xr_T60, x_5dB, x_65dB], {'xi': xi_60}])
                break

        pbar1.update()
    pbar1.close()

    return T_20, T_30, T_60, tau, fvec, IRs, IR2s, Rev_int_curves, line_fits_T20, line_fits_T30, line_fits_T60

def degree_of_nonlinearity(actual_level, estimated_level):
    '''
    Degree of nonlinearity according to Annex B2 of ISO 3382-2:2008
    Inputs:
    - actual level: level in db of the decay curve.
    - estimated level: level in db of the estimated linear regression curve for the same sample
    Output:
    - xi: the non-linearity parameter, in thousand percentile.
    '''
    L_avg = np.sum(actual_level) / len(actual_level)
    r2 = (np.sum((estimated_level - L_avg) ** 2)) / (np.sum((actual_level - L_avg) ** 2))
    xi = 1000 * (1 - r2)
    return xi

def degree_of_curvature(T_20, T_30):
    '''
    Degree of curvature according to Annex B3 of ISO 3382-2:2008
    '''
    return 100 * ((T_30 / T_20) - 1)
