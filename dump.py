import os

__author__ = 'adel'
import numpy as np
from lmfit import minimize, Parameters, Parameter, report_fit
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter
from plotter_D_A import plot_multiple, plot_unique, plot_double

"Manip post-processing"
"""
TODO:
ISOTROPIC TREATMENT:
    Comparaison between Volume measured by image treatment and the spring position ==> Rho*g*V-kx = 0
    fit V(P) with the theoretical formula (Church) and deduce the E (3D elastic modulus).
    Verify the difference between Isotropic forward and Isotropic backward.
BUCKLING:
    fit Volume as: V(t) = V0 * exp(-tau_v*t)*cos(w_v*t+shift)+d
    derivate x over time to get x',x"
    fit m*x" = archimed(fitted)+k*x -alpha*x' and deduce alpha. Compare it to alpha sphere.
    integrate x over time
    deduce thrust (supposing alpha constant)
    """
# 5


K_THRESHOLD_MAX = 10
K_THRESHOLD_MIN = 6
DELTA_T = 0.1
COUNTER = 0


def iso_writer(rootfolder, list_fit):
    wfile = open(rootfolder + "isotropic_fit_output" + '.txt', "w")
    wfile.writelines("Experiment_name\tSlope\tOffset\t")
    for el in list_fit:
        wfile.writelines("\n")
        wfile.writelines(el[0] + '\t')
        wfile.writelines("{0:.5e}".format(el[1]) + '\t')
        wfile.writelines("{0:.5e}".format(el[2]) + '\t')
    wfile.close()


def parameters_writer(rootfolder, list_fit):
    wfile = open(rootfolder, "w")
    for el in list_fit[0].keys():
        wfile.writelines(el + "\t")
    wfile.writelines(el + "\n")
    for d in list_fit:
        for val in d.keys():
            wfile.writelines("{0:.3e}".format(d[val]) + "\t")
        wfile.writelines("\n")
    wfile.close()


def column(matrix, i):
    return [row[i] for row in matrix]


def read_results_file(filename, timeCut=None):
    print filename
    try:
        results = np.genfromtxt(filename, dtype=None, delimiter='\t')
    except ValueError:
        results = np.genfromtxt(filename, dtype=None, delimiter='\t', usecols=np.arange(0, 10))
    results_keys = list(results[0])
    results = results[1:]
    dict_results = {}
    index_cut = 0
    try:
        dict_results[results_keys[0]] = np.array(column(results, 0), dtype=np.float)
    except ValueError:
        dict_results[results_keys[0]] = np.array(column(results, 0), dtype=np.string_)
    if timeCut is not None:
        index_cut = np.amin(np.argwhere(dict_results[results_keys[0]] > timeCut))
        dict_results[results_keys[0]] = dict_results[results_keys[0]][index_cut:] - dict_results[results_keys[0]][
            index_cut]
    for i in xrange(len(results_keys)):
        try:
            dict_results[results_keys[i]] = np.array(column(results, i), dtype=np.float)[index_cut:]
        except ValueError:
            dict_results[results_keys[i]] = np.array(column(results, i), dtype=np.string_)[index_cut:]
    return dict_results


def isotropic(filename, v0, x0, rho, g=9.81):
    try:
        results = read_results_file(filename)
        Volume_array_normalized = v0 - results['External_volume']
        spring_position_array_normalized = x0 - results['Ylow']
        params = Parameters()
        params.add('A', value=1)
        params.add('B', value=0)
        try:
            result = minimize(residual_isotropic, params, args=(spring_position_array_normalized, Volume_array_normalized))
        except TypeError:
            result = minimize(residual_isotropic, params, args=(spring_position_array_normalized, Volume_array_normalized),method="nelder")
        v = result.params.valuesdict()
        x_th = np.arange(np.amin(Volume_array_normalized), np.amax(Volume_array_normalized), 0.0000001)
        y_th = v['A'] * x_th + v['B']
        # report_fit(result.params, min_correl=0.5)
        filename_list = filename.split("\\")
        txt = filename_list[-1].split("_")
        txt = "_".join(txt[0:-1])
        k = -rho * g / v['A']
        print "Stiffness: " + str(k)
        list_results = [txt, k, v['B']]
        return list_results
    except:
        return None


def residual_isotropic(parameters, deltaX, deltaV):
    v = parameters.valuesdict()
    deltaX_fit = v['A'] * deltaV + v['B']

    return deltaX_fit - deltaX


def volume_fit(time, volumeexp,far_field):
    params = Parameters()
    params.add('amp', value=-1e-6)
    params.add('tau_short', value=0.03, min=0.001)
    params.add('omega', value=400, min=1, max=600)
    # params.add('omega',   expr='omega0*sqrt(1-(1/(tau_short*omega0)**2))')
    params.add('shift', value=1, min=0, max=2 * np.pi)
    params.add('V0', value=1e-5, min=0.0000001)
    params.add('Vinf', value=5e-5, min=0)
    params.add('tau_long', value=0.1, min=0.0001)

    result = minimize(residual_volume, params, args=(time, volumeexp,far_field), method='leastsq', maxfev=int(5e4), xtol=1e-8,
                      ftol=1e-8)
    v = result.params.valuesdict()
    t = time
    V_th = residual_volume(params, t, volumeexp,far_field, fit=False)
    # plt.plot(time - t0, volumeexp, 'ro', label="Experimental")
    # plt.plot(t_th - t0, V_th, 'b-', label="Theoretical")
    # plt.title("Volume")
    # plt.legend()
    # plt.show()
    report_fit(result.params, min_correl=0.5)
    return v, V_th


def residual_volume(parameters, t, deltav,far_field =True, fit=True,):
    v = parameters.valuesdict()
    t0 = t[0]
    if not far_field:
        deltav_theoretical = v['amp'] * np.exp(-((t - t0) / v['tau_short'])) * np.sin(v['omega'] * (t - t0) + v['shift']) + \
                         v['Vinf'] * (1 - (1 - v['V0'] / v['Vinf']) * np.exp(-(t - t0) / v['tau_long']))
    else:
        deltav_theoretical =v['Vinf'] * (1 - (1 - v['V0'] / v['Vinf']) * np.exp(-(t - t0) / v['tau_long']))

    if fit:
        return deltav_theoretical - deltav
    else:
        return deltav_theoretical


def r_fit(time, R, volume_params,far_field):
    t_short = volume_params['tau_short']
    omega = volume_params['omega']
    t_long = volume_params['tau_long']
    params = Parameters()
    params.add('amp', value=-0.02)
    params.add('tau_short', value=t_short, vary=True)
    params.add('omega', value=omega, vary=True)
    # params.add('omega',   expr='omega0*sqrt(1-(1/(tau_short*omega0)**2))')
    params.add('shift', value=1, min=0, max=2 * np.pi)
    params.add('Rinf', value=-0.01, )
    params.add('R0', value=-0.05, )
    params.add('tau_long', value=t_long, vary=True)

    result = minimize(residual_r, params, args=(time, R,far_field), method='leastsq', maxfev=int(5e4), xtol=1e-8,
                      ftol=1e-8)
    report_fit(result.params, min_correl=0.5)
    v = result.params.valuesdict()
    t = time
    R_th = residual_r(params, time, R, fit=False)

    return v, R_th


def residual_r(parameters, t, r,far_field=True, fit=True):
    v = parameters.valuesdict()
    t0 = t[0]
    if not far_field:
        r_theoretical = v['amp'] * np.exp(-((t - t0) / v['tau_short'])) * np.sin(v['omega'] * (t - t0) + v['shift']) + v['Rinf'] * (1 - (1 -v['R0'] /v['Rinf']) * np.exp(-(t - t0) /v['tau_long']))
    else:
        r_theoretical =  v['Rinf'] * (1 - (1 -v['R0'] /v['Rinf']) * np.exp(-(t - t0) /v['tau_long']))
    if fit:
        return r_theoretical - r
    else:
        return r_theoretical


def filterer(yg, win_size, pol_order, delta, double_filter=False, show_flag=False):
    if not double_filter:
        yg_tild = savgol_filter(yg, win_size, pol_order,mode='nearest')
        yg_tild_dot = savgol_filter(yg, win_size, pol_order, deriv=1, delta=delta,mode='nearest')
        yg_tild_dot_dot = savgol_filter(yg, win_size, pol_order, deriv=2, delta=delta)
        if show_flag:
            plt.plot(yg_tild)
            plt.show()
            plt.plot(yg_tild_dot)
            plt.show()
            plt.plot(yg_tild_dot_dot)
            plt.show()
        return yg_tild, yg_tild_dot, yg_tild_dot_dot
    else:
        yg_tild = savgol_filter(yg, win_size, pol_order)
        yg_tild_dot = savgol_filter(yg_tild, win_size, pol_order, deriv=1, delta=delta)
        yg_tild_dot_dot = savgol_filter(yg_tild_dot, win_size, pol_order, deriv=1, delta=delta)
        plt.plot(yg_tild)
        plt.show()
        plt.plot(yg_tild_dot)
        plt.show()
        plt.plot(yg_tild_dot_dot)
        plt.show()
        return yg_tild, yg_tild_dot, yg_tild_dot_dot


def far_field_fit(time, yg, volume_params, R_params, rho, k, g=9.81, mass=0, manip_name=None):
    # yg = savgol_filter(yg,11,3)
    Finf = volume_params['Vinf'] * rho * g
    F0 = volume_params['V0'] * rho * g
    deltaF = Finf - F0
    Rinf = R_params['Rinf']
    R0 = R_params['R0']
    deltaR = Rinf - R0
    print deltaF
    params = Parameters()
    Amin = 0.05
    params.add('t0', value=0, vary=False)
    params.add('q0', value=0.1)
    params.add('q0p', value=-0.1)
    params.add('k', value=k, min =k-1, max=k+1vary=False, min=1)
    params.add('massTot', value=0.036, min=0.01,max =1)
    params.add('deltaprim', value=0.1)
    params.add('deltaF', value=deltaF, vary=False)
    params.add('Finf', value=Finf, vary=False)
    params.add('deltaR', value=deltaR, vary=False)
    params.add('Rinf', value=Rinf, vary=False)
    params.add('tau_v', value=volume_params['tau_long'], vary=False)
    params.add('tau_R', value=R_params['tau_long'], vary=False)
    params.add('C1', expr="Finf/k", )
    params.add('C2', expr="deltaF/(massTot*((1/tau_v)**2-(2*lamda/tau_v)+omega0**2))")
    params.add('D1', expr="Rinf", )
    params.add('D2', expr="(omega0**2)*deltaR/((1/tau_R)**2-(2*lamda/tau_R)+omega0**2)")
    params.add('omega0', expr="(k/massTot)**0.5")
    params.add('lamda', value=1)
    params.add('deltaprim', expr='lamda**2-omega0**2')
    result = minimize(residual_far_field, params, args=(time, yg), method='leastsq', maxfev=int(5e4), xtol=1e-8,
                      ftol=1e-8)
    report_fit(result.params, min_correl=0.5)
    print result.success
    v = params.valuesdict()
    yg_th = residual_far_field(params, time, yg, fit=False)
    v["added_mass"] =0
    v['real_mass'] = mass
    v["added_mass"] = v['massTot'] - v['real_mass']

    # yg_f,ygdot_f,ygdotdot_f = filterer(yg,win_size,pol_order,double_filter)
    return v, yg_th


def residual_far_field(parameters, t, yg, fit=True):
    if type(parameters) is not dict:
        v = parameters.valuesdict()
    if v['deltaprim'] > 0:
        v['sh'] = (v['q0p'] + v['lamda'] * v['q0']) / v['deltaprim'] ** 0.5
        yg_th_general = np.exp(-v['lamda'] * (t - v['t0'])) * (
        v['q0'] * np.cosh(v['deltaprim'] ** (0.5) * (t - v['t0'])) + v['sh'] * np.sinh(
            v['deltaprim'] ** (0.5) * (t - v['t0']))) + v['C1'] + v['C2'] * np.exp(-(t) / v['tau_v']) + v['D1'] + v[
                                                                                                                      'D2'] * np.exp(
            -(t) / v['tau_R'])
    elif v['deltaprim'] == 0:
        v['sh'] = v['q0p'] + v['lamda'] * v['q0']
        yg_th_general = np.exp(-v['lamda'] * (t - v['t0'])) * (v['q0'] + v['sh'] * (t - v['t0'])) + v['C1'] + v[
                                                                                                                  'C2'] * np.exp(
            -(t) / v['tau_v']) + v['D1'] + v['D2'] * np.exp(-(t) / v['tau_R'])
    else:
        v['sh'] = (v['q0p'] + v['lamda'] * v['q0']) / (-v['deltaprim']) ** 0.5
        yg_th_general = np.exp(-v['lamda'] * (t - v['t0'])) * (
        v['q0'] * np.cos((-v['deltaprim']) ** (0.5) * (t - v['t0'])) + v['sh'] * np.sin(
            (-v['deltaprim']) ** (0.5) * (t - v['t0']))) + v['C1'] + v['C2'] * np.exp(-(t) / v['tau_v']) + v['D1'] + v[
                                                                                                                         'D2'] * np.exp(
            -(t) / v['tau_R'])
    if fit:
        return yg_th_general - yg
    else:
        return yg_th_general


def far_field(filename, time_cut, v0, yg0,ly0, rho, k, g=9.81,mass=0,manip_name="manip",filter=False):
    time_cut += DELTA_T
    results_delayed_dict = read_results_file(filename, time_cut)
    time_delayed = results_delayed_dict["time"] - results_delayed_dict["time"][0]
    volume_array_normalized = v0 - results_delayed_dict['External_volume']
    yl = ly0-results_delayed_dict['Ylow']
    yg = yg0 - results_delayed_dict['Gravity_Y']
    R = yg-yl

    if filter:
        FILTER_WIN = 11
        print 'filtering\n'
        volume_array_normalized=filterer(volume_array_normalized,FILTER_WIN,3,1)
        volume_array_normalized=volume_array_normalized[0]
        R=filterer(R,FILTER_WIN,3,1)
        R = R[0]
        yg=filterer(yg,FILTER_WIN,3,1)
        yg=yg[0]

    volume_params,v_th = volume_fit(time_delayed, volume_array_normalized,far_field=True)
    R_params,R_th = r_fit(time_delayed, R, volume_params,far_field=True)
    far_field_param, yg_th = far_field_fit(time_delayed, yg, volume_params, R_params, rho, k,mass=mass)
    image_name = filename.split(".t")[0]
    Vlim = (5e-6,2e-5)
    Rlim = (-1e-2,0)
    Yglim = (-0.03,-0.01)
    plot_double(time_delayed, yg, yg_th, "Experimental", "Fit", image_name + "_far_field_fit.png","Yg far field fit","Time(s)","Yg(m)",ylim=Yglim)
    plot_double(time_delayed, volume_array_normalized, v_th, "Experimental", "Fit", image_name + "_Volume.png","Volume fit","Time(s)","V(m^3)",ylim = Vlim)
    plot_double(time_delayed, R,R_th, "Experimental", "Fit", image_name + "_R.png","Radius fit","Time(s)","R(m)",ylim=Rlim)
    plt.title("Far field fit")
    return volume_params,R_params,far_field_param


def rolling(filename):
    results = np.genfromtxt(filename, dtype=None, delimiter="\t", skip_header=1)


def mass_residual(params, ygpp, deltay, k):
    v = params.valuesdict()
    Fv = v['mass'] * ygpp + k * deltay
    Fv = Fv[:, 0]
    # print "["+str(Fv)+","+str(v['mass'])+"]"
    return Fv


def free_oscillation(file, rho, k, mass, R=0.025):
    results_dict = read_results_file(file)
    time = results_dict["time"]
    yg = results_dict['Gravity_Y']
    delta = (-time[0] + time[1000]) / 1000.
    yg, ygp, ygpp = filterer(yg, 101, 3, delta)
    deltay = yg - yg[-1]

    index = np.argwhere(abs(ygp) < 1e-4)
    print index.size
    params = Parameters()
    params.add('mass', value=mass, min=0.001)
    mi = minimize(mass_residual, params, args=(ygpp[index], deltay[index], k))
    report_fit(mi.params, min_correl=0.5)
    v = params.valuesdict()
    Fv = v['mass'] * ygpp + k * deltay
    print  v['mass']
    Fv_array = np.ndarray((Fv.shape[0], 2), dtype=Fv.dtype)
    Fv_array[:, 0] = time
    Fv_array[:, 1] = Fv
    V_array = np.ndarray((Fv.shape[0], 2), dtype=Fv.dtype)
    V_array[:, 0] = time
    V_array[:, 1] = ygp
    return Fv_array, V_array


def manip_analyser(rootfolder, nbrOfExperiments, rho, mass=0, free_oscillation_flag=False,filter = False):
    if not free_oscillation_flag:
        subfolder = rootfolder + "Isotropic_B\\"
        file_list = os.listdir(subfolder)
        textfile_list = []
        for i in xrange(len(file_list)):
            if file_list[i].endswith(".txt"):
                textfile_list.append(file_list[i])
        list_fit = []
        for i in xrange(len(textfile_list)):
            if i < 2 * nbrOfExperiments[0]:
                k = i / nbrOfExperiments[0]
            else:
                k = ((i - (2 * nbrOfExperiments[0])) / nbrOfExperiments[1]) + 2
            r = isotropic(subfolder + textfile_list[i], ext_V_reference[k], ly_reference[k], rho)
            if r is not None:
                list_fit.append(r)
        iso_writer(rootfolder, list_fit)
        k = 0
        i = 0
        for el in list_fit:
            if el[1] < K_THRESHOLD_MAX and el[1] > K_THRESHOLD_MIN:
                k += el[1]
                i += 1
        k = k / float(i)
        print 'Mean stiffness: ' + str(k)

        subfolder = rootfolder + "Buckling\\"
        file_list = os.listdir(subfolder)
        textfile_list = []
        V_list=[]
        R_list = []
        far_field_list = []

        for i in xrange(len(file_list)):
            if file_list[i].endswith(".txt"):
                textfile_list.append(file_list[i])
        for i in xrange(0, nbrOfExperiments[0]):
            v,r,y=far_field(subfolder + textfile_list[i], timStartArraybuckling[i], ext_V_reference[i], Yg_reference[i],ly_reference[i], rho,
                     k, manip_name=textfile_list[i].split(".t")[0],mass=mass,filter= filter)
            V_list.append(v)
            R_list.append(r)
            far_field_list.append(y)

        parameters_writer(rootfolder+"Volume_fit.txt",V_list)
        parameters_writer(rootfolder+"R_fit.txt",R_list)
        parameters_writer(rootfolder+"Yg_fit.txt",far_field_list)

    if free_oscillation_flag:
        k = 8.2
        array2Ds_list_x = []
        array2Ds_list_y = []
        array_v_list = []
        array_v2_list = []
        subfolder_in = rootfolder + "\\Input_fit\\"
        subfolder_out = rootfolder + "\\Output_fit\\"
        file_list = os.listdir(subfolder_in)
        textfile_list = []
        for i in xrange(len(file_list)):
            if file_list[i].endswith(".txt"):
                textfile_list.append(file_list[i])
        for i in xrange(len(file_list)):
            array2D, array2D_v = free_oscillation(subfolder_in + textfile_list[i], rho, k, mass[i])
            array2Ds_list_x.append(array2D[:, 0])
            array2Ds_list_y.append(array2D[:, 1])
            array_v_list.append(array2D_v[:, 1])
            array_v2_list.append(array2D_v[:, 1] ** 2)
            output_filename = textfile_list[i].split(".t")[0]
            manip_name = output_filename.split("_")
            manip_name = '_'.join(manip_name[0:1])
            output_filename_g = subfolder_out + output_filename

            plot_unique(array2D[:, 0], array2D[:, 1], output_filename_g + "Fv(t)" + ".PNG", manip_name, "r-",
                        "Fv(t) " + manip_name, "time(s)", "Fv(N)", show_flag=False)
            plot_unique(array2D_v[:, 1], array2D[:, 1], output_filename_g + "Fv(V)" + ".PNG", manip_name, "b^",
                        "Fv(V) " + manip_name, "V(m/s)", "Fv(N)", show_flag=False)
            plot_unique(array2D_v[:, 1] ** 2, array2D[:, 1], output_filename_g + "Fv(V^2)" + ".PNG", manip_name, "bo",
                        "Fv(V^2) " + manip_name, "V^2(m2/s2)", "Fv(N)", show_flag=False)
            plot_unique(array2D[:, 0], array2D_v[:, 1], output_filename_g + "V(t)" + ".PNG", manip_name, "b*",
                        "V(t) " + manip_name, "time(s)", "V(m/s)", show_flag=False)

        plot_multiple(array2Ds_list_x, array2Ds_list_y, 3, subfolder_out + "all_Fv(t).PNG",
                      "F_v(t) for the different balls", "time(s)", "Fv(N)", show_flag=False)
        plot_multiple(array_v_list, array2Ds_list_y, 3, subfolder_out + "all_Fv(V).PNG",
                      "F_v(V) for the different balls", "V(m/s)", "Fv(N)", show_flag=False)
        plot_multiple(array_v2_list, array2Ds_list_y, 3, subfolder_out + "all_Fv(V^2).PNG",
                      "F_v(V^2) for the different balls", "V^2(m2/s2)", "Fv(N)", show_flag=False)
        plot_multiple(array2Ds_list_x, array_v_list, 3, subfolder_out + "all_V(t).PNG", "V(t) for the different balls",
                      "time(s)", "V(m/s", show_flag=False)


if __name__ == "__main__":
    Yg_reference = [2.57520e-02, 2.58761e-02, 2.59152e-02, 2.57457e-02, 2.57457e-02]
    ly_reference = [5.10853e-02, 5.10960e-02, 5.12694e-02, 5.10834e-02, 5.10834e-02]
    ext_V_reference = [6.54484e-05, 6.54422e-05, 6.54462e-05, 6.54477e-05, 6.54480e-05]
    timeStartbuckling = [1520, 1120, 173, 1921, 535]
    timStartArraybuckling = np.array(timeStartbuckling) * 0.0002
    mass_ball = [0.017, 0.034, 0.042]
    # #2
    # ly_reference = [5.24448e-02,5.22848e-02,5.27078e-02,5.27705e-02,5.28126e-02]
    # ext_V_reference =[6.54160e-05,6.54171e-05,6.54485e-05,6.52726e-05,6.54427e-05]
    # timeStartbuckling = [1660,1680,1052,1130,160]
    # timStartArraybuckling = np.array(timeStartbuckling)*0.0002
    # #
    # # 6.5
    # Yg_reference = [2.63940e-02, 2.75485e-02, 2.76895e-02, 2.72654e-02, 2.74623e-02]
    # ly_reference = [5.14205e-02, 5.25599e-02, 5.27431e-02, 5.23628e-02, 5.25493e-02]
    # ext_V_reference = [6.54476e-05, 6.54452e-05, 6.54495e-05, 6.54458e-05, 6.54447e-05]
    # timeStartbuckling = [1670, 1780, 1555, 1100, 1840]
    # timStartArraybuckling = np.array(timeStartbuckling) * 0.0002

    filename = 'D:\\Documents\\Footage\\videos\\Manip_finales\\Glycerol_100\\Ressort\\5_25\\output_general\\'
    # filename = 'D:\\Documents\\Footage\\videos\\Manip_finales\\Glycerol_100\\Ressort\\Oscillations_libres\\'
    print filename
    rho = 1180
    manip_analyser(filename, (5, 4), rho, mass=mass_ball[1], free_oscillation_flag=False,filter = False)
    # D:\Documents\Footage\videos\Manip_finales\Glycerol_100\Ressort\Oscillations_libres
