__author__ = 'Adel'
import numpy as np
import lmfit
import matplotlib.pyplot as plt

PATH_FIT_B = "./Input/13_03_2015/Buckling/Ball_shape/results_fit.txt"
PATH_VENTOUSE_B= "./Input/13_03_2015/Buckling/Ventouse/results_ventouse.txt"


PRESSURE_THRESH = 50000.0
THICKNESS = 0.00429
INITIAL_DIAMETER = 0.05135
MASS = 0.039
GRAVITY = 9.81
PATM = 101325
STIFFNESS = 7.46
INDEX_BUCKLING = 930
INDEX_BUCKLING_2 = 100000
INDEX_DEBUCKLING = 100000
SINGULAR = True
OUTPUT_QS = './Output/2015_03_13/Debuckling/'
OUTPUT_G  = './Output/2015_03_13/'
DATE = '2015_03_13'
font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
R0_REF = 51.35/1000.0
R0PIX_REF = 400.0
CONV = R0_REF/R0PIX_REF
D = 4.29/1000.0

A_ERROR_D = 0.06/1000.0
A_ERROR_R_REF = 0.01/1000.0
A_ERROR_R_PIX = 1.0
A_ERROR_Pext = 500.0
RESOLUTION_ERROR = 50


def read_ventouse_file(path,multiple=None,conversion = 1.0):
    ventouse_data =  {}
    t = []
    p = []
    x = []
    if multiple is None:
        fileData = open(path, "r")
        line = fileData.readline()

        data = fileData.readline().rstrip("\n").split("\t")
        while data[0] != "":
            t.append(float(data[0]))
            p.append(float(data[1]))
            x.append(float(data[2]))
            data = fileData.readline().rstrip("\n").split("\t")
    else:
        t = []
        p = []
        x = []
        ventouse_data =  {}
        for pth in multiple:
            fileData = open(path+pth, "r")
            line = fileData.readline()
            data = fileData.readline().rstrip("\n").split("\t")
            while data[0] != "":
                t.append(float(data[0]))
                p.append(float(data[1]))
                x.append(float(data[2]))
                data = fileData.readline().rstrip("\n").split("\t")
    t = np.array(t)
    p = np.array(p)
    x = np.array(x)*conversion
    ventouse_data ['t'] = t
    ventouse_data ['p'] = p
    ventouse_data ['x'] = x
    return ventouse_data


def read_fit_file(path, multiple=None,conversion = 1.0):
    fit_data =  {}
    t = []
    p = []
    max_height = []
    max_width = []
    max_width_theta = []
    x_low = []
    y_low = []
    gravity_x = []
    gravity_y = []
    external_volume = []
    if multiple is None:
        fileData = open(path, "r")
        line = fileData.readline()
        data = fileData.readline().rstrip("\n").split("\t")
        while data[0] != "":
            t.append(float(data[0]))
            p.append(float(data[1]))
            max_height.append(float(data[2]))
            max_width.append(float(data[3]))
            max_width_theta.append(float(data[4]))
            x_low.append(float(data[5]))
            y_low.append(float(data[6]))
            gravity_x.append(float(data[7]))
            gravity_y.append(float(data[8]))
            external_volume.append(float(data[9]))
            data = fileData.readline().rstrip("\n").split("\t")
    else:
        for pth in multiple:
            fileData = open(path+pth, "r")
            line = fileData.readline()
            data = fileData.readline().rstrip("\n").split("\t")
            while data[0] != "":
                t.append(float(data[0]))
                p.append(float(data[1]))
                max_height.append(float(data[2]))
                max_width.append(float(data[3]))
                max_width_theta.append(float(data[4]))
                x_low.append(float(data[5]))
                y_low.append(float(data[6]))
                gravity_x.append(float(data[7]))
                gravity_y.append(float(data[8]))
                external_volume.append(float(data[9]))
                data = fileData.readline().rstrip("\n").split("\t")

    t = np.array(t)
    p = np.array(p)
    max_height= np.array(max_height)*conversion
    max_width= np.array(max_width)*conversion
    max_width_theta= np.array(max_width_theta)
    x_low= np.array(x_low)*conversion
    y_low= np.array(y_low)*conversion
    gravity_x= np.array(gravity_x)*conversion
    gravity_y= np.array(gravity_y)*conversion
    external_volume= np.array(external_volume)*conversion**3
    fit_data ['t'] = t
    fit_data ['p'] = p
    fit_data ['max_height'] = max_height
    fit_data ['max_width'] = max_width
    fit_data ['max_width_theta'] = max_width_theta
    fit_data['x_low'] = x_low
    fit_data['y_low'] = y_low
    fit_data['gravity_x'] = gravity_x
    fit_data['gravity_y'] = gravity_y
    fit_data['external_volume'] = external_volume
    return fit_data


def compute_rt__abserror(max_width):
    return (CONV*A_ERROR_R_PIX + ((CONV*max_width)/R0PIX_REF)*A_ERROR_R_REF+((CONV*max_width)*R0_REF/R0PIX_REF**2)*A_ERROR_R_PIX)


def compute_vext_abserror(max_width):
    rt__abserror = compute_rt__abserror(max_width)
    vext_abserror = 4.*np.pi*max_width**2 *rt__abserror
    return vext_abserror

def compute_shV_abserror():
    R0_E = 4*np.pi*(R0_REF**2+(R0_REF-D)**2)*A_ERROR_R_REF
    D_E =4*np.pi*(R0_REF-D)**2*A_ERROR_D
    shV_abserror = 4*np.pi*((R0_REF**2+(R0_REF-D)**2)*A_ERROR_R_REF+(R0_REF-D)**2*A_ERROR_D)
    return shV_abserror

def compute_vmedian_relerror(max_width,vmedian):
    sv_error = compute_shV_abserror()
    # vext_error = compute_vext_abserror(max_width)
    # vmedian_t_error = vext_error + 0.5*sv_error
    vmedian0 = vmedian[0]+0.000001
    deltaVM = 0
    return deltaVM

def compute_deltaP_abserror(max_width,vint,pint):
    sv_error = compute_shV_abserror()
    vext_error = compute_vext_abserror(max_width)
    vint_error = vext_error+sv_error
    A =(vint_error[0]/vint[0])+(vint_error/vint)
    pint_error = A/4*pint
    deltaP_error = A_ERROR_Pext+pint_error
    return deltaP_error
def parse_data (index, fit_data,ventouse_data):

    quasi_static = {}
    dynamic = {}
    for k in fit_data.keys():
        quasi_static[k] = fit_data[k][0:index]
        dynamic [k] = fit_data[k][index:-1]
    for k in ventouse_data.keys():
        quasi_static[k] = ventouse_data[k][0:index]
        dynamic [k] = ventouse_data[k][index:-1]
    return quasi_static,dynamic


def residLinearK(parameters,p,x, externalVolume, mass,name, gravity = 9.81):
    v = parameters.valuesdict()
    deltaX_theorique = -(1./v['k'])*(((1000*externalVolume*gravity)-(mass*gravity))-((1000*externalVolume[0]*gravity)-(mass*gravity)))
    resid = (x-x[0])-deltaX_theorique
    t = np.arange(0,deltaX_theorique.shape[0],1)
    plt.plot(externalVolume,x-x[0],'-r',label = r'$\Delta X$ experimental '+name)
    plt.plot(externalVolume,deltaX_theorique,'ob',label = r'$\Delta X$ theoretical')
    plt.title('Evolution of '+r'$\Delta X$ '+name+'with pressure')
    plt.xlabel('Volume (m3)')
    plt.ylabel(r'$displacement$')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(OUTPUT_QS+'k_quasi_static_fit'+name+DATE+'.png', bbox_inches='tight')

    plt.clf()
    return resid


def residLinearE(parameters,dpc, dvm, name):
    v = parameters.valuesdict()
    deltaV_theorique = (1./v['E'])*dpc
    resid = dvm-deltaV_theorique
    t = np.arange(0,deltaV_theorique.shape[0],1)
    plt.plot(dpc,dvm,'-r',label = r'$\Delta V/V$ experimental')
    plt.plot(dpc,deltaV_theorique,'ob',label = r'$\Delta V/V$ theoretical')
    plt.title('Evolution of volume with '+r'$\Delta Pc$')
    plt.xlabel(r'$\Delta P$'+' (Pa)')
    plt.ylabel(r'$\Delta V/V$')
    plt.xticks(np.arange(0.0,80000.0,10000.0))
    plt.yticks(np.arange(0.0,0.2,0.05))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(OUTPUT_QS+'dv_dp_quasi_static_fit'+name+DATE+'.png', bbox_inches='tight')

    plt.clf()
    return resid


def residFree(parameters,t,x):
    v = parameters.valuesdict()
    theoretical_X = v["A"]+v["B"]*np.exp(-(t-v["t0"])/v["tau_ff"])*np.cos(v["w"]*(t-v['t0']))
    return x-theoretical_X



def residBuckling(parameters,t,x):
    v = parameters.valuesdict()
    theoretical_X = v["A"]+v["B"]*np.exp(-(t-v["t0"])/v["tau_ff"])*np.cos(v["w"]*(t-v['t0']))-v["C"]*np.exp(-(t-v["t0"])/v["tau_vol"])
    return x-theoretical_X


def compute_fit(t, p, x, free=False, t0=0.0):

    poscillation = lmfit.Parameters()
    poscillation.add('A', value=10.0, min=0.0)
    poscillation.add('B', value=10.0, min=0.0)
    poscillation.add('t0', value=t0)
    poscillation.add('tau_ff', value=10.0, min=0.0)
    poscillation.add('w', value=10.0, min=0.0)
    poscillation.add('C', value=10.0, min=0.0)
    poscillation.add('tau_vol', value=10.0, min=0.0)

    if free:
        mi_oscillation = lmfit.minimize(residFree,poscillation,args=(t,x))
        vDict = poscillation.valuesdict()
        X = vDict["A"]+(vDict["B"]*np.exp(-(t-vDict["t0"])/vDict["tau_ff"])*np.cos(vDict["w"]*(t-vDict['t0'])))
        Xpoint = -1*vDict["B"]*np.exp(-(t-vDict["t0"])/vDict["tau_ff"])*((1/vDict["tau_ff"])*np.cos(vDict["w"]*(t-vDict['t0']))+vDict["w"]*np.sin(vDict["w"]*(t-vDict['t0'])))
    else:
        mi_oscillation = lmfit.minimize(residBuckling,poscillation,args=(t,x))
        vDict = poscillation.valuesdict()
        X = vDict["A"]+(vDict["B"]*np.exp(-(t-vDict["t0"])/vDict["tau_ff"])*np.cos(vDict["w"]*(t-vDict['t0'])))-vDict["C"]*np.exp(-(t-vDict["t0"])/vDict["tau_vol"])
        Xpoint = -1*vDict["B"]*np.exp(-(t-vDict["t0"])/vDict["tau_ff"])*((1/vDict["tau_ff"])*np.cos(vDict["w"]*(t-vDict['t0']))+vDict["w"]*np.sin(vDict["w"]*(t-vDict['t0'])))

    return mi_oscillation, poscillation, X, Xpoint


def deltaVmedian(externalVolume,shellVolume,maxwidth,V0=True,givenVolume=None):
    medianVolume = externalVolume - shellVolume/2.
    if givenVolume is None:
        if V0:
            deltaVmedian = (medianVolume[0]-medianVolume)/medianVolume[0]
            A_error_deltaVmedian =  compute_vmedian_relerror(maxwidth,medianVolume)*deltaVmedian
            return deltaVmedian,medianVolume[0],A_error_deltaVmedian
        else:
            deltaVmedian = (medianVolume[-1]-medianVolume)/medianVolume[-1]
            A_error_deltaVmedian =  compute_vmedian_relerror(maxwidth,medianVolume)*deltaVmedian
            return deltaVmedian,medianVolume[-1],A_error_deltaVmedian
    else:
        deltaVmedian = (givenVolume-medianVolume)/givenVolume
        A_error_deltaVmedian =  compute_vmedian_relerror(maxwidth,medianVolume)*deltaVmedian
        return deltaVmedian,givenVolume,A_error_deltaVmedian



def deltaP (p,externalVolume,shellVolume,maxwidth,patm= 103000,V0=True,givenVolume=None):
    p_ext = p+ patm
    R_error_pext = A_ERROR_Pext/p

    internalVolume = externalVolume - shellVolume
    if givenVolume is None:
        if V0:
            p_int = patm* internalVolume[0]/internalVolume
            v0 = internalVolume[0]
        else:
            p_int = patm* internalVolume[-1]/internalVolume
            v0 = internalVolume[-1]


    else:
        p_int = patm* givenVolume/internalVolume
        v0 = givenVolume
    # A_error_P = compute_deltaP_abserror(maxwidth,internalVolume,p_int)
    A_error_P = 0
    return p_ext-p_int,v0,A_error_P
def shellVolume(thickness,initialDiameter,externalVolumeT0):
    shellv = externalVolumeT0 - ((4.*np.pi/3.)*((initialDiameter/2)-thickness)**3)
    return shellv


def multiple_comparison(listToFiles):
    vs = []
    dpcList=[]
    dvmList=[]
    xList = []
    y_lowList = []
    gravity_yList = []
    for e in listToFiles:
        fit_data = read_fit_file(e['path_fit'],e['sub_pathes_fit'])
        ventouse_data = read_ventouse_file(e['path_ventouse'],e['sub_pathes_ventouse'])
        quasi_static, dynamic = parse_data(e['index_buckling'],fit_data,ventouse_data)
        sv = shellVolume(THICKNESS,INITIAL_DIAMETER,quasi_static['external_volume'][0])
        dvm = deltaVmedian(quasi_static['external_volume'],sv)
        dp = deltaP(quasi_static['p'],quasi_static['external_volume'],sv)
        vs.append(quasi_static['external_volume'])
        xList.append(quasi_static['x'])
        y_lowList.append(quasi_static['y_low'])
        gravity_yList.append(quasi_static['gravity_y'])
        dpcList.append(dpc)
        dvmList.append(dvm)
    plt.title('Evolution of volume with '+r'$\Delta Pc$')
    plt.xlabel('Pressure (Pa)')
    plt.ylabel(r'$\Delta V/V$')
    colors = ['b','g','r','c','m','y','k','w']
    linestyle = ['o','^','-']
    for i in xrange(0,len(dvmList)):
        plt.plot(dpcList[i],dvmList[i],colors[i%8]+linestyle[i%3],label = r'$\Delta V/V = f(\Delta Pc)$ '+listToFiles[i]['name'])
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('./Output/comparison_dv_dp_quasi_static.png', bbox_inches='tight')
    plt.show()
    plt.clf()

    plt.title('Comparison of all displacements (gravity_y) in function of '+r'$\Delta Pc$')
    plt.xlabel('Volume (m3)')
    plt.ylabel('displacement '+r'$\Delta X$'+' (m)')
    for i in xrange(0,len(vs)):
        plt.plot(vs[i],gravity_yList[i]-gravity_yList[i][0],colors[i%8]+linestyle[i%3],label = r'$\Delta V/V = f(\Delta Pc)$ '+listToFiles[i]['name'])
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('./Output/comparison_all_displacements_gravity_y.png', bbox_inches='tight')
    plt.show()
    plt.clf()
    plt.title('Comparison of all displacements (y_low) in function of '+r'$\Delta Pc$')
    plt.xlabel('Volume (m3)')
    plt.ylabel('displacement '+r'$\Delta X$'+' (m)')
    for i in xrange(0,len(vs)):
        plt.plot(vs[i],y_lowList[i]-y_lowList[i][0],colors[i%8]+linestyle[i%3],label = r'$\Delta V/V = f(\Delta Pc)$ '+listToFiles[i]['name'])
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('./Output/comparison_all_displacements_y_low.png', bbox_inches='tight')
    plt.show()
    plt.clf()
    plt.title('Comparison of all displacements (x) in function of '+r'$\Delta Pc$')
    plt.xlabel('Volume (m3)')
    plt.ylabel('displacement '+r'$\Delta X$'+' (m)')
    for i in xrange(0,len(vs)):
        plt.plot(vs[i],xList[i]-xList[i][0],colors[i%8]+linestyle[i%3],label = r'$\Delta V/V = f(\Delta Pc)$ '+listToFiles[i]['name'])
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('./Output/comparison_all_displacements_x.png', bbox_inches='tight')
    plt.show()
    plt.clf()
def residLinearX0p(parameters,x,p):
    v = parameters.valuesdict()
    theoreticalx0p = v['x0b']+ v['alpha']*p
    return x-theoreticalx0p
def conversion_corrector(diameter_before,diameter_after,conversion_ratio):
    new_conversion_ratio = (diameter_after/diameter_before)*conversion_ratio
    return new_conversion_ratio/conversion_ratio


plt.rc('figure',figsize=(16,12))
if __name__ == "__main__":
    if SINGULAR:
        sub_pathes_fit = ['results_fit_p0.txt','results_fit_p100.txt','results_fit_p200.txt','results_fit_p300.txt','results_fit_p400.txt','results_fit_p500.txt','results_fit_p600.txt','results_fit_p700.txt','results_fit_p800.txt']
        sub_pathes_ventouse = ['results_ventouse_p0.txt','results_ventouse_p100.txt','results_ventouse_p200.txt','results_ventouse_p300.txt','results_ventouse_p400.txt','results_ventouse_p500.txt','results_ventouse_p600.txt','results_ventouse_p700.txt','results_ventouse_p800.txt']
        conv = 1.0
        fit_data_b = read_fit_file(PATH_FIT_B,multiple=None)
        ventouse_data_b = read_ventouse_file(PATH_VENTOUSE_B,multiple=None)
        

        
        # get the moment when we have the buckling
        # Pressure


        # separation of phases
        quasi_static_b, dynamic_b = parse_data(INDEX_BUCKLING,fit_data_b,ventouse_data_b)
        # plt.plot(dynamic_b2['t'],dynamic_b2['x'])
        # plt.show()


        sv = shellVolume(THICKNESS,INITIAL_DIAMETER,quasi_static_b['external_volume'][0])
        dvm_b,vm0,A_ERROR_V = deltaVmedian(quasi_static_b['external_volume'],sv,quasi_static_b['max_width'],V0 = True)
        dp_b,vi0,A_ERROR_P = deltaP(quasi_static_b['p'],quasi_static_b['external_volume'],sv,quasi_static_b['max_width'],V0=True)

        begi = 500
        endi =800
        begi = 0
        endi =-1

        ref_len = quasi_static_b['x']+((quasi_static_b['external_volume']*1000*GRAVITY-MASS*GRAVITY)/STIFFNESS)
        
        
        
        # plt.plot(quasi_static_db['t'],Xdb_corrected,'-b',label='Xdb(t)')
        # plt.plot(quasi_static_b2['t'],Xb_corrected,'-r',label='Xb(t)')
        
        

        

        dvm_b,vm0,A_ERROR_V = deltaVmedian(quasi_static_b['external_volume'],sv,quasi_static_b['max_width'],V0 = True)
        dp_b,vi0,A_ERROR_P = deltaP(quasi_static_b['p'],quasi_static_b['external_volume'],sv,quasi_static_b['max_width'],V0=True)
        
        
       
        # plt.plot(quasi_static_db['t'],quasi_static_db['external_volume'],'-b',label='Vdb')
        # plt.plot(fit_data_b2['t'],V_b2,'or',label='Vb2')
        # plt.plot(quasi_static_b['t'],quasi_static_b['external_volume'],'-g',label='Vb')
        # plt.title('Vb + vdb')
        # plt.xlabel('time')
        # plt.ylabel('V')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        # plt.savefig(OUTPUT_G+'v_dp'+DATE+'.png', bbox_inches='tight')
        # plt.show()
        # plt.clf()


        


        dvm =dvm_b
        dp = dp_b
        # plinearE = lmfit.Parameters()
        # plinearE.add('E', value = 100000, )
        # print 'E 9P'
        # mi_linear = lmfit.minimize(residLinearE,plinearE,args=(dp,dvm,'9p'))
        # lmfit.printfuncs.report_fit(plinearE)
        dv_dp = open(OUTPUT_QS+'dv_dp.txt','w')
        dv_dp.writelines("i"+"\t"+"dp"+"\t"+"dvm")
        for i in xrange(dvm.shape[0]):
            dv_dp.writelines("\n")
            dv_dp.writelines("{0:.5e}".format(i)+'\t')
            dv_dp.writelines("{0:.5e}".format(dp[i])+'\t')
            # dv_dp.writelines("{0:.5e}".format(A_ERROR_P[i])+'\t')
            dv_dp.writelines("{0:.5e}".format(dvm[i])+'\t')
            # dv_dp.writelines("{0:.5e}".format(A_ERROR_V[i])+'\t')
        # XERR = A_ERROR_P
        plt.plot(dvm,dp,'-or',label = r'$\Delta V/V = f(\Delta Pc)$')
        # plt.errorbar(dp[index_resolution], dvm[index_resolution],xerr=XERR[index_resolution], yerr=A_ERROR_V[index_resolution], fmt='o')
        plt.title('Evolution of volume with '+r'$\Delta Pc$')
        plt.ylabel('Pressure (Pa)')
        plt.xlabel(r'$\Delta V/V$')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(OUTPUT_QS+'dp_dv_150615'+'.png', bbox_inches='tight')
        plt.show()
        plt.clf()

        # plt.plot(dp,A_ERROR_P,'-r',label = r'$ERROR P = f(\Delta Pc)$')
        # plt.xlabel('Pressure (Pa)')
        # plt.ylabel('ERROR P')
        # plt.legend()
        # plt.show()
        # plt.clf()
        #
        # dvm = deltaVmedian(fit_data['external_volume'],sv)
        # dpc = deltaPc(fit_data['p'],fit_data['external_volume'],sv)
        # dv_dp = open(OUTPUT_G+'dv_dp.txt','w')
        # dv_dp.writelines("i"+"\t"+"dpc"+"\t"+"dvm")
        # for i in xrange(dvm.shape[0]):
        #     dv_dp.writelines("\n")
        #     dv_dp.writelines("{0:.5e}".format(i)+'\t')
        #     dv_dp.writelines("{0:.5e}".format(dpc[i])+'\t')
        #     dv_dp.writelines("{0:.5e}".format(dvm[i])+'\t')
        # plt.plot(dpc,dvm,'-r',label = r'$\Delta V/V = f(\Delta Pc)$')
        # plt.title('Evolution of volume with '+r'$\Delta Pc)$')
        # plt.xlabel('Pressure (Pa)')
        # plt.ylabel(r'$\Delta V/V$')
        # plt.legend()
        # plt.savefig(OUTPUT_G+'dv_dp_total.png', bbox_inches='tight')
        # plt.clf()
        # Fit Sherical
        # sv = shellVolume(THICKNESS,INITIAL_DIAMETER,spherical_volume_approximation[0])
        # dvm = deltaVmedian(spherical_volume_approximation,sv)
        # dpc = deltaPc(quasi_static['p'],spherical_volume_approximation,sv)
        # plinearE = lmfit.Parameters()
        # plinearE.add('E', value = 100000, min = 10000)
        # print 'E spherical'
        # mi_linear = lmfit.minimize(residLinearE,plinearE,args=(dpc,dvm,'spherical'))
        # lmfit.printfuncs.report_fit(plinearE)
        #
        # plt.plot(dpc,dvm,'-r',label = r'$\Delta V/V = f(\Delta Pc)$')
        # plt.title('Evolution of volume with '+r'$\Delta Pc$ (using spherical approximation)')
        # plt.xlabel('Pressure (Pa)')
        # plt.ylabel(r'$\Delta V/V$')
        # plt.legend()
        # plt.savefig(OUTPUT_QS+'dv_dp_quasi_static_spherical.png', bbox_inches='tight')
        # plt.clf()
    else:
        dict_13_03 = {'path_fit':'./Input/13_03_2015/Buckling/Ball_shape/results_fit.txt','path_ventouse':'./Input/13_03_2015/Buckling/Ventouse/results_ventouse.txt','name':'13_03_2015','index_buckling':900,'sub_pathes_fit':None,'sub_pathes_ventouse':None}
        dict_31_03 = {'path_fit':'./Input/31_03_2015/Ball_shape/results_fit.txt','path_ventouse':'./Input/31_03_2015/Ventouse/results_ventouse.txt','name':'31_03_2015','index_buckling':1000,'sub_pathes_fit':None,'sub_pathes_ventouse':None}
        dict_07_04 = {'path_fit':'./Input/07_03_1315/Ball_shape/results_fit.txt','path_ventouse':'./Input/07_03_1315/Ventouse/results_ventouse.txt','name':'07_03_1315','index_buckling':850,'sub_pathes_fit':None,'sub_pathes_ventouse':None}
        sub_pathes_fit = ['results_fit_p0.txt','results_fit_p100.txt','results_fit_p200.txt','results_fit_p300.txt','results_fit_p400.txt','results_fit_p500.txt','results_fit_p600.txt','results_fit_p700.txt']
        sub_pathes_ventouse = ['results_ventouse_p0.txt','results_ventouse_p100.txt','results_ventouse_p200.txt','results_ventouse_p300.txt','results_ventouse_p400.txt','results_ventouse_p500.txt','results_ventouse_p600.txt','results_ventouse_p700.txt']
        dict_09_04 = {'path_fit':'./Input/09_03_1315/Buckling/Ball_shape/','path_ventouse':'./Input/09_03_1315/Buckling/Ventouse/','name':'09_03_1315','index_buckling':850,'sub_pathes_fit':sub_pathes_fit,'sub_pathes_ventouse':sub_pathes_ventouse}
        sub_pathes_fit = ['results_fit_p0.txt','results_fit_p100.txt','results_fit_p200.txt','results_fit_p400.txt','results_fit_p500.txt','results_fit_p600.txt','results_fit_p700.txt']
        sub_pathes_ventouse = ['results_ventouse_p0.txt','results_ventouse_p100.txt','results_ventouse_p200.txt','results_ventouse_p400.txt','results_ventouse_p500.txt','results_ventouse_p600.txt','results_ventouse_p700.txt']
        dict_13_03 = {'path_fit':'./Input/13_03_2015/Buckling/Ball_shape/','path_ventouse':'./Input/13_03_2015/Buckling/Ventouse/','name':'13_03_2015','index_buckling':850,'sub_pathes_fit':sub_pathes_fit,'sub_pathes_ventouse':sub_pathes_ventouse}
        list_To_Files = []
        list_To_Files.append(dict_13_03)
        list_To_Files.append(dict_31_03)
        list_To_Files.append(dict_07_04)
        list_To_Files.append(dict_09_04)
        list_To_Files.append(dict_13_03)
        multiple_comparison(list_To_Files)
__author__ = 'adel'
