__author__ = 'adel'
import numpy as np
import matplotlib.pyplot as plt
PATH_IN = './Output/2015_04_20/Buckling/Dynamic/fitting_parameters2015_04_20_raw.txt'
PATH_OUT = './Output/2015_04_20/Buckling/Dynamic/physical_parameters2015_04_20_raw.txt'
STIFFNESS = 7.46
SHAPE_CHANGING  = True
OUTPUT_D = './Output/2015_04_20/Buckling/Dynamic/'
# plt.rc('figure',figsize=(16,12))

def read_output(path_in):

    f = open(PATH_IN,'r')
    if SHAPE_CHANGING:

        As=[]
        Bs=[]
        Cs=[]
        Tau_ffs=[]
        Tau_vols=[]
        Omegas = []
        Phis = []
        Betas =[]
        T1s=[]
        Time_ranges = []
        f.readline()
        data = f.readline().rstrip("\n").split("\t")
        while data[0] != "":
            Time_ranges.append((data[0]))
            As.append(float(data[1]))
            Phis.append(float(data[2]))
            Tau_ffs.append(float(data[3]))
            Omegas.append(float(data[4]))
            Tau_vols.append(float(data[5]))
            Cs.append(float(data[7]))
            Bs.append(float(data[6]))
            Betas.append(float(data[8]))
            T1s.append(float(data[9]))
            data = f.readline().rstrip("\n").split("\t")
        As = np.array(As)
        Bs = np.array(Bs)
        Cs = np.array(Cs)
        Omegas = np.array(Omegas)
        Tau_ffs = np.array(Tau_ffs)
        Tau_vols =np.array(Tau_vols)
        Phis = np.array(Phis)
        Betas = np.array(Betas)
        T1s = np.array(T1s)
        dict_input = {'A':As,'B':Bs,'C':Cs,'Omega':Omegas,'Phi':Phis,'tau_ff':Tau_ffs,'tau_vol':Tau_vols,'time_range':Time_ranges,'beta':Betas,'T1':T1s}
        return dict_input
    else:
        As=[]
        Bs=[]
        Tau_ffs=[]
        Omegas = []
        Phis = []
        Time_ranges = []
        f.readline()
        data = f.readline().rstrip("\n").split("\t")
        while data[0] != "":
            Time_ranges.append((data[0]))
            As.append(float(data[2]))
            Phis.append(float(data[3]))
            Tau_ffs.append(float(data[4]))
            Omegas.append(float(data[5]))
            Bs.append(float(data[1]))
            data = f.readline().rstrip("\n").split("\t")
        As = np.array(As)
        Bs = np.array(Bs)
        Omegas = np.array(Omegas)
        Tau_ffs = np.array(Tau_ffs)
        Phis = np.array(Phis)
        dict_input = {'A':As,'B':Bs,'Omega':Omegas,'Phi':Phis,'tau_ff':Tau_ffs,'time_range':Time_ranges,}
        return dict_input
def process_output(dict_input,path_out):
    dict_processed = {}

    if SHAPE_CHANGING:
        dict_processed['time_range']= dict_input['time_range']
        dict_processed['Omega0'] =np.sqrt(dict_input['Omega']**2 + 1./dict_input['tau_ff']**2)
        dict_processed['m_eff'] = (STIFFNESS/dict_processed['Omega0']**2)
        dict_processed['tau_ff']= dict_input['tau_ff']
        dict_processed['Omega']= dict_input['Omega']
        dict_processed['Phi'] = dict_input['Phi']
        dict_processed['G1'] = dict_input['B']*STIFFNESS
        dict_processed['G2'] = dict_input['C']*(dict_processed['m_eff']*((1./dict_input['tau_vol']**2)-(2./(dict_input['tau_vol']*dict_input['tau_ff'])))+STIFFNESS)
        dict_processed['tau_vol']= dict_input['tau_vol']
        dict_processed['F0']= dict_input['beta']*dict_processed['m_eff']
        dict_processed['T1']= dict_input['T1']
    else:
        dict_processed['time_range']= dict_input['time_range']
        dict_processed['Omega0'] =np.sqrt(dict_input['Omega']**2 + 1./dict_input['tau_ff']**2)
        dict_processed['m_eff'] = (STIFFNESS/dict_processed['Omega0']**2)
        dict_processed['tau_ff']= dict_input['tau_ff']
        dict_processed['Omega']= dict_input['Omega']
        dict_processed['Phi'] = dict_input['Phi']

    write_output(dict_processed,path_out)
    return dict_processed
def write_output(dict_process, path_out):
    f = open(path_out,'w')
    f.writelines(["%s\t" % item  for item in list(dict_process.keys())])
    for k in xrange(dict_process['tau_ff'].shape[0]):
        f.writelines("\n")
        for val in dict_process:
            if val!= 'time_range':
                f.writelines("{0:.3e}".format(dict_process[val][k])+"\t")
            else:
                f.writelines(dict_process[val][k]+"\t")
def plot_deltaV(dict_process):
    time = np.arange(0.0,20.0,0.01)
    archimedes = []
    for i in xrange(0,dict_process['G1'].shape[0]):
        delta_Farchimed = dict_process['G1'][i]+dict_process['G2'][i]*(np.exp(-time/dict_process['tau_vol'][i]))
        archimedes.append(delta_Farchimed)
        plt.plot(time,delta_Farchimed,label = 'F Archimedes for tau_vol= '+str(dict_process['tau_vol'][i]))
        plt.xlabel('time (s)')
        plt.ylabel('F archimedes (N)')
        plt.title ('F archimedes (t) '+dict_process['time_range'][i])
        plt.legend()
        plt.savefig(OUTPUT_D+'farchimedes_t'+dict_process['time_range'][i]+'.png',bbox_inches='tight')
        plt.show()
        plt.clf()
    for i in xrange(len(archimedes)):
        plt.plot(time,archimedes[i],label = 'F Archimedes '+dict_process['time_range'][i])

    plt.xlabel('time (s)')
    plt.ylabel('F archimedes (N)')
    plt.title ('Comparison of different fitted archimedes forces')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(OUTPUT_D+'Comparison of different fitted archimedes forces.png',bbox_inches='tight')
    plt.clf()


if __name__ == '__main__':
    d_input = read_output(PATH_IN)
    d_process = process_output(d_input,PATH_OUT)
    plot_deltaV(d_process)