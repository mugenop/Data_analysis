import numpy as np
from matplotlib import pyplot as plt
import os


# ##5
Yg_reference = [2.57520e-02,2.58761e-02,2.59152e-02,2.57457e-02,2.57457e-02]
ly_reference = [5.10853e-02,5.10960e-02,5.12694e-02,5.10834e-02,5.10834e-02]
ext_V_reference =[6.54484e-05,6.54422e-05,6.54462e-05,6.54477e-05,6.54480e-05]
timeStartbuckling = [1500,1000,163,1821,435]
timeProgressionRolling = [False,False,True,True,True]
IndexPressure = 340
timStartArraybuckling = np.array(timeStartbuckling)*0.0002
timeStopArraybuckling = timStartArraybuckling+0.1
timeStartdebuckling = [5850,5850,180,161,1665]
timStartArraydebuckling = np.array(timeStartdebuckling)*0.0002
timeStopArraydebuckling = timStartArraydebuckling+0.1
pressure_ref = 101300
#



# # #6.5
# Yg_reference = [2.63940e-02,2.75485e-02,2.76895e-02,2.72654e-02,2.74623e-02]
# ly_reference = [5.14205e-02,5.25599e-02,5.27431e-02,5.23628e-02,5.25493e-02]
# ext_V_reference =[6.54476e-05,6.54452e-05,6.54495e-05,6.54458e-05,6.54447e-05]
# timeStartbuckling = [1530,1650,1435,975,1718]
# timeProgressionRolling = [False,False,True,True,True]
# # IndexPressure = 340
# timStartArraybuckling = np.array(timeStartbuckling)*0.0002
# timeStopArraybuckling = timStartArraybuckling+0.1
# timeStartdebuckling = [2765,3985,4240,3925,2830]
# timStartArraydebuckling = np.array(timeStartdebuckling)*0.0002
# timeStopArraydebuckling = timStartArraydebuckling+0.1
# pressure_ref = 101300
#
#
# #2
# Yg_reference = [2.62422e-02,2.62123e-02,2.64830e-02,2.65178e-02,2.65126e-02]
# ly_reference = [5.24448e-02,5.22848e-02,5.27078e-02,5.27705e-02,5.28126e-02]
# ext_V_reference =[6.54160e-05,6.54171e-05,6.54485e-05,6.52726e-05,6.54427e-05]
# timeStartbuckling = [1530,1550,922,1000,30]
# timeProgressionRolling = [False,False,True,True,True]
# timStartArraybuckling = np.array(timeStartbuckling)*0.0002
# timeStopArraybuckling = timStartArraybuckling+0.1
# timeStartdebuckling = [4255,4280,4480,2330,1370]
# timStartArraydebuckling = np.array(timeStartdebuckling)*0.0002
# timeStopArraydebuckling = timStartArraydebuckling+0.1
# pressure_ref = 101300

def file_reader(name,type):
    f = open(name,"r")
    lines = f.readlines()
    print name
    if type=="results":
        timeList = []
        pressureList = []
        deltaPList = []
        lowPointList = []
        gravity_yList = []
        externalVList = []
        elements = lines[0].split('\t')
        image_name = False
        if elements[2]=="image_name":
            image_name =True
        if image_name:
            for e in lines[1:]:
                if len(e)<10:
                    pass
                else:
                    elements = e.split('\t')
                    center_of_gravity =float(elements[6])
                    volume = float(elements[7])
                    ylow = float(elements[4])
                    test01 = center_of_gravity>=0
                    #if not test01:
                        #print "center of gravity negative"
                    test02 = center_of_gravity< 0.1
                    #if not test02:
                        #print "center of gravity higher than thresh"
                    test03 = volume>1e-05
                    #if not test03:
                        #print "volume lower than thresh"
                    test04 = volume<10e-05
                    #if not test04:
                        #print "volume higher than thresh"
                    test05 = ylow>0
                    #if not test05:
                        #print "ylow negative"
                    test06 = ylow<0.1
                    #if not test06:
                        #print "ylow higher than thresh"
                    test = test01 and test02 and test03 and test04 and test05 and test06
                    if test:
                        timeList.append(float(elements[0]))
                        pressureList.append(float(elements[1]))
                        deltaPList.append(float(3))
                        lowPointList.append(float(elements[4]))
                        gravity_yList.append(float(elements[6]))
                        externalVList.append(float(elements[7]))
        else:
            for e in lines[1:]:
                if len(e)<10:
                    pass
                else:
                    elements = e.split('\t')
                    center_of_gravity =float(elements[5])
                    volume = float(elements[6])
                    ylow = float(elements[3])
                    test01 = center_of_gravity>=0
                    #if not test01:
                        #print "center of gravity negative"
                    test02 = center_of_gravity< 0.1
                    #if not test02:
                        #print "center of gravity higher than thresh"
                    test03 = volume>1e-05
                    #if not test03:
                        #print "volume lower than thresh"
                    test04 = volume<10e-05
                    #if not test04:
                        #print "volume higher than thresh"
                    test05 = ylow>0
                    #if not test05:
                        #print "ylow negative"
                    test06 = ylow<0.1
                    #if not test06:
                        #print "ylow higher than thresh"
                    test = test01 and test02 and test03 and test04 and test05 and test06
                    if test:
                        timeList.append(float(elements[0]))
                        pressureList.append(float(elements[1]))
                        deltaPList.append(float(2))
                        lowPointList.append(float(elements[3]))
                        gravity_yList.append(float(elements[5]))
                        externalVList.append(float(elements[6]))
        return np.array(timeList),np.array(pressureList),np.array(deltaPList),np.array(lowPointList),np.array(gravity_yList),np.array(externalVList)

def plot_multiple(datax,datay,numberOfexperiments,outputFilename,title,xlabel,ylabel,show_flag=False,isotropic= False):
    if not isotropic:
        color_codes=["mo","co","ro","bo","go"]
        fig = plt.figure(figsize=(17,17), dpi=100)
        fig.autolayout = True
        for i in xrange(len(datax)):
            labelname = "Manip"+str(i+1)
            plt.plot(datax[i],datay[i],color_codes[i],label=labelname,markersize = 10)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(title)
        lgd = plt.legend(bbox_to_anchor=(1, 1), loc=2, borderpad=1.,fancybox=True,markerscale=1,fontsize=18)
        plt.grid()
        plt.savefig(outputFilename, bbox_inches='tight')
        if show_flag:
            plt.show()
        plt.close()
    else:
        color_codes = ["m","c","r","b","g"]
        sign_codes=["o","*","^","s","d"]
        experiment_parts = ["iso_b","iso_db","Equi_b","Equi_db","rolling"]
        fig = plt.figure(figsize=(17, 17), dpi=200)
        fig.autolayout = True
        nbr_of_parts = numberOfexperiments[0]
        for i in xrange(len(datax)):
            if i>=2*numberOfexperiments[0]:
                nbr_of_parts = numberOfexperiments[1]
                manip_nbr = 3+(i-(2*numberOfexperiments[0]))/nbr_of_parts
                part_nbr = (i-(2*numberOfexperiments[0]))%nbr_of_parts
            else:
                manip_nbr = 1+(i/nbr_of_parts)
                part_nbr = i%nbr_of_parts
            labelname = str(manip_nbr)+"_"+experiment_parts[part_nbr]
            try:
                plt.plot(datax[i],datay[i],color_codes[manip_nbr-1]+sign_codes[part_nbr],label=labelname,markersize = 10)

            except:
                print 'here'
                print i
                print nbr_of_parts
                print manip_nbr
                print part_nbr
                print len(datax)
                print len(datay)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(title)
        lgd = plt.legend(bbox_to_anchor=(1, 1), loc=2, borderpad=1.,fancybox=True,markerscale=1,fontsize=18)
        plt.grid()
        plt.savefig(outputFilename, bbox_inches='tight')
        if show_flag:
            plt.show()
        plt.close()

def plot_unique(datax,datay,outputFilename,graphlabel,graphtype,title,xlabel,ylabel,show_flag=False):
    # try:
    #     datax_min = np.amin(datax)
    # except:
    #     print datax
    #     print datay
    #     print outputFilename+" problem"
    # datax_max = np.amax(datax)
    # step_x = (datax_max-datax_min)/float(datax.shape[0])
    # datay_min = np.amin(datay)
    #
    # datay_max = np.amax(datay)


    # step_y = (datay_max-datay_min)/float(datay.shape[0])
    # print outputFilename
    # print datax_min-step_x,datax_max+step_x
    # print datay_min-step_y,datay_max+step_y
    #
    fig = plt.figure(figsize=(17, 17), dpi=200)
    fig.autolayout = True

    # plt.xlim(datax_min-step_x,datax_max+step_x)
    # plt.ylim(datay_min-step_y,datay_max+step_y)
    plt.plot(datax,datay,graphtype,label=graphlabel,markersize = 10)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    lgd = plt.legend(bbox_to_anchor=(1, 1), loc=2, borderpad=1.,fancybox=True,markerscale=1,fontsize=18)
    plt.grid()
    plt.savefig(outputFilename, bbox_inches='tight')
    if show_flag:
        plt.show()
    plt.close()

def plot_double(datax,datay1,datay2,graphlabel1,graphlabel2,outputFilename,title,xlabel,ylabel,show_flag=False, ylim=None,xlim=None):
    fig = plt.figure(figsize=(17, 17), dpi=200)
    fig.autolayout = True
    plt.plot(datax,datay1,'r^',label=graphlabel1,markersize = 10)
    plt.plot(datax,datay2,'b-',label=graphlabel2,markersize = 10)
    if ylim is not None:
            plt.ylim(ylim)
    if xlim is not None:
            plt.xlim(xlim)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    lgd = plt.legend(bbox_to_anchor=(1, 1), loc=2, borderpad=1.,fancybox=True,markerscale=1,fontsize=18)
    plt.grid()
    plt.savefig(outputFilename, bbox_inches='tight')
    if show_flag:
        plt.show()
    plt.close()

