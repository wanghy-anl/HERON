
# Copyright 2020, Battelle Energy Alliance, LLC
# ALL RIGHTS RESERVED
"""
  Example class for validators.
"""
import numpy as np
# import pickle as pk
# import os
import sys
import xml.etree.ElementTree as ET

from utils import InputData, InputTypes

from .Validator import Validator

class Para_RefGov_SESBOPTES_MW(Validator):
  """
    A Reference Governor Validator for dispatch decisions.
    Accepts parameterized A,B,C,D matrices from external XML file, and validate 
    the dispatch power (BOP & SES, unit=MW), and 
    the next stored energy level (TES, unit=MWh)

    Haoyu Wang, ANL-NSE, Jan 21, 2021
  """
  # ---------------------------------------------
  # INITIALIZATION
  @classmethod
  def get_input_specs(cls):
    print("Haoyu TEST: get_input_specs\n")
    """
      Set acceptable input specifications.
      @ In, None
      @ Out, specs, InputData, specs
    """
    specs = Validator.get_input_specs()
    specs.name = 'Para_RefGov_SESBOPTES_MW'
    specs.description = r"""Uses a demonstration-only validator that constrains the change
          for any resource to a constant ``delta''."""
    specs.addSub(InputData.parameterInputFactory('delta', contentType=InputTypes.FloatType,
        descr=r"""the maximum absolute change in any resource between successive time steps."""))
    specs.addSub(InputData.parameterInputFactory('tolerance', contentType=InputTypes.FloatType,
        descr=r"""the strictness with which the constraint should be enforced. Note that some small
              numerical exception is expected."""))
    return specs

  def __init__(self):
    print("Haoyu TEST: __init__\n")
    """
      Constructor.
      @ In, None
      @ Out, None
    """
    self.name = 'BaseValidator'
    self._tolerance = 1.003e-6
    self._unitInfo = { # a two-layer dictionary, containing the matrix, lower/upper constraints, and XLast to start with
      'BOP':{
        'MatrixFile':"D:\\GitProjects\\wanghy_fork\\HERON\\tests\\integration_tests\\validator_FARM\\BOP_SES_TES_Mat_Unit_MW\\DMDcCxCoeff_BOP_para.xml",
        'Targets_Min':[7.25E8, 51.8E5], # Target Variable minimum values, BOP outputs [BOP.portElec_b.W, BOP_steamTurbine_portHP_p], units [W,Pascal]
        'Targets_Max':[14.5E8, 53.8E5], # Target Variable maximum values, BOP outputs [BOP.portElec_b.W, BOP_steamTurbine_portHP_p], units [W,Pascal]
        'XLast':[], # the final system state vector
        'vp_hist':[], # the history of v, absolute value, MW
        'y_hist':[] # the history of y, absolute value
        },
      'SES':{
        'MatrixFile':"D:\\GitProjects\\wanghy_fork\\HERON\\tests\\integration_tests\\validator_FARM\\BOP_SES_TES_Mat_Unit_MW\\DMDcCxCoeff_SES_para.xml",
        'Targets_Min':[13.5E6, 1089.], # Target Variable minimum values, SES outputs [SES.portElec_b.W, SES.GTunit.Tf], units [W,K]
        'Targets_Max':[50.0E6, 1873.], # Target Variable maximum values, SES outputs [SES.portElec_b.W, SES.GTunit.Tf], units [W,K]
        'XLast':[], # the final system state vector
        'vp_hist':[], # the history of v, absolute value, MW
        'y_hist':[] # the history of y, absolute value
        },  
      'TES':{
        'Initial_Level': 5.0E5, # initial stored energy (MWh), to be erased when new API becomes available
        'MatrixFile':"D:\\GitProjects\\wanghy_fork\\HERON\\tests\\integration_tests\\validator_FARM\\BOP_SES_TES_Mat_Unit_MW\\DMDcCxCoeff_TES_para.xml",
        'Targets_Min':[2.5, 2.5], # Target Variable minimum values
        'Targets_Max':[55., 55.], # Target Variable maximum values
        'XLast':[], # the final system state vector
        'vp_hist':[], # the history of v, absolute value, MW
        'y_hist':[] # the history of y, absolute value
        }       
      }

  def read_input(self, inputs):
    print("Haoyu TEST: read_input\n")
    """
      Loads settings based on provided inputs
      @ In, inputs, InputData.InputSpecs, input specifications
      @ Out, None
    """
    pass

  # ---------------------------------------------
  # API
  def validate(self, components, dispatch, times, meta):
    print("Haoyu TEST: validate\n")
    """
      Method to validate a dispatch activity.
      @ In, components, list, HERON components whose cashflows should be evaluated
      @ In, activity, DispatchState instance, activity by component/resources/time
      @ In, times, np.array(float), time values to evaluate; may be length 1 or longer
      @ In, meta, dict, extra information pertaining to validation
      @ Out, errs, list, information about validation failures
    """
    # errs will be returned to dispatcher. errs contains all the validation errors calculated in below
    errs = [] # TODO best format for this?

    # get time interval
    Tr_Update_hrs = float(times[1]-times[0])
    Tr_Update_sec = Tr_Update_hrs*3600.

    # loop through the <Component> items in HERON
    for comp, info in dispatch._resources.items(): 
      print("Haoyu Debug, comp=",comp, str(comp)) # e.g. comp= <HERON Component "SES""> <HERON Component "SES"">
      # loop through the items defined in the __init__ function
      for unit in self._unitInfo: 
        print("Haoyu Debug, CompInfo, unit=",unit) # e.g. CompInfo, unit= SES
        # Identify the profile as defined in the __init__ function
        if str(unit) not in str(comp):
          # If the "unit" and "comp" do not match, go to the next "unit" in loop
          continue
        else: # when the str(unit) is in the str(comp) (e.g. "SES" in "<HERON Component "SES"">")
          self._unitInfo[unit]['vp_hist']=[]; self._unitInfo[unit]['y_hist']=[]
          """ Read State Space XML file (generated by Raven parameterized DMDc) """
          MatrixFile = self._unitInfo[unit]['MatrixFile']
          Tss, n, m, p, para_array, UNorm_list, XNorm_list, XLast_list, YNorm_list, A_list, B_list, C_list, eig_A_array = read_parameterized_XML(MatrixFile)
          
          """ MOAS steps Limit """
          g = int(Tr_Update_sec/Tss)+1 # numbers of steps to look forward, , type = <class 'int'>
          # print("Haoyu Debug, g= ",g)
          """ Keep only the profiles with YNorm within the [y_min, y_max] range """
          y_min = self._unitInfo[unit]['Targets_Min']
          y_max = self._unitInfo[unit]['Targets_Max']     
          para_array, UNorm_list, XNorm_list, XLast_list, YNorm_list, A_list, B_list, C_list, eig_A_array = check_YNorm_within_Range(
            y_min, y_max, para_array, UNorm_list, XNorm_list, XLast_list, YNorm_list, A_list, B_list, C_list, eig_A_array)
          if YNorm_list == []:
            sys.exit('ERROR:  No proper linearization point (YNorm) found in Matrix File. \n\tPlease provide a state space profile linearized within the [y_min, y_max] range\n')
          max_eigA_id = eig_A_array.argmax()
          A_m = A_list[max_eigA_id]; B_m = B_list[max_eigA_id]; C_m = C_list[max_eigA_id]; D_m = np.zeros((p,m)) # all zero D matrix
          
          # loop through the resources in info (only one resource here - electricity)
          for res in info:
            if str(res) == "electricity":
              # loop through the time index (tidx) and time in "times"
              if self._unitInfo[unit]['XLast']==[]:
                x_sys = np.zeros(n)
              
              for tidx, time in enumerate(times):
                # Copy the system state variable
                x_KF = x_sys
                """ Get the r_value, original actuation value """
                current = float(dispatch.get_activity(comp, res, times[tidx]))
                # check if TES: power = (curr. MWh energy - prev. MWh energy)/interval Hrs
                
                if comp.get_interaction().is_type('Storage') and tidx == 0:
                  init_level = comp.get_interaction().get_initial_level(meta)
                  
                if str(unit) == "TES":
                  # Initial_Level = float(self._unitInfo[unit]['Initial_Level'])
                  Initial_Level = float(init_level)
                  if tidx == 0: # for the first hour, use the initial level. charging yields to negative r_value
                    r_value = -(current - Initial_Level)/Tr_Update_hrs
                  else: # for the other hours
                    # r_value = -(current - float(dispatch.get_activity(comp, res, times[tidx-1])))/Tr_Update_hrs
                    r_value = -(current - Allowed_Level)/Tr_Update_hrs
                else: # when not TES, 
                  r_value = current # measured in MW
                
                """ Find the correct profile according to r_value"""
                profile_id = (np.abs(para_array - r_value)).argmin()
                
                # Retrive the correct A, B, C matrices
                A_d = A_list[profile_id]; B_d = B_list[profile_id]; C_d = C_list[profile_id]; D_d = np.zeros((p,m)) # all zero D matrix
                # Retrive the correct y_0, r_0 and X
                y_0 = YNorm_list[profile_id]; r_0 = float(UNorm_list[profile_id])
                
                # Build the s, H and h for MOAS
                s = [] # type == <class 'list'>
                for i in range(0,p):
                  s.append([abs(y_max[i] - y_0[i])])
                  s.append([abs(y_0[i] - y_min[i])])
                # print(s)
                H, h = fun_MOAS_noinf(A_d, B_d, C_d, D_d, s, g) # H and h, type = <class 'numpy.ndarray'>
                # print("H", H); print("h", h)
                # first v_RG: consider the step "0" - step "g"
                v_RG = fun_RG_SISO(0, x_KF, r_value-r_0, H, h, p) # v_RG: type == <class 'numpy.ndarray'>
                
                """ 2nd adjustment """
                # MOAS for the steps "g+1" - step "2g"
                Hm, hm = fun_MOAS_noinf(A_m, B_m, C_m, D_m, s, g)
                # Calculate the max/min for v, ensuring the hm-Hxm*x(g+1) always positive for the next g steps.
                v_max, v_min = fun_2nd_gstep_calc(x_KF, Hm, hm, A_m, B_m, g)

                if v_RG < v_min:
                  v_RG = v_min
                elif v_RG > v_max:
                  v_RG = v_max
                
                v_value = v_RG + r_0 # absolute value of electrical power (MW)
                v_value = float(v_value)

                # Update x_sys, and keep record in vp_hist and yp_hist within this hour
                for i in range(int(Tr_Update_sec/Tss)):
                  self._unitInfo[unit]['vp_hist'].append(v_value)
                  y_sim = np.dot(C_d,x_sys)
                  self._unitInfo[unit]['y_hist'].append(y_sim+y_0)
                  x_sys = np.dot(A_d,x_sys)+np.dot(B_d,v_RG)

                # Convert to V1:
                                
                if str(unit) == "TES":
                  if tidx == 0: # for the first hour, use the initial level
                    Allowed_Level = Initial_Level - v_value*Tr_Update_hrs # Allowed_Level: predicted level due to v_value
                  else: # for the other hours
                    Allowed_Level = Allowed_Level - v_value*Tr_Update_hrs
                  V1 = Allowed_Level
                else: # when not TES, 
                  V1 = v_value
                
                # print("Haoyu Debug, unit=",str(unit),", t=",time, ", curr=",current, ", V1= ",V1, ", r(MW)=",r_value, ", vp(MW)=", v_value)
                # print("Haoyu Debug, unit=",str(unit),", t=",time, ", curr=",current, ", V1= ",V1, ", delta=", (V1-current))
                print("Haoyu Debug, unit=",str(unit),", t=",time, ", curr= %.8g, V1= %.8g, delta=%.8g" %(current, V1, (V1-current)))

                # Write up any violation to the errs:
                if abs(current - V1) > self._tolerance*max(abs(current),abs(V1)):
                  # violation
                  errs.append({'msg': f'Reference Governor Violation',
                              'limit': V1,
                              'limit_type': 'lower' if (current < V1) else 'upper',
                              'component': comp,
                              'resource': res,
                              'time': time,
                              'time_index': tidx,
                              })

    if errs == []: # if no validation error:
      print(" ")
      print("*********************************************************************")
      print("*** Haoyu Debug, Validation Success, Print for offline processing ***")
      print("*********************************************************************")
      print(" ")
      t_hist = np.arange(0,len(self._unitInfo['TES']['vp_hist'])*Tss,Tss)
      for unit in self._unitInfo:
        y_hist = np.array(self._unitInfo[unit]['y_hist']).T
        # print(str(unit),y_hist)
        for i in range(len(t_hist)):
          print(str(unit), ",t,",t_hist[i],",vp,",self._unitInfo[unit]['vp_hist'][i],",y1,",y_hist[0][i], ",y1min,",self._unitInfo[unit]['Targets_Min'][0],",y1max,",self._unitInfo[unit]['Targets_Max'][0],",y2,",y_hist[1][i], ",y2min,",self._unitInfo[unit]['Targets_Min'][1],",y2max,",self._unitInfo[unit]['Targets_Max'][1])



          # for res in info:
          #   print("Haoyu Debug, info - res = ", str(res), type(str(res))) # info - res =  electricity
          #   for t, time in enumerate(times):
          #     current = dispatch.get_activity(comp, res, time) # "current is measured in MW
          #     print("Haoyu Debug, t=", times[t], "t=", time, "current activity=",current)
          #     # r_value_RG = current*1E6 - r_0[0]
          #     if str(unit) == "BOP":
          #       v_max = 1130.838; v_min = 1083.233 # units: MW
          #     elif str(unit) == "SES":
          #       v_max = 50.000; v_min = 19.317 # units: MW
          #     elif str(unit) == "TES":
          #       v_max = 5.1E5; v_min = 4.9E5

          #     if current > v_max:
          #       V1 = v_max
          #     elif current < v_min:
          #       V1 = v_min
          #     else:
          #       V1 = current
              
              
          #     if abs(current - V1) > self._tolerance:
          #       # violation
          #       errs.append({'msg': f'Reference Governor Violation',
          #                     'limit': V1,
          #                     'limit_type': 'lower' if (current < V1) else 'upper',
          #                     'component': comp,
          #                     'resource': res,
          #                     'time': time,
          #                     'time_index': t,
          #                     })

    return errs


class RefGov_SESBOP_W(Validator):
  """
    Example class for validating dispatch decisions.
    Arbitrarily, uses percent based ramping limits
  """
  # ---------------------------------------------
  # INITIALIZATION
  @classmethod
  def get_input_specs(cls):
    print("Haoyu TEST: get_input_specs\n")
    """
      Set acceptable input specifications.
      @ In, None
      @ Out, specs, InputData, specs
    """
    specs = Validator.get_input_specs()
    specs.name = 'RefGov_SESBOP_W'
    # TODO left for convenience
    return specs

  def __init__(self):
    print("Haoyu TEST: __init__\n")
    
    """
      Constructor.
      @ In, None
      @ Out, None
    """
    self.name = 'BaseValidator'
    # self._allowable = 2
    self._tolerance = 1e-10
    self._CompInfo = { # a two-layer dictionary, containing the matrix, lower/upper constraints, and MOAS steps for each component
      'SES':{
        'MatrixFile':"D:\\GitProjects\\wanghy_fork\\HERON\\tests\\integration_tests\\validator_FARM\\BOP_SES_Matrices_unit_W\\DMDcCxCoeff_SES.xml",
        'Targets_Min':[12.0E6, 1089.], # Target Variable minimum values
        'Targets_Max':[40.0E6, 1873.], # Target Variable maximum values
        'MOASsteps': 10, # the projection steps for RefGov
        'XLast':[] # the final system state vector
        },
      'BOP':{
        'MatrixFile':"D:\\GitProjects\\wanghy_fork\\HERON\\tests\\integration_tests\\validator_FARM\\BOP_SES_Matrices_unit_W\\DMDcCxCoeff_BOP.xml",
        'Targets_Min':[7.5E8, 65.8E5], # Target Variable minimum values
        'Targets_Max':[1.5E9, 67.8E5], # Target Variable maximum values
        'MOASsteps': 10, # the projection steps for RefGov
        'XLast':[] # the final system state vector
        }
      }
    
    # self._serialized_file = "D:\\GitProjects\\wanghy_fork\\HERON\\tests\\integration_tests\\validator_RefGov\\RefGov_SES.pk" # this must be taken from the HERON input file (read_input)
    # self._component = "SES" # name of the component must be associated to this validator (read_input)
    # self._prod_variable  = "P1"
    # obj = open(self._serialized_file,"rb+")
    # print(type(obj))
    # self.ref_gov = pk.load(open(self._serialized_file,"rb+"))
    # self.ref_gov = pk.load(self._serialized_file)

  def read_input(self, inputs):
    print("Haoyu TEST: read_input\n")
    """
      Loads settings based on provided inputs
      @ In, inputs, InputData.InputSpecs, input specifications
      @ Out, None
    """
    pass

  # ---------------------------------------------
  # API
  def validate(self, components, dispatch, times, meta):
    print("Haoyu TEST: validate\n")
    """
      Method to validate a dispatch activity.
      @ In, components, list, HERON components whose cashflows should be evaluated
      @ In, activity, DispatchState instance, activity by component/resources/time
      @ In, times, np.array(float), time values to evaluate; may be length 1 or longer
      @ In, meta, dict, extra information pertaining to validation
      @ Out, errs, list, information about validation failures
    """
    # self.ref_gov = pk.load(open(self._serialized_file,"rb+"))

    errs = [] # TODO best format for this?
    for comp, info in dispatch._resources.items():
      print("Haoyu Debug, comp=",comp, str(comp)) # comp= <HERON Component "SES""> <HERON Component "SES"">
      for x in self._CompInfo:
        print("Haoyu Debug, CompInfo, x=",x) # Haoyu Debug, CompInfo, x= SES / Haoyu Debug, CompInfo, x= BOP
        # print(self._CompInfo[x])
        if str(x) not in str(comp):
          # If the x and comp do not match, go to the next in loop
          # print('NOT MATCH')
          continue
        else:
          # If the x and comp match, perform RG steps:
          # print('MATCH!!!')
          # Step 1: Load XML file containing matrices
          # UNorm, XNorm, XLast, YNorm, Atilde, Btilde, Ctilde = 
          r_0,     XNorm, XLast, y_0,   A_d,    B_d,    C_d    = read_XML_input(self._CompInfo[x]['MatrixFile'])
          # print("r_0=",r_0)
          n = len(A_d); m = len(B_d[0]); p = len(C_d)  # n: dim of x; m: dim of u; p: dim of y. Type = <class 'int'>
          D_d = np.zeros((p,m)) # all zero D matrix
          
          # for res in info:
          #   for t, time in enumerate(times):
          #     current = dispatch.get_activity(comp, res, time) # "current is measured in MW
          #     # print("Haoyu Debug, current=",current)
          #     # r_value_RG = current*1E6 - r_0[0]
          #     if str(x) == "BOP":
          #       v_max = 1450; v_min = 1400 # units: MW
          #     elif str(x) == "SES":
          #       v_max = 40; v_min = 15 # units: MW
              
          #     if current > v_max:
          #       V1 = v_max
          #     elif current < v_min:
          #       V1 = v_min
          #     else:
          #       V1 = current
              


            

          # Step 2: get the y_min, y_max, and g
          y_min = self._CompInfo[x]['Targets_Min'] # [20.0E6, 1270.]  # minimum values
          y_max = self._CompInfo[x]['Targets_Max'] # [50.0E6, 1410.]  # maximum value

          g = int(self._CompInfo[x]['MOASsteps'])  # numbers of steps to look forward

          # Step 3: Calculate Maximal Output Admissible Set (MOAS) 
          s = [] # type == <class 'list'>
          for i in range(0, p):
            s.append([y_max[i] - y_0[i]])
            s.append([y_0[i] - y_min[i]])
            # print(s)
          H, h = fun_MOAS(A_d, B_d, C_d, D_d, s, g) # H and h, type = <class 'numpy.ndarray'>
          # print("H:\n", H); print("h:\n", h)

          # Step 4: perform RG
      
          for res in info:
            # print("Haoyu Debug, res=",res)
            for t, time in enumerate(times):
              print("Haoyu Debug, res=", res, 't=', t, 'time=', time) # Haoyu Debug, res= electricity t= 0 time= 0.0 / t= 20 time= 2.0
              current = dispatch.get_activity(comp, res, time) # "current is measured in MW
              # print("Haoyu Debug, current=",current)
              r_value_RG = current*1E6 - r_0[0]

              # Check the XLast: if not updated, use the initial one from XML file
              if self._CompInfo[x]['XLast']==[]:
                X_Last_RG = XLast - XNorm
              else:
                X_Last_RG = np.asarray(self._CompInfo[x]['XLast'])
              # print("X_Last_RG",type(X_Last_RG),np.shape(X_Last_RG),X_Last_RG)
              

              """ Call the Reference Governor to mild the r_value """
              v_RG, v_min, v_max = fun_RG_SISO_vBound(0, X_Last_RG, r_value_RG, H, h, p)  # v_RG: type == <class 'numpy.ndarray'>
              V1 = (v_RG + r_0[0])*1E-6 # V1 is measured in MW
              
              # Print out y[2], the hidden constraints got violated
              # print("comp={},t={}\n".format(str(x),t))
              HxxHvv=np.dot(H[:, 0:n],X_Last_RG)+np.dot(H[:, n],v_RG)
              # print("comp={}, y_0={}, t={}, v={}, Hx*x+Hv*v_RG=\n".format(str(x),y_0,t,(v_RG+r_0)),HxxHvv)
              # print(max(HxxHvv[0:4*p:]))
              Y1Y2 = np.asarray([max(HxxHvv[0::2*p]),max(HxxHvv[2::2*p])])+y_0
              # print(Y1Y2)
              # print("comp={}, t={}, V1={}, y_out={}, X_Last={}\n".format(str(x),t,V1,Y1Y2,X_Last_RG))
              # print("Hv*v_RG=\n",np.dot(H[:, n],v_RG))
              # print("y_0=\n",y_0)

              # Update the XLast
              # print(np.shape(np.dot(A_d,X_Last_RG)))
              # print(np.shape(np.dot(B_d,np.reshape(v_RG,-1))))
              self._CompInfo[x]['XLast']=np.dot(A_d,X_Last_RG)+np.dot(B_d,np.reshape(v_RG,-1))



              
              if abs(current - V1) > self._tolerance:
                # violation
                errs.append({'msg': f'Reference Governor Violation',
                              'limit': V1,
                              'limit_type': 'lower' if (current < V1) else 'upper',
                              'component': comp,
                              'resource': res,
                              'time': time,
                              'time_index': t,
                              })

    return errs

def read_XML_input(MatrixFileName):
  tree = ET.parse(MatrixFileName)
  root = tree.getroot()
  for child1 in root:
    # print(' ',child1.tag) # DMDrom
    for child2 in child1:
      # print('  > ', child2.tag) # ROM, DMDcModel

      for child3 in child2:
        # print('  >  > ', child3.tag) # UNorm, XNorm, XLast, Atilde, Btilde, Ctilde
        if child3.tag == 'UNorm':
          # print(child3.text)
          Temp_txtlist = child3.text.split('; ')
          Temp_floatlist = [float(item) for item in Temp_txtlist]
          UNorm = np.asarray(Temp_floatlist)
          # print(np.shape(self.UNorm))
        if child3.tag == 'XNorm':
          # print(child3.text)
          Temp_txtlist = child3.text.split('; ')
          Temp_floatlist = [float(item) for item in Temp_txtlist]
          XNorm = np.asarray(Temp_floatlist)
          # print(np.shape(self.XNorm))
        if child3.tag == 'XLast':
          # print(child3.text)
          Temp_txtlist = child3.text.split('; ')
          Temp_floatlist = [float(item) for item in Temp_txtlist]
          XLast = np.asarray(Temp_floatlist)
          # print(np.shape(self.XLast))
        if child3.tag == 'YNorm':
          # print(child3.text)
          Temp_txtlist = child3.text.split('; ')
          Temp_floatlist = [float(item) for item in Temp_txtlist]
          YNorm = np.asarray(Temp_floatlist)
          # print(np.shape(self.YNorm))

        for child4 in child3:
          # print('  >  >  > ', child4.tag) # real, imaginary, matrixShape, formatNote
          if child4.tag == 'real':
            Temp_txtlist = child4.text.split('; ')
            Temp_txtlist = [item.split(' ') for item in Temp_txtlist]
            Temp_floatlist = [[float(y) for y in x ] for x in Temp_txtlist]
            # print(Temp_txtlist)
            # print(Temp_floatlist)
            if child3.tag == 'Atilde':
              A_Re = np.asarray(Temp_floatlist)
            if child3.tag == 'Btilde':
              B_Re = np.asarray(Temp_floatlist)
            if child3.tag == 'Ctilde':
              C_Re = np.asarray(Temp_floatlist)

          if child4.tag == 'imaginary':
            Temp_txtlist = child4.text.split('; ')
            Temp_txtlist = [item.split(' ') for item in Temp_txtlist]
            Temp_floatlist = [[float(y) for y in x ] for x in Temp_txtlist]
            # print(Temp_txtlist)
            # print(Temp_floatlist)
            if child3.tag == 'Atilde':
              A_Im = np.asarray(Temp_floatlist)
            if child3.tag == 'Btilde':
              B_Im = np.asarray(Temp_floatlist)
            if child3.tag == 'Ctilde':
              C_Im = np.asarray(Temp_floatlist)

  Atilde = A_Re
  Btilde = B_Re
  Ctilde = C_Re
  return UNorm, XNorm, XLast, YNorm, Atilde, Btilde, Ctilde

def fun_MOAS(A, B, C, D, s, g):
  p = len(C)  # dimension of y
  T = np.linalg.solve(np.identity(len(A)) - A, B)
  """ Build the S matrix"""
  S = np.zeros((2 * p, p))
  for i in range(0, p):
    S[2 * i, i] = 1.0
    S[2 * i + 1, i] = -1.0
  Kx = np.dot(S, C)
  # print("Kx", Kx)
  Lim = np.dot(S, (np.dot(C, T) + D))
  # print("Lim", Lim)
  Kr = np.dot(S, D)
  # print("Kr", Kr)
  """ Build the core of H and h """
  H = np.concatenate((0 * Kx, Lim), axis=1)
  h = s
  NewBlock = np.concatenate((Kx, Kr), axis=1)
  H = np.concatenate((H, NewBlock))
  h = np.concatenate((h, s))

  """ Build the add-on blocks of H and h """
  i = 0
  while i < g:
    i = i + 1
    Kx = np.dot(Kx, A)
    Kr = Lim - np.dot(Kx, T)

    NewBlock = np.concatenate((Kx, Kr), axis=1)
    H = np.concatenate((H, NewBlock))
    h = np.concatenate((h, s))
    """ To Insert the ConstRedunCheck """

  return H, h

def fun_RG_SISO_vBound(v_0, x, r, H, h, p):
  n = len(x)  # dimension of x
  x = np.vstack(x)  # x is horizontal array, must convert to vertical for matrix operation
  # because v_0 and r are both scalar, so no need to vstack
  Hx = H[:, 0:n];  Hv = H[:, n:]
  alpha = h - np.dot(Hx, x) - np.dot(Hv, v_0)  # alpha is the system remaining vector
  beta = np.dot(Hv, (r - v_0))  # beta is the anticipated response vector with r
  # print("alpha = \n",alpha)
  # print("Hv = \n", Hv)
  # print("beta = \n", beta)

  """ Calculate the vBounds """
  v_st = [] # smaller than
  v_bt = [] # bigger than
  # for the first 2p rows (final steady state corresponding to constant v), keep the max/min.
  for k in range(0, 2 * p):
    if Hv[k] > 0:
      v_st.append(alpha[k] / Hv[k] + v_0)
    elif Hv[k] < 0:
      v_bt.append(alpha[k] / Hv[k] + v_0)
  # for the following rows, adjust the max/min when necessary.
  for k in range(2 * p, len(alpha)):
    if Hv[k] > 0 and alpha[k] > 0 and alpha[k] < beta[k]:
      v_st.append(alpha[k] / Hv[k] + v_0)
    elif Hv[k] < 0 and alpha[k] > 0 and alpha[k] < beta[k]:
      v_bt.append(alpha[k] / Hv[k] + v_0)
  v_max = float(min(v_st))
  v_min = float(max(v_bt))
  # print("r =",r,"type=",type(r))
  # print("v_min=",v_min,"type=",type(v_min))
  # print("v_max=",v_max,"type=",type(v_max))

  if r > v_max:
    v = v_max
  elif r < v_min:
    v = v_min
  else:
    v = r

  # v = np.asarray(v).flatten()
  # print("v=",v,"type=",type(v))

  return v, v_min, v_max

def read_parameterized_XML(MatrixFileName):
  tree = ET.parse(MatrixFileName)
  root = tree.getroot()
  para_array = []; UNorm_list = []; XNorm_list = []; XLast_list = []; YNorm_list =[]
  A_Re_list = []; B_Re_list = []; C_Re_list = []; A_Im_list = []; B_Im_list = []; C_Im_list = []  
  for child1 in root:
    # print(' ',child1.tag) # DMDrom
    for child2 in child1:
      # print('  > ', child2.tag) # ROM, DMDcModel
      for child3 in child2:
        # print('  >  > ', child3.tag) # dmdTimeScale, UNorm, XNorm, XLast, Atilde, Btilde, Ctilde
        if child3.tag == 'dmdTimeScale':
          # print(child3.text)
          Temp_txtlist = child3.text.split(' ')
          Temp_floatlist = [float(item) for item in Temp_txtlist]
          TimeScale = np.asarray(Temp_floatlist)
          TimeInterval = TimeScale[1]-TimeScale[0]
          # print(TimeInterval) #10.0
        if child3.tag == 'UNorm':
            for child4 in child3:
              # print('  >  >  > ', child4.tag)
              # print('  >  >  > ', child4.attrib)
              para_array.append(float(child4.attrib['ActuatorParameter']))
              Temp_txtlist = child4.text.split(' ')
              Temp_floatlist = [float(item) for item in Temp_txtlist]
              UNorm_list.append(np.asarray(Temp_floatlist))
            para_array = np.asarray(para_array)
            # print(para_array)
            # print(UNorm_list)
            # print(np.shape(self.UNorm))
        if child3.tag == 'XNorm':
          for child4 in child3:
            Temp_txtlist = child4.text.split(' ')
            Temp_floatlist = [float(item) for item in Temp_txtlist]
            XNorm_list.append(np.asarray(Temp_floatlist))
          # print(XNorm_list)
          # print(np.shape(self.XNorm))
        if child3.tag == 'XLast':
          for child4 in child3:
            Temp_txtlist = child4.text.split(' ')
            Temp_floatlist = [float(item) for item in Temp_txtlist]
            XLast_list.append(np.asarray(Temp_floatlist))
          # print(XLast_list)
          # print(np.shape(self.XLast))
        if child3.tag == 'YNorm':
          for child4 in child3:
            Temp_txtlist = child4.text.split(' ')
            Temp_floatlist = [float(item) for item in Temp_txtlist]
            YNorm_list.append(np.asarray(Temp_floatlist))
          # print(YNorm_list)
          # print(YNorm_list[0])
          # print(np.shape(YNorm_list))
          # print(np.shape(self.YNorm))
        for child4 in child3:
          for child5 in child4:
            # print('  >  >  > ', child5.tag) # real, imaginary, matrixShape, formatNote                     
            if child5.tag == 'real':
              Temp_txtlist = child5.text.split(' ')
              Temp_floatlist = [float(item) for item in Temp_txtlist]
              # print(Temp_txtlist)
              # print(Temp_floatlist)
              if child3.tag == 'Atilde':
                A_Re_list.append(np.asarray(Temp_floatlist))
              if child3.tag == 'Btilde':
                B_Re_list.append(np.asarray(Temp_floatlist))
              if child3.tag == 'Ctilde':
                C_Re_list.append(np.asarray(Temp_floatlist))

            if child5.tag == 'imaginary':
              Temp_txtlist = child5.text.split(' ')
              Temp_floatlist = [float(item) for item in Temp_txtlist]
              # print(Temp_txtlist)
              # print(Temp_floatlist)
              if child3.tag == 'Atilde':
                A_Im_list.append(np.asarray(Temp_floatlist))
              if child3.tag == 'Btilde':
                B_Im_list.append(np.asarray(Temp_floatlist))
              if child3.tag == 'Ctilde':
                C_Im_list.append(np.asarray(Temp_floatlist))

  # print(A_Re_list)
  # print(C_Im_list)
  n = len(XNorm_list[0]) # dimension of x
  m = len(UNorm_list[0]) # dimension of u
  p = len(YNorm_list[0]) # dimension of y

  # Reshape the A, B, C lists
  for i in range(len(para_array)):
      A_Re_list[i]=np.reshape(A_Re_list[i],(n,n)).T
      A_Im_list[i]=np.reshape(A_Im_list[i],(n,n)).T
      B_Re_list[i]=np.reshape(B_Re_list[i],(m,n)).T
      B_Im_list[i]=np.reshape(B_Im_list[i],(m,n)).T
      C_Re_list[i]=np.reshape(C_Re_list[i],(n,p)).T
      C_Im_list[i]=np.reshape(C_Im_list[i],(n,p)).T

  # print(A_Re_list[19])
  # print(B_Re_list[19])
  # print(C_Re_list[19])

  A_list = A_Re_list
  B_list = B_Re_list
  C_list = C_Re_list

  eig_A_array=[]
  # eigenvalue of A
  for i in range(len(para_array)):
      w,v = np.linalg.eig(A_list[i])
      eig_A_array.append(max(w))
  eig_A_array = np.asarray(eig_A_array)
  # print(eig_A_array)
  
  return TimeInterval, n, m, p, para_array, UNorm_list, XNorm_list, XLast_list, YNorm_list, A_list, B_list, C_list, eig_A_array

def check_YNorm_within_Range(y_min, y_max, para_array, UNorm_list, XNorm_list, XLast_list, YNorm_list, A_list, B_list, C_list, eig_A_array):
  UNorm_list_ = []; XNorm_list_ = []; XLast_list_ = []; YNorm_list_ =[]
  A_list_ = []; B_list_ = []; C_list_ = []; para_array_ = []; eig_A_array_ =[]

  for i in range(len(YNorm_list)):
    state = True
    for j in range(len(YNorm_list[i])):
      if YNorm_list[i][j] < y_min[j] or YNorm_list[i][j] > y_max[j]:
        state = False
    if state == True:
      UNorm_list_.append(UNorm_list[i])
      XNorm_list_.append(XNorm_list[i])
      XLast_list_.append(XLast_list[i])
      YNorm_list_.append(YNorm_list[i])
      A_list_.append(A_list[i])
      B_list_.append(B_list[i])
      C_list_.append(C_list[i])
      para_array_.append(para_array[i])
      eig_A_array_.append(eig_A_array[i])

  para_array_ = np.asarray(para_array_); eig_A_array_ = np.asarray(eig_A_array_)
  return para_array_, UNorm_list_, XNorm_list_, XLast_list_, YNorm_list_, A_list_, B_list_, C_list_, eig_A_array_

def fun_MOAS_noinf(A, B, C, D, s, g):
  p = len(C)  # dimension of y
  T = np.linalg.solve(np.identity(len(A))-A, B)
  """ Build the S matrix"""
  S = np.zeros((2*p, p))
  for i in range(0,p):
    S[2*i, i] = 1.0
    S[2*i+1, i] = -1.0
  Kx = np.dot(S,C)
  # print("Kx", Kx)
  Lim = np.dot(S,(np.dot(C,T) + D))
  # print("Lim", Lim)
  Kr = np.dot(S,D)
  # print("Kr", Kr)
  """ Build the core of H and h """
  # H = np.concatenate((0*Kx, Lim),axis=1); h = s
  # NewBlock = np.concatenate((Kx, Kr),axis=1)
  # H = np.concatenate((H,NewBlock)); h = np.concatenate((h,s))
  H = np.concatenate((Kx, Kr),axis=1); h = s

  """ Build the add-on blocks of H and h """
  i = 0
  while i < g :
    i = i + 1
    Kx = np.dot(Kx, A)
    Kr = Lim - np.dot(Kx,T)

    NewBlock = np.concatenate((Kx,Kr), axis=1)
    H = np.concatenate((H,NewBlock)); h = np.concatenate((h,s))
    """ To Insert the ConstRedunCheck """

  return H, h

def fun_RG_SISO(v_0, x, r, H, h, p):
  n = len(x) # dimension of x
  x = np.vstack(x) # x is horizontal array, must convert to vertical for matrix operation
  # because v_0 and r are both scalar, so no need to vstack
  Hx = H[:, 0:n]; Hv = H[:, n:]
  alpha = h - np.dot(Hx,x) - np.dot(Hv,v_0) # alpha is the system remaining vector
  beta = np.dot(Hv, (r-v_0)) # beta is the anticipated response vector with r

  kappa = 1
  for k in range(0,len(alpha)):
    if 0 < alpha[k] and alpha[k] < beta[k]:
      kappa = min(kappa, alpha[k]/beta[k])     
    else:
      kappa = kappa
  v = np.asarray(v_0 + kappa*(r-v_0)).flatten()

  return v

def fun_2nd_gstep_calc(x, Hm, hm, A_m, B_m, g):
  n = len(x) # dimension of x
  # x = np.vstack(x) # x is horizontal array, must convert to vertical for matrix operation
  # because v_0 and r are both scalar, so no need to vstack
  Hxm = Hm[:, 0:n]; Hvm = Hm[:, n:]

  T = np.linalg.solve(np.identity(n)-A_m, B_m)
  Ag = np.identity(n)
  for k in range(g+1):
      Ag = np.dot(Ag,A_m)
  
  alpha = hm - np.dot(Hxm, np.dot(Ag, np.vstack(x)))
  beta = np.dot(Hxm, np.dot((np.identity(n)-Ag),T))
  # print(np.shape(alpha))
  # print(np.shape(beta))
  v_st = []; v_bt = []
  for k in range(0,len(alpha)):
    if beta[k]>0:
      v_st.append(alpha[k]/beta[k])
    elif beta[k]<0:
      v_bt.append(alpha[k]/beta[k])
  # print('v_smaller_than,\n',v_st)
  v_max = np.asarray(min(v_st))
  v_min = np.asarray(max(v_bt))
  return v_max, v_min
