# Copyright 2017 Battelle Energy Alliance, LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Created on June 04, 2022
@author: khnguy22 NCSU
updated on June 29, 2023
@author: Juan C. Luque-Gutierrez NCSU


comments: Interface for SIMULATE3 loading pattern optimzation
"""
import os
import numpy
from statistics import mean # Needed for Avg Kinf at EOC Calculation. 

class SimulateData:
  """
  Class that parses output of SIMULATE3 for a multiple run
  Partially copied from MOF work
  """
  def __init__(self,filen):
    """
      Constructor
      @ In, filen, string or dict, file name to be parsed, read one file at a time?
      @ Out, None
    """
    self.data = {} # dictionary to store the data from model output file
    self.lines = open(os.path.abspath(os.path.expanduser(filen)),"r").readlines() # raw data from model output
    # retrieve data
    self.data['axial_mesh'] = self.axialMeshExtractor()
    self.data['keff'] = self.coreKeffEOC()
    self.data['FDeltaH'] = self.maxFDH()
    self.data["kinf"] = self.kinfEOC()
    self.data["boron"] = self.boronEOC()
    self.data["cycle_length"] = self.EOCEFPD()
    self.data["PinPowerPeaking"] = self.pinPeaking()
    self.data["exposure"] = self.burnupEOC()
    self.data["assembly_power"] = self.assemblyPeakingFactors()
    self.data["fuel_type"] = self.fa_type()
    self.data["fuelcost_invEFPD"] = self.FuelCostEFPD()
    self.data["ArtObjOne"] = self.ArtificialObjectiveOne()
    self.data["ArtObjTwo"] = self.ArtificialObjectiveTwo()
    self.data["ConstCompl"] = self.ConstraintCompliance()
    self.data["Avg_Kinf_type_EOC"] = self.Avg_Kinf_EOC()

#    self.data["pin_peaking"] = self.pinPeaking()
    # this is a dummy variable for demonstration with MOF
    # check if something has been found
    if all(v is None for v in self.data.values()):
      raise IOError("No readable outputs have been found!")
#------------------------------------------------------------------------------------------
  #function to retrivedata
  def getPin(self):
    """
      Retrive total number of pins
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                            {'info_ids':list(of ids of data),
                              'values': list}
    """
    outputDict = None
    for line in self.lines:
      if line.find('Assembly Core Maps . . . .')>=0:
        temp = line.strip().split('(')
        temp = temp[1].split(',')[0]
        break
    outputDict = {'info_ids':['pin_number'], 'values': [int(temp)] }
    return outputDict

  def axialMeshExtractor(self):
    """
      Extracts the axial mesh used in the SIMULATE output file.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                            {'info_ids':list(of ids of data),
                              'values': list}
    """
    outputDict = None
    reverseAxialPositions = [] #SIMULATE reversed lists axial meshes from bottom to top
    searchingHeights = False
    for line in self.lines:
      if "** Studsvik CMS Steady-State 3-D Reactor Simulator **" in line:
        searchingHeights = False
      if "Grid Location Information" in line:
        searchingHeights = False
        break
      if searchingHeights:
        line = line.replace("-","")
        elems = line.strip().split()
        if elems:
          reverseAxialPositions.append(float(elems[-1]))
      if "Axial Nodal Boundaries (cm)" in line:
          searchingHeights = True
    #The bot/top axial node in the reflectors are not considered
    reverseAxialPositions.pop(0)
    reverseAxialPositions.pop(-1)

    forwardAxialPositions = []
    for position in reverseAxialPositions:
      forwardAxialPositions.insert(0,position)

    outputDict = {'info_ids':['no_axial_node'],
                  'values': [len(forwardAxialPositions)] }

    return outputDict

  def getCoreWidth(self):
    """
      Retrive core width
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                            {'info_ids':list(of ids of data),
                              'values': list}
    """
    outputDict = None
    for line in self.lines:
      if line.find("'DIM.PWR'")>=0:
        temp = line.strip().split(' ')
        temp = temp[1].split('/')[0]
        break
    outputDict = {'info_ids':['core_width'], 'values': [int(temp)] }
    return outputDict


  def coreKeffEOC(self):
    """
      Extracts the core K-effective value from the provided simulate file lines.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                            {'info_ids':list(of ids of data),
                              'values': list}
    """
    keffList = []
    outputDict = None
    for line in self.lines:
      if "K-effective . . . . . . . . . . . . ." in line:
        elems = line.strip().split()
        keffList.append(float(elems[-1]))
    if not keffList:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['eoc_keff'], 'values': [keffList[-1]] }
    return outputDict


  def assemblyPeakingFactors(self):
    """
      Extracts the assembly radial power peaking factors as a dictionary
      with the depletion step in GWD/MTU as the dictionary keys.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                            {'info_ids':list(of ids of data),
                              'values': list}
    """
    radialPowerDictionary = {}
    searching_ = False
    outputDict = None
    for line in self.lines:
      if "Case" in line and "GWd/MT" in line:
        elems = line.strip().split()
        depl = elems[-2]
        if depl in radialPowerDictionary:
          pass
        else:
          radialPowerDictionary[depl] = {}
      if "**   H-     G-     F-     E-     D-     C-     B-     A-     **" in line:
        searching_ = False

      if searching_:
        elems = line.strip().split()
        if elems[0] == "**":
          posList = elems[1:-1]
        else:
          radialPowerDictionary[depl][elems[0]] = {}
          for i,el in enumerate(elems[1:-1]):
            radialPowerDictionary[depl][elems[0]][posList[i]] = float(el)

      if "PRI.STA 2RPF  - Assembly 2D Ave RPF - Relative Power Fraction" in line:
        searching_ = True

    if not radialPowerDictionary:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      maxPeaking = 0.0
      for depl in radialPowerDictionary:
        for row in radialPowerDictionary[depl]:
          for col in radialPowerDictionary[depl][row]:
            maxPeaking = max(radialPowerDictionary[depl][row][col],maxPeaking)
      outputDict = {'info_ids':['FA_peaking'], 'values': [maxPeaking] }

    return outputDict

  def EOCEFPD(self):
    """
      Returns maximum of EFPD values for cycle exposure in the simulate
      file.

      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    list_ = []
    outputDict = None
    for line in self.lines:
      if "Cycle Exp." in line:
        if "EFPD" in line:
          elems = line.strip().split()
          spot = elems.index('EFPD')
          list_.append(float(elems[spot-1]))

    if not list_:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['MaxEFPD'], 'values': [list_[-1]]}

    return outputDict

  def maxFDH(self):
    """
      Returns maximum of F-delta-H values for each cycle exposure in the simulate
      file.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    list_ = []
    outputDict = None
    for line in self.lines:
      if "F-delta-H" in line:
        elems = line.strip().split()
        spot = elems.index('F-delta-H')
        list_.append(float(elems[spot+1]))

    if not list_:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['MaxFDH'], 'values': [max(list_)] }

    return outputDict

  def pinPeaking(self):
    """
      Returns maximum value of pin peaking values, Fq, for each cycle exposure in the simulate
      file.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    outputDict = None
    list_ = []
    for line in self.lines:
      if "Max-3PIN" in line:
        elems = line.strip().split()
        spot = elems.index('Max-3PIN')
        list_.append(float(elems[spot+1]))

    if not list_:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['pin_peaking'], 'values': [max(list_)] }

    return outputDict

  def boronEOC(self):
    """
      Returns EOC and max boron values in PPM at each depletion step.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    boronList = []
    outputDict = None
    for line in self.lines:
      if "Boron Conc." in line and "ppm" in line:
        elems = line.strip().split()
        spot = elems.index('ppm')
        boronList.append(float(elems[spot-1]))

    if not boronList:
      return ValueError("NO values returned. Check SIMULATE file executed correctly")
    else:
      outputDict = {'info_ids':['eoc_boron', 'max_boron'],
                    'values': [boronList[-1], max(boronList)] }

    return outputDict

  def kinfEOC(self):
    """
      Returns a list of kinf values from Simulate3.
      Only work for PWR
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    kinfList = []
    searchingForKinf = False
    outputDict = None
    for line in self.lines:
      elems = line.strip().split()
      if not elems:
        pass
      else:
        if searchingForKinf:
          if elems[0] == '1':
            kinfList.append(float(elems[1]))
            searchingForKinf = False
        if "PRI.STA 2KIN  - Assembly 2D Ave KINF - K-infinity" in line:
          searchingForKinf = True

    if not kinfList:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['eoc_kinf'], 'values': [ kinfList[-1]] }

    return outputDict

  def relativePower(self):
    """
      Extracts the Relative Core Power from the provided simulate file lines.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    relativePowers = []
    outputDict = None
    for line in self.lines:
      if "Relative Power. . . . . . .PERCTP" in line:
        p1 = line.index("PERCTP")
        p2 = line.index("%")
        searchSpace = line[p1:p2]
        searchSpace = searchSpace.replace("PERCTP","")
        relativePowers.append(float(searchSpace))

    if not relativePowers:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['relative power'], 'values': [relativePowers] }

    return outputDict

  def relativeFlow(self):
    """
      Extracts the Relative Core Flow rate from the provided simulate file lines.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    relativeFlows = []
    outputDict = None
    for line in self.lines:
      if "Relative Flow . . . . . .  PERCWT" in line:
        p1 = line.index("PERCWT")
        p2 = line.index("%")
        searchSpace = line[p1:p2]
        searchSpace = searchSpace.replace("PERCWT","")
        relativeFlows.append(float(searchSpace))

    if not relativeFlows:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['relative flow'], 'values': [relativeFlows] }

    return outputDict

  def thermalPower(self):
    """
      Extracts the operating thermal power in MW from the provided simulate file lines.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    powers = []
    outputDict = None
    for line in self.lines:
      if "Thermal Power . . . . . . . . CTP" in line:
        elems = line.strip().split()
        spot = elems.index('MWt')
        powers.append(float(elems[spot-1]))

    if not powers:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['thermal power'], 'values': [powers] }

    return outputDict

  def coreFlow(self):
    """
      Returns the core coolant flow in Mlb/hr from the provided simulate file lines.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    flows = []
    outputDict = None
    for line in self.lines:
      if "Core Flow . . . . . . . . . . CWT" in line:
        elems = line.strip().split()
        spot = elems.index("Mlb/hr")
        flows.append(float(elems[spot-1]))

    if not flows:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['core flow'], 'values': [flows] }

    return outputDict

  def inletTemperatures(self):
    """
      Returns the core inlet temperatures in degrees Fahrenheit from the
      provided simulate file lines.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    temperatures = []
    outputDict = None
    for line in self.lines:
      if "Inlet . . . .TINLET" in line:
        p1 = line.index("K")
        p2 = line.index("F")
        searchSpace = line[p1:p2]
        searchSpace = searchSpace.replace("K","")
        temperatures.append(float(searchSpace))

    if not temperatures:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['inlet temperatures'], 'values': [temperatures] }

    return outputDict

  def pressure(self):
    """
      Returns the core exit pressure in PSIA.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    pressure = []
    outputDict = None
    for line in self.lines:
      if "Core Exit Pressure  . . . . .  PR" in line:
        p1 = line.index("bar")
        p2 = line.index("PSIA")
        searchSpace = line[p1:p2]
        searchSpace = searchSpace.replace("bar","")
        pressure.append(float(searchSpace))

    if not pressure:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['pressure'], 'values': [pressure] }

    return outputDict

  def burnupEOC(self):
    """
      Extracts the cycle burnups at a each state point within the depletion.
      @ In, None
      @ Out, outputDict, dict, the dictionary containing the read data (None if none found)
                        {'info_ids':list(of ids of data),
                          'values': list}
    """
    burnups = []
    for line in self.lines:
      if "Cycle Exp." in line:
        elems = line.strip().split()
        spot = elems.index('GWd/MT')
        burnups.append(float(elems[spot-1]))
    if not burnups:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['exposure'], 'values': [burnups[-1]] }

    return outputDict

  def fa_type(self):
    '''
    Extracts the fuel type and calculates the fuel cost based on the amount and enrichment of each fuel type.
    '''
    #fuel_type = []
    FAlist = []
    for line in self.lines:
      if "'FUE.TYP'" in line:
        p1 = line.index(",")
        p2 = line.index("/")
        search_space = line[p1:p2]
        search_space = search_space.replace(",","")
        tmp= search_space.split()
        for ii in tmp:
          FAlist.append(float(ii))
    FAtype = list(set(FAlist))
    FAlist_A = FAlist[0]
    FAlist_B = FAlist[1:9] + FAlist[9:73:9]
    FAlist_C = FAlist[10:18] + FAlist[19:27] + FAlist[28:36] + FAlist[37:45] + FAlist[46:54] + FAlist[55:63] + FAlist[64:72] + FAlist[73:81]
    FAcount_A = [float(fa == FAlist_A) for fa in FAtype]
    FAcount_B = [float(FAlist_B.count(fa)*2) for fa in FAtype]
    FAcount_C = [float(FAlist_C.count(fa)*4) for fa in FAtype]
    FAcount = [FAcount_A[j] + FAcount_B[j] + FAcount_C[j] for j in range(len(FAtype))]
    FA_type_dict = {int(FAtype[i]):FAcount[i] for i in range(len(FAtype))}
    print(FAcount)
    #stop
    #Considering that: FA type 0 is empty, type 1 reflector, type 2 2% enrichment, types 3 and 4 2.5% enrichment, and types 5 and 6 3.2% enrichment. The cost of burnable is not being considered
    #if len(FAcount) == 7:
    #  fuel_cost = (FAcount[0] + FAcount[1])*0 + FAcount[2]*2.69520839 + (FAcount[3] + FAcount[4])*3.24678409 + (FAcount[5] + FAcount[6])*4.03739539
    #else:
    #  fuel_cost = (FAcount[0] + FAcount[1])*0 + FAcount[2]*2.69520839 + (FAcount[3] + FAcount[4])*3.24678409 + (FAcount[5])*4.03739539
    #print(fuel_cost)
    #fuel_type.append(float(search_space))
    #stop
    
    # Dictionary with the unit cost for each FA type.

    # FA type 0 = empty         -> M$ 0.0
    # FA type 1 = reflector     -> M$ 0.0
    # FA type 2 = 2.00 wt%      -> M$ 2.69520839
    # FA type 3 = 2.50 wt%      -> M$ 3.24678409
    # FA type 4 = 2.50 wt% + Gd -> M$ 3.24678409
    # FA type 5 = 3.20 wt%      -> M$ 4.03739539
    # FA type 6 = 3.20 wt% + Gd -> M$ 4.03739539
    # The cost of burnable poison is not being considered.
    
    cost_dict = {
      0: 0,
      1: 0,
      2: 2.69520839,
      3: 3.24678409,
      4: 3.24678409,
      5: 4.03739539,
      6: 4.03739539
    }
    
    fuel_cost = 0
    for fuel_type, fuel_count in FA_type_dict.items():
      fuel_cost += fuel_count * cost_dict[fuel_type]
    
    if not fuel_cost:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['fuel_cost'], 'values': [fuel_cost]}
    return outputDict

  def FuelCostEFPD(self):
    """
    Artificial variable, product of fuel cost and cycle length. This artificial variable is used to explore a potential solution to the max cycle length min fuel cost problem. 

    """
    list_efpd = []
    outputDict = None
    for line in self.lines:
      if "Cycle Exp." in line:
        if "EFPD" in line:
          elems = line.strip().split()
          spot = elems.index('EFPD')
          list_efpd.append(float(elems[spot-1]))
          
    print(list_efpd[-1])

    FAlist = []
    for line in self.lines:
      if "'FUE.TYP'" in line:
        p1 = line.index(",")
        p2 = line.index("/")
        search_space = line[p1:p2]
        search_space = search_space.replace(",","")
        #print(search_space)
        tmp= search_space.split()
        for ii in tmp:
          FAlist.append(float(ii))
    FAtype = list(set(FAlist))
    FAlist_A = FAlist[0]
    FAlist_B = FAlist[1:9] + FAlist[9:73:9]
    FAlist_C = FAlist[10:18] + FAlist[19:27] + FAlist[28:36] + FAlist[37:45] + FAlist[46:54] + FAlist[55:63] + FAlist[64:72] + FAlist[73:81]
    FAcount_A = [float(fa == FAlist_A) for fa in FAtype]
    FAcount_B = [float(FAlist_B.count(fa)*2) for fa in FAtype]
    FAcount_C = [float(FAlist_C.count(fa)*4) for fa in FAtype]
    FAcount = [FAcount_A[j] + FAcount_B[j] + FAcount_C[j] for j in range(len(FAtype))]
    #print(FAcount)
    #stop
    #Considering that: FA type 0 is empty, type 1 reflector, type 2 2% enrichment, types 3 and 4 2.5% enrichment, and types 5 and 6 3.2% enrichment. The cost of burnable is not being considered
    if len(FAcount) == 7:
      fuel_cost = (FAcount[0] + FAcount[1])*0 + FAcount[2]*2.69520839 + (FAcount[3] + FAcount[4])*3.24678409 + (FAcount[5] + FAcount[6])*4.03739539
    else:
      fuel_cost = (FAcount[0] + FAcount[1])*0 + FAcount[2]*2.69520839 + (FAcount[3] + FAcount[4])*3.24678409 + (FAcount[5])*4.03739539
    print(fuel_cost)
    
    fuelcost_invEFPD = fuel_cost/list_efpd[-1]
    
    print(fuelcost_invEFPD)
    
    if not fuelcost_invEFPD:
      return ValueError("Verify interface script to see if there is any error")
    
    else:
      outputDict = {'info_ids':['fuelcost_invEFPD'], 'values':[fuelcost_invEFPD]}
      return outputDict

  def ArtificialObjectiveOne(self):
    """ 
    Artificial objective to include penalty. This variable is associated with the Cycle Length for min-min optimization (therefore, objective -Cycle Lenght)
    """
    
    list_efpd = []
    outputDict = None
    for line in self.lines:
      if "Cycle Exp." in line:
        if "EFPD" in line:
          elems = line.strip().split()
          spot = elems.index('EFPD')
          list_efpd.append(float(elems[spot-1]))
    EFPD = list_efpd[-1]
    print(f"ARTOBJONE - EFPD is equal to: {EFPD}")
    
    
    boronList = []
    for line in self.lines:
      if "Boron Conc." in line and "ppm" in line:
        elems = line.strip().split()
        spot = elems.index('ppm')
        boronList.append(float(elems[spot-1]))
    max_boron = max(boronList)
    print(f"ARTOBJONE - max_boron is equal to: {max_boron}")  
    
    g1 = 1.3 - max_boron/1000
    
    if g1 < 0:
      g1 = abs(g1)
    else:
      g1 = 0
    print(f"This is g1: {g1}")
    
    pin_peaking_list = []
    for line in self.lines:
      if "Max-3PIN" in line:
        elems = line.strip().split()
        spot = elems.index('Max-3PIN')
        pin_peaking_list.append(float(elems[spot+1]))
    
    pin_peaking = max(pin_peaking_list)
    print(f"ARTOBJONE - pin_peaking is equal to: {pin_peaking}")    
    
    g2 = 2.1 - pin_peaking
    
    if g2 < 0:
      g2 = abs(g2)
    else:
      g2 = 0
    print(f"This is g2: {g2}")
    
    list_FDH = []
    for line in self.lines:
      if "F-delta-H" in line:
        elems = line.strip().split()
        spot = elems.index('F-delta-H')
        list_FDH.append(float(elems[spot+1]))
    
    MaxFDH = max(list_FDH)
    print(f"ARTOBJONE - MaxFDH is equal to: {MaxFDH}")    
    
    g3 = 1.48 - MaxFDH
    
    if g3 < 0:
      g3 = abs(g3)
    else:
      g3 = 0
    print(f"This is g3: {g3}")
    
    w = (4*g1 + g2 + g3)*700 #P1 Penalty for Artificial Objective 1 = 70. One value assigned for sum of degree of violation of the three constraints.
    
    ArtObjOne = -list_efpd[-1] + w
    print(f"THIS IS THE ARTIFICIAL OBJECTIVE ONE:{ArtObjOne}")
    
    if not ArtObjOne:
      return ValueError("Verify that all values involved have been parsed correctly")
    
    else:
      outputDict = {'info_ids':['ArtObjOne'], 'values':[ArtObjOne]}
    return outputDict

  def ArtificialObjectiveTwo(self):
    """
    Artificial objective to include penalty. This variable is associated with Fuel Cost for min-min optimization.
    """
    FAlist = []
    for line in self.lines:
      if "'FUE.TYP'" in line:
        p1 = line.index(",")
        p2 = line.index("/")
        search_space = line[p1:p2]
        search_space = search_space.replace(",","")
        #print(search_space)
        tmp= search_space.split()
        for ii in tmp:
          FAlist.append(float(ii))
    FAtype = list(set(FAlist))
    FAlist_A = FAlist[0]
    FAlist_B = FAlist[1:9] + FAlist[9:73:9]
    FAlist_C = FAlist[10:18] + FAlist[19:27] + FAlist[28:36] + FAlist[37:45] + FAlist[46:54] + FAlist[55:63] + FAlist[64:72] + FAlist[73:81]
    FAcount_A = [float(fa == FAlist_A) for fa in FAtype]
    FAcount_B = [float(FAlist_B.count(fa)*2) for fa in FAtype]
    FAcount_C = [float(FAlist_C.count(fa)*4) for fa in FAtype]
    FAcount = [FAcount_A[j] + FAcount_B[j] + FAcount_C[j] for j in range(len(FAtype))]
    #print(FAcount)
    #stop
    #Considering that: FA type 0 is empty, type 1 reflector, type 2 2% enrichment, types 3 and 4 2.5% enrichment, and types 5 and 6 3.2% enrichment. The cost of burnable is not being considered
    if len(FAcount) == 7:
      fuel_cost = (FAcount[0] + FAcount[1])*0 + FAcount[2]*2.69520839 + (FAcount[3] + FAcount[4])*3.24678409 + (FAcount[5] + FAcount[6])*4.03739539
    else:
      fuel_cost = (FAcount[0] + FAcount[1])*0 + FAcount[2]*2.69520839 + (FAcount[3] + FAcount[4])*3.24678409 + (FAcount[5])*4.03739539
    print(f"ARTOBJTWO - fuel_cost is equal to: {fuel_cost}")

    boronList = []
    for line in self.lines:
      if "Boron Conc." in line and "ppm" in line:
        elems = line.strip().split()
        spot = elems.index('ppm')
        boronList.append(float(elems[spot-1]))
    max_boron = max(boronList)
    print(f"ARTOBJTWO - max_boron is equal to: {max_boron}")  
    
    g1 = 1.3 - max_boron/1000
    
    if g1 < 0:
      g1 = abs(g1)
    else:
      g1 = 0
    print(f"This is g1: {g1}")
    
    pin_peaking_list = []
    for line in self.lines:
      if "Max-3PIN" in line:
        elems = line.strip().split()
        spot = elems.index('Max-3PIN')
        pin_peaking_list.append(float(elems[spot+1]))
    
    pin_peaking = max(pin_peaking_list)
    print(f"ARTOBJTWO - pin_peaking is equal to: {pin_peaking}")    
    
    g2 = 2.1 - pin_peaking
    
    if g2 < 0:
      g2 = abs(g2)
    else:
      g2 = 0
    print(f"This is g2: {g2}")
    
    list_FDH = []
    for line in self.lines:
      if "F-delta-H" in line:
        elems = line.strip().split()
        spot = elems.index('F-delta-H')
        list_FDH.append(float(elems[spot+1]))
    
    MaxFDH = max(list_FDH)
    print(f"ARTOBJTWO - MaxFDH is equal to: {MaxFDH}")    
    
    g3 = 1.48 - MaxFDH
    
    if g3 < 0:
      g3 = abs(g3)
    else:
      g3 = 0
    print(f"This is g3: {g3}")
    
    w = (4*g1 + g2 + g3)*700 #P2 Penalty for Artificial Objective 1 = 70. One value assigned for sum of degree of violation of the three constraints.
     
    ArtObjTwo = fuel_cost + w
    print(f"THIS IS THE ARTIFICIAL OBJECTIVE TWO: {ArtObjTwo}")

    if not ArtObjTwo:
      return ValueError("Verify that all values involved have been parsed correctly")
    
    else:
      outputDict = {'info_ids':['ArtObjTwo'], 'values':[ArtObjTwo]}
    return outputDict

  def ConstraintCompliance(self):
      
    """
    Function to count constraints compliance of each genome in each generation. 
    """
    Constraints = 0
    OutputDict = None 
    boronList = []
    for line in self.lines:
      if "Boron Conc." in line and "ppm" in line:
        elems = line.strip().split()
        spot = elems.index('ppm')
        boronList.append(float(elems[spot-1]))
    max_boron = max(boronList)    
    
    g1 = 1.3 - max_boron/1000
        
    if g1 < 0:
      Constraints += 1
    else:
      Constraints += 0
    
    pin_peaking_list = []
    for line in self.lines:
      if "Max-3PIN" in line:
        elems = line.strip().split()
        spot = elems.index('Max-3PIN')
        pin_peaking_list.append(float(elems[spot+1]))
    
    pin_peaking = max(pin_peaking_list)   
    
    g2 = 2.1 - pin_peaking
    
    if g2 < 0:
      Constraints += 1
    else:
      Constraints += 0
    
    list_FDH = []
    for line in self.lines:
      if "F-delta-H" in line:
        elems = line.strip().split()
        spot = elems.index('F-delta-H')
        list_FDH.append(float(elems[spot+1]))
    
    MaxFDH = max(list_FDH)
    
    g3 = 1.48 - MaxFDH
    
    if g3 < 0:
      Constraints += 1
    else:
      Constraints += 0

    if g1 + g2 + g3 == 0:
      Constraints = 0
    else:
      Constraints = Constraints

    print(f"The number of constraints violated is: {Constraints}")
    
    if not Constraints:
      return ValueError("Verify that constraints have been defined correctly")
    
    else:
      outputDict = {'info_ids':['Constraints'], 'values':[Constraints]}
    return outputDict

  def Avg_Kinf_EOC(self):
    """
    Return the average K-inf value for each fuel type at End-of-Cycle. 
    """
    kinf_map_dictionary = {}
    searching_ = False
    for line in self.lines:
      if "Case" in line and "GWd/MT" in line:
        elems = line.strip().split()
        depl = elems[-2]
        if depl in kinf_map_dictionary:
          pass
        else:
          kinf_map_dictionary[depl] = {}
      if "**   H-     G-     F-     E-     D-     C-     B-     A-     **" in line:
        searching_ = False
        
      if searching_:
        elems = line.strip().split()
        if elems[0] == "**":
          pos_list = elems[1:-1]
        else:
          kinf_map_dictionary[depl][elems[0]] = {}
          for i, el in enumerate(elems[1:-1]):
            kinf_map_dictionary[depl][elems[0]][pos_list[i]] = float(el)
            
      if  "PRI.STA 2KIN  - Assembly 2D Ave KINF - K-infinity" in line:
        searching_ = True
    print(f"Printing for debugging purposes: {kinf_map_dictionary}")
    
    EOC_depl = 0
    for key in kinf_map_dictionary.keys():
      depl = float(key)
      if depl > float(EOC_depl) and kinf_map_dictionary[key]:
        EOC_depl = key
    kinf_map_dictionary_EOC = kinf_map_dictionary[EOC_depl]
    
    fuel_type = []
    dict = {}
    
    FAlist = []
    for line in self.lines:
      if "'FUE.TYP'" in line:
        p1 = line.index(",")
        p2 = line.index("/")
        search_space = line[p1:p2]
        search_space = search_space.replace(",","")
        tmp = search_space.split()
        for i in tmp:
          FAlist.append(float(i))
    FAtypes = list(set(FAlist))
    
    fa_width = 15
    fa_center = fa_width // 2
    for i in range(len(FAlist)):
      row = i//9 + fa_center
      col = i%9 + fa_center
      key = (row, col)
      type = FAlist[i]
      kinf = None
      if type > 1:
        kinf = kinf_map_dictionary_EOC[str(row + 1)][str(col + 1)]
      dict[key] = {
        "type": type,
        "k-inf": kinf
      }
    for x in range(fa_width):
      for y in range(fa_width):
        target_x = x
        target_y = y
        if x < fa_center:
          target_x = (fa_center * 2) - x
        if y < fa_center:
          target_y = (fa_center * 2) - y
        dict[(x, y)] = dict[(target_x, target_y)]
    
    for x in range(fa_width):
      for y in range(fa_width):
        print(int(dict[(x, y)]["type"]), end="")
      print()
      
    kinfs_by_type = {}
    for key, value in dict.items():
      pos_x, pos_y = key
      type = value["type"]
      kinf = value["k-inf"]
      if type > 1:
        if type not in kinfs_by_type.keys():
          kinfs_by_type[type] = []
        kinfs_by_type[type].append(kinf)
    
    averages = {k:mean(v) for k, v in kinfs_by_type.items()}
    Avg_K_FA2 = averages[2]
    print(f"This is average K-inf at EOC for FA type 2: {Avg_K_FA2}")
    if not Avg_K_FA2:
      return ValueError("No values returned. Check Simulate File executed correctly")
    else:
      outputDict = {'info_ids':['avg k at eoc'], 'values': [Avg_K_FA2]}
    
      
  def writeCSV(self, fileout):
    """
      Print Data into CSV format
      @ In, fileout, str, the output file name
      @ Out, None
    """
    fileObject = open(fileout.strip()+".csv", mode='wb+') if not fileout.endswith('csv') else open(fileout.strip(), mode='wb+')
    headers=[]
    nParams = numpy.sum([len(data['info_ids']) for data in self.data.values() if data is not None and type(data) is dict])
    outputMatrix = numpy.zeros((nParams,1))
    index=0
    for data in self.data.values():
      if data is not None and type(data) is dict:
        headers.extend(data['info_ids'])
        for i in range(len(data['info_ids'])):
          outputMatrix[index]= data['values'][i]
          index=index+1
    numpy.savetxt(fileObject, outputMatrix.T, delimiter=',', header=','.join(headers), comments='')
    fileObject.close()
