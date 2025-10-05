from _general import Enviroment
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random
from itertools import count
import threading
import time

class KineticalCalculator:
    def __init__(self , accuracy = 1e-3):
        self.accuracy = accuracy
        self.fitted = False
        self.calculations = 0
    def fit(self , enviroment):
        if type(enviroment) == Enviroment :
            self.enviroment = enviroment
        else:
            raise ValueError("The input should be an instance of Enviroment class")
        self.rate_constants = enviroment.rate_constants
        self.reactions_by_index = enviroment.reaction_by_index 
        self.stoichiometric_coefficient_by_reaction = enviroment.stoichiometric_coefficient_by_reaction
        self.rate_dependency_by_reaction = enviroment.rate_dependency_by_reaction
        self.number_of_reactions = len(enviroment)
        self.concentrations = []
        for compound in enviroment.compounds_concentration :
            self.concentrations.append(compound["concentration"])
        self.fitted = True
    def calculate(self  , time , checkpoint_time = [] , plot = False  ,concentration_below_zero = "SetToZero"):
        '''
        concentration_below_zero = "SetToZero" / "DoNotChange" / "GoNegetive"
        '''
        plt.figure()
        colors = []
        if plot:
            for i in self.enviroment.compounds :
                colors.append(((random.randint(0, 95)/100) , (random.randint(0, 95)/100) , (random.randint(0, 95)/100)))
                plt.xlabel("time")
                plt.ylabel("concentration")
                
        if not self.fitted :
            raise NameError("You should fir the model to an enviromt object before calculation")
        compounds_lenght = len(self.concentrations)
        checkpoints = [["time" , self.enviroment.compounds_unicode_formula]]
        def calculate_concentration_change():
            concentration_change = [0] * compounds_lenght
            for rxn_index in range(self.number_of_reactions) :
                rf = self.accuracy * self.rate_constants[rxn_index][0]
                counter = 0
                for compound in self.reactions_by_index[rxn_index][0]:
                    if not(self.concentrations[compound] == 0 and self.rate_dependency_by_reaction[rxn_index][0][counter] < 0):
                        rf = rf * (self.concentrations[compound] ** self.rate_dependency_by_reaction[rxn_index][0][counter])
                    else: 
                        rf = 0
                    counter += 1
                counter = 0
                rb = self.accuracy * self.rate_constants[rxn_index][1]
                for compound in self.reactions_by_index[rxn_index][1]:
                    if not(self.concentrations[compound] == 0 and self.rate_dependency_by_reaction[rxn_index][1][counter] < 0):
                        rb = rb * (self.concentrations[compound] ** self.rate_dependency_by_reaction[rxn_index][1][counter])
                    else:
                        rb = 0
                    counter += 1
                counter = 0
                for reactant in self.reactions_by_index[rxn_index][0]:
                    concentration_change[reactant] += (rb - rf) * self.stoichiometric_coefficient_by_reaction[rxn_index][0][counter] 
                    counter += 1
                counter = 0
                for product in self.reactions_by_index[rxn_index][1]:
                    concentration_change[product] += (rf - rb) * self.stoichiometric_coefficient_by_reaction[rxn_index][1][counter]
                    counter += 1
            return concentration_change
        t = 0
        for i in range(int(time/self.accuracy+1)) :
            delta_concentration = calculate_concentration_change()
            for j in range(len(delta_concentration)):
                concentration = self.concentrations[j] + delta_concentration[j]
                if concentration > 0 :
                    self.concentrations[j] = concentration
                elif concentration_below_zero == "DoNotChange":
                    pass
                elif concentration_below_zero == "GoNegetive":
                    self.concentrations[j] = concentration
                elif concentration_below_zero == "SetToZero":
                    self.concentrations[j] = 0
                else:
                    raise ValueError("concentration_below_zero should have one of the following values: SetToZero\\DoNotChange\\GoNegetive")
            if plot :
                for k in range(len(self.concentrations)):
                    plt.plot([t , t-self.accuracy],[self.concentrations[k] , self.concentrations[k]-delta_concentration[k]] , color = colors[k])
            for checkpoint_t in checkpoint_time:
                
                if t <= checkpoint_t < t + self.accuracy:
                    checkpoints.append([checkpoint_t , self.concentrations.copy()])
                    
            t += self.accuracy
        if plot :
            for k in range(len(self.concentrations)):
                plt.plot([0 , 0],[0 , 0] , color = colors[k], label = self.enviroment.compounds[k].unicode_formula)
            plt.legend()
            plt.show(block = False)
            
            print("type exit to continue:")
            exited = False
            while not exited:
                input_text = input()
                if input_text == "exit":
                    plt.close()
                    break
                else:
                    print("invalid_input")
        checkpoints.append([time , self.concentrations])
        
        
        return checkpoints
    
    def calculate_responsively(self, checkpoint_time = [] ,  animation_update_interval = 0.1 , plot = False , concentration_below_zero = "SetToZero"):
        colors = []
        if plot:
            for i in self.enviroment.compounds :
                colors.append(((random.randint(0, 95)/100) , (random.randint(0, 95)/100) , (random.randint(0, 95)/100)))
                plt.xlabel("time")
                plt.ylabel("concentration")
                
        if not self.fitted :
            raise NameError("You should fir the model to an enviromt object before calculation")
        compounds_lenght = len(self.concentrations)
        checkpoints = [["time" , self.enviroment.compounds_unicode_formula]]
        def calculate_concentration_change():
            concentration_change = [0] * compounds_lenght
            for rxn_index in range(self.number_of_reactions) :
                rf = self.accuracy * self.rate_constants[rxn_index][0]
                counter = 0
                for compound in self.reactions_by_index[rxn_index][0]:
                    if not(self.concentrations[compound] == 0 and self.rate_dependency_by_reaction[rxn_index][0][counter] < 0):
                        rf = rf * (self.concentrations[compound] ** self.rate_dependency_by_reaction[rxn_index][0][counter])
                    else: 
                        rf = 0
                    counter += 1
                counter = 0
                rb = self.accuracy * self.rate_constants[rxn_index][1]
                for compound in self.reactions_by_index[rxn_index][1]:
                    if not(self.concentrations[compound] == 0 and self.rate_dependency_by_reaction[rxn_index][1][counter] < 0):
                        rb = rb * (self.concentrations[compound] ** self.rate_dependency_by_reaction[rxn_index][1][counter])
                    else:
                        rb = 0
                    counter += 1
                counter = 0
                for reactant in self.reactions_by_index[rxn_index][0]:
                    concentration_change[reactant] += (rb - rf) * self.stoichiometric_coefficient_by_reaction[rxn_index][0][counter] 
                    counter += 1
                counter = 0
                for product in self.reactions_by_index[rxn_index][1]:
                    concentration_change[product] += (rf - rb) * self.stoichiometric_coefficient_by_reaction[rxn_index][1][counter]
                    counter += 1
            return concentration_change
        time = count()
        def animate(i):
            t = self.accuracy * next(time)
            delta_concentration = calculate_concentration_change()
            for j in range(len(delta_concentration)):
                concentration = self.concentrations[j] + delta_concentration[j]
                if concentration > 0 :
                    self.concentrations[j] = concentration
                elif concentration_below_zero == "DoNotChange":
                    pass
                elif concentration_below_zero == "GoNegetive":
                    self.concentrations[j] = concentration
                elif concentration_below_zero == "SetToZero":
                    self.concentrations[j] = 0
                else:
                    raise ValueError("concentration_below_zero should have one of the following values: SetToZero\\DoNotChange\\GoNegetive")
            if plot :
                for k in range(len(self.concentrations)):
                    plt.plot([t , t-self.accuracy],[self.concentrations[k] , self.concentrations[k]-delta_concentration[k]] , color = colors[k])
            for checkpoint_t in checkpoint_time:
                
                if t <= checkpoint_t < t + self.accuracy:
                    checkpoints.append([checkpoint_t , self.concentrations.copy()])


        ani = FuncAnimation(plt.gcf() , animate , interval = animation_update_interval , cache_frame_data=False)
        if plot :
            for k in range(len(self.concentrations)):
                plt.plot([0 , 0],[0 , 0] , color = colors[k], label = self.enviroment.compounds[k].unicode_formula)

            plt.legend()
            plt.show(block = False)
            print("type exit to continue / stop to pause the plotting / resume to resume the plotting:")
            exited = False
            while not exited:
                input_text = input()
                if input_text == "exit":
                    ani.pause()
                    plt.close()
                    break
                elif input_text == "stop":
                    ani.pause()
                elif input_text == "resume":
                    ani.resume()
                else:
                    print("invalid_input")
        checkpoints.append([time , self.concentrations])
        return checkpoints