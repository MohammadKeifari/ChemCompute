import re
class Compound: 
    def __init__(self , formula  , phase_point_list=None , mp=None, bp=None ,scription=True):
        '''
        Creates a Compound Object based on formula, physical properties and chemical properties

        Arg:
        formula(string): Compound chemical formula
        scription(boolean) : If True the formula will be superscipted and subscripted --> unicode_formula
        phase(list of dict): A list of phase point dictionaries. Each dictionary should have:
            - phase(str): defines the phase of the compund at that temperature. must be one of:
                - "s" : solid
                - "g" : gas
                - "l" : liquid
                - "aq" : aqueous 
            - temperature: the temperature of compound in that phase
        mp(int or float): melting point
        bp(int or float): boiling point
        To determine the phase phase_point_list will have priority over mp and bp
        '''
        superscript_characters=["\u2070" ,"\u00b9" ,"\u00b2" ,"\u00b3" ,"\u2074" 
                                ,"\u2075" ,"\u2076" ,"\u2077" ,"\u2078" ,"\u2079" 
                                ,"\u207a" , "\u207b"]
        subscript_characters = ["\u2080" , "\u2081", "\u2082", "\u2083", "\u2084"
                                , "\u2085", "\u2086", "\u2087", "\u2088", "\u2089" ]
        phases = ["g" , "l" , "s" , "aq"]

        self.formula = formula
        if scription:
            def formula_to_unicode_formula(formula):
                unicode_formula = ""
                if "+" in formula :
                    splitted_formula = formula.split("+")
                    for char in splitted_formula[0]:
                        if char.isdigit():
                            unicode_formula += subscript_characters[int(char)]
                        else:
                            unicode_formula += char
                    unicode_formula += superscript_characters[10]
                    for char in splitted_formula[1]:
                        if char.isdigit():
                            unicode_formula += superscript_characters[int(char)]
                        else:
                            unicode_formula += char
                elif "-" in formula :
                    splitted_formula = formula.split("-")
                    for char in splitted_formula[0]:
                        if char.isdigit():
                            unicode_formula += subscript_characters[int(char)]
                        else:
                            unicode_formula += char
                    unicode_formula += superscript_characters[11]
                    for char in splitted_formula[1]:
                        if char.isdigit():
                            unicode_formula += superscript_characters[int(char)]
                        else:
                            unicode_formula += char
                else:
                    for char in formula:
                        if char.isdigit():
                            unicode_formula += subscript_characters[int(char)]
                        else:
                            unicode_formula += char
                return unicode_formula
            self.unicode_formula = formula_to_unicode_formula(formula)
        else:
            self.unicode_formula = formula
        self.phase_point_list = []
        if phase_point_list != None:
            for phase_point in phase_point_list : 
                if phase_point["phase"] in phases:
                    self.phase_point_list.append(phase_point)
                else:
                    raise ValueError("The acceptable inputs for phase are s / l / g / aq")
        self.mp = mp
        self.bp = bp

    def phase(self , temperature):
        '''
        Returns the phase of a compound in the specific temperature

        Arg:
        temperature(int or float): The temperature you want to know the phase of the compound in

        Returns:
        phase:  the phase of a compound in the temperature
        '''
        phase_point_list_temperatures = [phase_point["temperature"] for phase_point in self.phase_point_list]
        if temperature in phase_point_list_temperatures:
            return (self.phase_point_list[phase_point_list_temperatures.index(temperature)])["phase"]
        elif self.bp != None and self.mp != None :
            if temperature <= self.mp :
                return "s" 
            elif self.bp >= temperature > self.mp:
                return "l"
            elif temperature > self.bp :
                return "g"
        elif self.bp == None and self.mp != None :
            if temperature <= self.mp :
                return "s"
            else:
                return "l"
        elif self.bp != None and self.mp == None :
            if temperature <= self.bp :
                return "l"
            else:
                return "g"
        elif self.bp == None and self.mp == None :
            return None
    def __str__(self):
        return self.unicode_formula
    def __eq__(self, value):
        if self.unicode_formula == value.unicode_formula :
            return True
        else:
            return False
class Reaction:
    def __init__(self):
        pass
    @classmethod
    def from_complex_syntax(self, reaction, K=1, kf=1, kb=1):
        '''
        Creates reaction object from an string with any of the following formats:
        A & 2B & ... > 3C & 2D-1 & ... \n
        The prefix number represents the stoichiometric_coefficient and the default value is 1 \n
        The suffix number represents rate_dependency and the default vaule is 1 \n
        example: Fe(CN)6-3 & Ce+2 > Fe(CN)6-4 & Ce+3

        To specify the phase use the following syntax:
        A.g & 2B.s & ... > 3c.aq & 2D-1.g & ...
        Args:
            reaction(string): The equation formula with the mentioned structure.

            K(int or float): The equilibrium constant of reaction

            kf(int or float): The forward reaction rate constant

            kb(int or float): The backward reaction rate constant

        '''
        self.K = K
        self.kf = kf
        self.kb = kb
        reformed_reaction = reaction.replace(" ","").split(">")
        splited_to_component_reaction = [component.split("&") for component in reformed_reaction] 
        component_counter = 0
        self.equation = reaction
        self.reactants = []
        self.products = []
        for component in splited_to_component_reaction:
            counter = 0
            for section in component:
                if not re.match(r'^(?:\d+(?:\.\d+)?_[A-Za-z0-9+.\-()]+_\d+(?:\.\d+)?|\d+(?:\.\d+)?_[A-Za-z0-9+.\-()]+|[A-Za-z0-9+.\-()]+_\d+(?:\.\d+)?|[A-Za-z0-9+.\-()]+)$' , section) :
                     raise ValueError
                
                splitted_section = section.split("_")
                lenght = len(splitted_section)
                if lenght == 3:
                    compound_info = {
                    "stoichiometric_coefficient" : float(splitted_section[0]),
                    "substance" : splitted_section[1],
                    "rate_dependency" : float(splitted_section[2])
                    }
                    if component_counter == 0:
                        self.reactants.append(compound_info)
                    elif component_counter == 1:
                        self.products.append(compound_info)
                elif lenght == 2:
                    if re.match(r'^\d+(?:\.\d+)?_[A-Za-z0-9+.\-()]+$' , section) :
                        compound_info = {
                            "stoichiometric_coefficient" : float(splitted_section[0]),
                            "substance" : splitted_section[1],
                            "rate_dependency" : 1
                        }
                        if component_counter == 0:
                            self.reactants.append(compound_info)
                        elif component_counter == 1:
                            self.products.append(compound_info)
                    elif re.match(r'^[A-Za-z0-9+.\-()]+_\d+(?:\.\d+)?$' , section) :
                        compound_info = {
                            "stoichiometric_coefficient" : 1,
                            "substance" : splitted_section[0],
                            "rate_dependency" : float(splitted_section[1])
                        }
                        if component_counter == 0:
                            self.reactants.append(compound_info)
                        elif component_counter == 1:
                            self.products.append(compound_info)
                elif lenght == 1:
                    compound_info = {
                    "stoichiometric_coefficient" : 1,
                    "substance" : splitted_section[0],
                    "rate_dependency" : 1
                    }
                    if component_counter == 0:
                        self.reactants.append(compound_info)
                    elif component_counter == 1:
                        self.products.append(compound_info)
            component_counter +=1    
        counter = 0
        
        for section in self.reactants :
            reactant = section["substance"]
            if re.match(r'^.*\.(s|g|l)$', reactant):  
                
                self.reactants[counter].update([("phase" , reactant[len(reactant)-1])])
                self.reactants[counter]["substance"] = Compound(reactant[0:len(reactant)-2])
            elif re.match(r'^.*\.aq$', reactant):
                self.reactants[counter].update([("phase" , "aq")])
                self.reactants[counter]["substance"] = Compound(reactant[0:len(reactant)-3])
            else:
                self.reactants[counter].update([("phase" , None)])
                self.reactants[counter]["substance"] = Compound(reactant)
            counter += 1
        
        counter = 0
        for section in self.products :
            product = section["substance"]
            if re.match(r'^.*\.(s|g|l)$', product):  
                self.products[counter].update([("phase" , product[len(product)-1])])
                self.products[counter]["substance"] = Compound(product[0:len(product)-2])
            elif re.match(r'^.*\.aq$', product):
                self.products[counter].update([("phase" , "aq")])
                self.products[counter]["substance"] = Compound(product[0:len(product)-3])
            else:
                self.products[counter].update([("phase" , None)])
                self.products[counter]["substance"] = Compound(product)
            counter += 1  
        self.compounds = []
        for compound in self.reactants :
            compound.update({"type" : "reactant"})
            self.compounds.append(compound)
        for compound in self.products :
            compound.update({"type" : "product"})
            self.compounds.append(compound)  
        return self
    @classmethod
    def from_simple_syntaxx(self , reaction , K=1, kf=1, kb=1):
        self.K = K
        self.kf = kf
        self.kb = kb
        self.equation = reaction
        self.reactants = []
        self.products = []
        reformed_reaction = reaction.replace(" ","").split(">")
        splited_to_component_reaction = [component.split("+") for component in reformed_reaction] 
    def __str__(self):
        return self.equation
    def __repr__(self):
        reaction_equation = ""
        counter = 0
        for compound in self.reactants :
            if int(compound["stoichiometric_coefficient"]) != 1:
                reaction_equation += ( " " + str(int(compound["stoichiometric_coefficient"])) + compound["substance"].unicode_formula) 
            else:
                reaction_equation += (" " + compound["substance"].unicode_formula)
            if compound["phase"] != None:
                reaction_equation += ( "(" + compound["phase"] + ")" )
            if counter < len(self.reactants) - 1:
                reaction_equation += " +"
            counter += 1    
        reaction_equation += " \u21cc"
        counter = 0
        for compound in self.products :
            if int(compound["stoichiometric_coefficient"]) != 1:
                reaction_equation += (" " + str(int(compound["stoichiometric_coefficient"])) + compound["substance"].unicode_formula)
            else :
                reaction_equation += (" " + compound["substance"].unicode_formula)
            if compound["phase"] != None:
                reaction_equation += ( "(" + compound["phase"] + ")" )
            if counter < len(self.products) - 1:
                reaction_equation += " +"   
            counter += 1  
        return reaction_equation
    def __add__(self , other):
        new_compounds_name = []
        new_reactants = []
        new_products = []
        for compound in (self.compounds + other.compounds ):
            compound_name = compound["substance"].formula
            if not compound_name in new_compounds_name :
                new_compounds_name.append(compound_name)
        for compound_name in new_compounds_name :
            stoichiometric_coefficient = 0
            for compound in (self.compounds + other.compounds ) :
                if compound["substance"].formula == compound_name :
                    if compound["type"] == "reactant" :
                        stoichiometric_coefficient += compound["stoichiometric_coefficient"]
                    elif compound["type"] == "product" :
                        stoichiometric_coefficient -= compound["stoichiometric_coefficient"]
            if stoichiometric_coefficient == 0:
                continue
            elif stoichiometric_coefficient > 0 :
                new_reactants.append(str(stoichiometric_coefficient) + "_" + compound_name)
            else :
                new_products.append(str(-stoichiometric_coefficient) +  "_" + compound_name)
        new_reaction = ""
        counter = 0
        for reactant in new_reactants :
            if counter < len(new_reactants) - 1:
                new_reaction += (reactant + " & ")
            else :
                new_reaction += (reactant)
            counter += 1
        new_reaction += " > "
        counter = 0
        for product in new_products :
            if counter < len(new_products) - 1:
                new_reaction += (product + " & ")
            else:
                new_reaction += (product)
            counter += 1
        try:   
            return Reaction(new_reaction)
        except:
            return None
    def __iadd__(self , other):
        return self.__add__(other)
class Enviromet():
    def _check_if_reaction(reaction):
        if isinstance(reaction , Reaction):
            return True
        else:
            raise ValueError
    def __init__(self , *reactions):
        self.reactions = []
        for reaction in reactions :
            if self._check_if_reaction(reaction):
                self.reactions.append(reaction)
    def __iadd__(self , reaction):
        if self._check_if_reaction(reaction):
            self.reactions.append(reaction)
