import re
class Compound: 
    """
    Represents a chemical compound with formula, physical properties, and optional superscript/subscript formatting.

    Attributes:
        formula (str): The chemical formula of the compound.
        unicode_formula (str): Unicode representation of the formula (with sub/superscripts if enabled).
        phase_point_list (list[dict]): A list of phase data points, each as {"temperature": float, "phase": str}.
        mp (float | None): Melting point of the compound (°C or K, depending on convention).
        bp (float | None): Boiling point of the compound.
    """

    def __init__(self , formula  , phase_point_list=None , mp=None, bp=None ,scription=True):
        """
        Initialize a Compound object based on its formula, phase information, and thermal properties.

        Args:
            formula (str): The compound's chemical formula.
            phase_point_list (list[dict], optional): List of phase points with the following keys:
                - "phase" (str): One of {"s", "l", "g", "aq"}.
                - "temperature" (float): The temperature associated with that phase.
            mp (float, optional): Melting point temperature.
            bp (float, optional): Boiling point temperature.
            scription (bool, optional): If True, converts the formula into Unicode with subscripts/superscripts.
        
        Raises:
            ValueError: If a phase in `phase_point_list` is not one of {"s", "l", "g", "aq"}.
        """
        superscript_characters=["\u2070" ,"\u00b9" ,"\u00b2" ,"\u00b3" ,"\u2074" 
                                ,"\u2075" ,"\u2076" ,"\u2077" ,"\u2078" ,"\u2079" 
                                ,"\u207a" , "\u207b"]
        subscript_characters = ["\u2080" , "\u2081", "\u2082", "\u2083", "\u2084"
                                , "\u2085", "\u2086", "\u2087", "\u2088", "\u2089" ]
        phases = ["g" , "l" , "s" , "aq"]

        self.formula = formula
        if scription:
            def formula_to_unicode_formula(formula):
                """Convert a normal formula string to Unicode format with subscripts/superscripts."""
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
        """
        Determine the physical phase of the compound at a given temperature.

        Args:
            temperature (float): Temperature to evaluate phase at.

        Returns:
            str | None: One of {"s", "l", "g", "aq"} or None if phase cannot be determined.
        """
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
        """Return the Unicode representation of the compound."""
        return self.unicode_formula
    def __eq__(self, value):
        """Compare compounds based on their Unicode formulas."""
        return self.unicode_formula == value.unicode_formula :
        
class Reaction:
    """
    Represents a reversible chemical reaction, with kinetic and equilibrium parameters.

    Attributes:
        reactants (list[dict]): Reactant data (coefficient, Compound, rate dependency).
        products (list[dict]): Product data (coefficient, Compound, rate dependency).
        K (float): Equilibrium constant.
        kf (float): Forward reaction rate constant.
        kb (float): Backward reaction rate constant.
        T (float): Reaction temperature.
        compounds (list[dict]): All species involved with concentrations and type.
    """
    def __init__(self, reactants , products , reactants_concentration , products_concentration , K=1, kf=1, kb=1 , T=298):
        """
        Initialize a Reaction object.

        Args:
            reactants (list[dict]): Each dict has:
                - "stoichiometric_coefficient" (float)
                - "compound" (Compound)
                - "rate_dependency" (float)
            products (list[dict]): Same structure as reactants.
            reactants_concentration (list[float], optional): Reactant concentrations.
            products_concentration (list[float], optional): Product concentrations.
            K (float, optional): Equilibrium constant.
            kf (float, optional): Forward rate constant.
            kb (float, optional): Backward rate constant.
            T (float, optional): Temperature (K).
        """
        self.K = K
        self.kf = kf
        self.kb = kb
        self.reactants = reactants
        self.products = products
        self.T = T
        self.compounds = []
        counter = 0 
        for compound in self.reactants :
            compound.update({"concentration" : reactants_concentration[counter]})
            reactant = compound.copy()
            reactant.update({"type" : "reactant"})
            self.compounds.append(reactant)
            counter += 1
        counter = 0
        for compound in self.products :
            compound.update({"concentration" : products_concentration[counter]})
            product = compound.copy()
            product.update({"type" : "product"})
            self.compounds.append(product)
            counter += 1
    @classmethod
    def from_string_complex_syntax(cls, reaction_str , concentrations = None , K=1, kf=1, kb=1 , T = 298):
        '''
        Creates reaction object from an string with any of the following formats:
        A & 2_B & ... > 3_C & 2_D_-1 & ... \n
        The prefix number represents the stoichiometric_coefficient and the default value is 1 \n
        The suffix number represents rate_dependency or order of the reaction with respect to the reactant or product and the default vaule is 1 \n
        example: Fe(CN)6-3 & Ce+2 > Fe(CN)6-4 & Ce+3
        Split the stoichiometric_coefficient and rate_dependency from the name of the compound with "_"
        Unlike from_string_simple_syntax here your compound name could include numbers and + , -

        To specify the phase use the following syntax:
        A.g & 2B.s & ... > 3c.aq & 2D.g-1 & ...
        Args:
            reaction(string): The equation formula with the mentioned structure.

            concentrations(list of int or float): the list of concentration of reactants and products in the order they come in the reaction from left to right
            
            K(int or float): The equilibrium constant of reaction

            kf(int or float): The forward reaction rate constant

            kb(int or float): The backward reaction rate constant

        '''

        reformed_reaction = reaction_str.replace(" ","").split(">")
        splited_to_component_reaction = [component.split("&") for component in reformed_reaction] 
        inputed_reactants = []
        inputed_products = []
        component_counter = 0
        for component in splited_to_component_reaction:
            counter = 0
            acceptable_pattern_for_section = re.compile(
                r'^(\d+(?:\.\d+)?_[A-Za-z0-9+.\-()]+_-?\d+(?:\.\d+)?|' 
                r'\d+(?:\.\d+)?_[A-Za-z0-9+.\-()]+|' 
                r'[A-Za-z0-9+.\-()]+_-?\d+(?:\.\d+)?|' 
                r'[A-Za-z0-9+.\-()]+)(\.s|\.g|\.l)?$'
            )
            for section in component:

                if not bool(acceptable_pattern_for_section.match(section)):
                     raise ValueError("You can't make a reaction from string with this expression")
                
                splitted_section = section.split("_")
                lenght = len(splitted_section)
                if lenght == 3:
                    compound_info = {
                    "stoichiometric_coefficient" : float(splitted_section[0]),
                    "compound" : splitted_section[1],
                    "rate_dependency" : float(splitted_section[2])
                    }

                elif lenght == 2:
                    if re.match(r'^\d+(?:\.\d+)?_[A-Za-z0-9+.\-()]+$' , section) :
                        compound_info = {
                            "stoichiometric_coefficient" : float(splitted_section[0]),
                            "compound" : splitted_section[1],
                            "rate_dependency" : 1
                        }

                    elif re.match(r'^[A-Za-z0-9+.\-()]+_\d+(?:\.\d+)?$' , section) :
                        compound_info = {
                            "stoichiometric_coefficient" : 1,
                            "compound" : splitted_section[0],
                            "rate_dependency" : float(splitted_section[1])
                        }

                elif lenght == 1:
                    compound_info = {
                    "stoichiometric_coefficient" : 1,
                    "compound" : splitted_section[0],
                    "rate_dependency" : 1
                    }
                if component_counter == 0:
                    inputed_reactants.append(compound_info)
                elif component_counter == 1:
                    inputed_products.append(compound_info)
            component_counter +=1    
        counter = 0
        
        for section in inputed_reactants :
            reactant = section["compound"]
            if re.match(r'^.*\.(s|g|l)$', reactant):  
                inputed_reactants[counter]["compound"] = Compound(formula=reactant[0:len(reactant)-2] , phase_point_list=[{"temperature" : T , "phase" : reactant[len(reactant)-1]}])
            elif re.match(r'^.*\.aq$', reactant):
                inputed_reactants[counter]["compound"] = Compound(formula=reactant[0:len(reactant)-3] , phase_point_list=[{"temperature" : T , "phase" : "aq"}])
            else:
                inputed_reactants[counter]["compound"] = Compound(formula=reactant)
            counter += 1
        
        counter = 0
        for section in inputed_products :
            product = section["compound"]
            if re.match(r'^.*\.(s|g|l)$', product):  
                inputed_products[counter]["compound"] = Compound(formula=product[0:len(product)-2] , phase_point_list=[{"temperature" : T , "phase" : product[len(product)-1]}])
            elif re.match(r'^.*\.aq$', product):
                inputed_products[counter]["compound"] = Compound(formula=product[0:len(product)-3] , phase_point_list=[{"temperature" : T , "phase" : "aq"}])
            else:
                inputed_products[counter]["compound"] = Compound(formula=product)
            counter += 1  
        if concentrations == None:
            concentrations = [0] * (len(inputed_reactants) + len(inputed_products))
        reactants_concentration = concentrations[:len(inputed_reactants)]
        products_concentrations = concentrations[len(inputed_reactants):]
        return cls(inputed_reactants , inputed_products, reactants_concentration , products_concentrations , K , kf , kb , T)
    @classmethod
    def from_string_simple_syntax(cls , reaction_str , concentrations=None , K=1, kf=1, kb=1 , T=298):
        '''
        Creates reaction object from an string with any of the following formats:
        A + 2B + ... > 3C + 2D-1 + ... \n
        The prefix number represents the stoichiometric_coefficient and the default value is 1 \n
        The suffix number represents rate_dependency or order of the reaction with respect to the reactant or product and the default vaule is 1 \n
        B + C > CB
        Unlike from_string_complex_syntax here your compound name shouldn't include numbers and + , -

        To specify the phase use the following syntax:
        A.g + 2B.s + ... > 3c.aq + 2D.g-1 + ...

        The general order for each compound should follow the format below:
        stoichiometric_coefficient(optional) + name(only alhpabet) + .phase(.s/.g/.l/.aq , optional) + rate_dependency(optional)
        
        Args:
            reaction(string): The equation formula with the mentioned structure.

            K(int or float): The equilibrium constant of reaction

            kf(int or float): The forward reaction rate constant

            kb(int or float): The backward reaction rate constant

            T(int or float): The temperature of reaction and the constants(K, kf ,kb)
        '''
        reformed_reaction = reaction_str.replace(" ","").split(">")
        splited_to_component_reaction = [component.split("+") for component in reformed_reaction] 
        component_counter = 0
        inputed_reactants = []
        inputed_products = []
        for component in splited_to_component_reaction:
            counter = 0
            acceptable_pattern_for_section = re.compile(
                r'^(?:'
                r'\d+(?:\.\d+)?[A-Za-z]+(?:\.[A-Za-z]+)?-?\d+(?:\.\d+)?|'  
                r'\d+(?:\.\d+)?[A-Za-z]+(?:\.[A-Za-z]+)?|'                 
                r'[A-Za-z]+(?:\.[A-Za-z]+)?-?\d+(?:\.\d+)?|'               
                r'[A-Za-z]+(?:\.[A-Za-z]+)?'                         
                r')$'   
            )
            for section in component:
                if not bool(acceptable_pattern_for_section.match(section)):
                     raise ValueError("You can't make a reaction from string with this expression")
                else:
                    start_of_name_index = 0
                    end_of_name_index = 0
                    def is_number(str):
                        try:
                            float(str) 
                            return True
                        except ValueError:
                            return False
                    # number + name + number
                    if re.match(r'^\d+(?:\.\d+)?[A-Za-z]+(?:\.[A-Za-z]+)?-?\d+(?:\.\d+)?$' , section) :
                        for endpoint in range(len(section)) :
                            if is_number(section[:endpoint]) and (not is_number(section[:endpoint + 1])) :
                                start_of_name_index = endpoint 
                                break
                        for startpoint in range(start_of_name_index , len(section)):
                            if is_number(section[startpoint:]) :
                                end_of_name_index = startpoint 
                                break
                        compound_info = {
                            "stoichiometric_coefficient" : float(section[:start_of_name_index]),
                            "compound" : section[start_of_name_index:end_of_name_index],
                            "rate_dependency" : float(section[end_of_name_index:])
                            }
                    # number + name
                    elif re.match(r'^\d+(?:\.\d+)?[A-Za-z]+(?:\.[A-Za-z]+)?$' , section):
                        end_of_name_index = len(section) 
                        for endpoint in range(len(section)) :
                            if is_number(section[:endpoint]) and (not is_number(section[:endpoint + 1])) :
                                start_of_name_index = endpoint 
                                break
                        compound_info = {
                            "stoichiometric_coefficient" : float(section[:start_of_name_index]),
                            "compound" : section[start_of_name_index:end_of_name_index],
                            "rate_dependency" : 1
                            }
                    # name + number
                    elif re.match(r'^[A-Za-z]+(?:\.[A-Za-z]+)*-?\d+(?:\.\d+)?$' , section):
                        start_of_name_index = 0
                        for startpoint in range(len(section)):
                            if is_number(section[startpoint:]) :
                                end_of_name_index = startpoint 
                                break
                        compound_info = {
                            "stoichiometric_coefficient" : 1,
                            "compound" : section[start_of_name_index:end_of_name_index],
                            "rate_dependency" : float(section[end_of_name_index:])
                            }
                    # name
                    else :
                        compound_info = {
                            "stoichiometric_coefficient" : 1,
                            "compound" : section,
                            "rate_dependency" : 1
                            }
                        
                if component_counter == 0:
                    inputed_reactants.append(compound_info)
                elif component_counter == 1:
                    inputed_products.append(compound_info)
            component_counter +=1  
        counter = 0
        
        for section in inputed_reactants :
            reactant = section["compound"]
            if re.match(r'^.*\.(s|g|l)$', reactant):  
                inputed_reactants[counter]["compound"] = Compound(formula=reactant[0:len(reactant)-2] , phase_point_list=[{"temperature" : T , "phase" : reactant[len(reactant)-1]}])
            elif re.match(r'^.*\.aq$', reactant):
                inputed_reactants[counter]["compound"] = Compound(formula=reactant[0:len(reactant)-3] , phase_point_list=[{"temperature" : T , "phase" : "aq"}])
            else:
                inputed_reactants[counter]["compound"] = Compound(formula=reactant)
            counter += 1
        counter = 0
        for section in inputed_products :
            product = section["compound"]
            if re.match(r'^.*\.(s|g|l)$', product):  
                inputed_products[counter]["compound"] = Compound(formula=product[0:len(product)-2] , phase_point_list=[{"temperature" : T , "phase" : product[len(product)-1]}])
            elif re.match(r'^.*\.aq$', product):
                inputed_products[counter]["compound"] = Compound(formula=product[0:len(product)-3] , phase_point_list=[{"temperature" : T , "phase" : "aq"}])
            else:
                inputed_products[counter]["compound"] = Compound(formula=product)
            counter += 1  

        if concentrations == None:
            concentrations = [0] * (len(inputed_reactants) + len(inputed_products))
        reactants_concentration = concentrations[:len(inputed_reactants)]
        products_concentrations = concentrations[len(inputed_reactants):]
        return cls(inputed_reactants , inputed_products, reactants_concentration , products_concentrations , K , kf , kb , T)
    def __str__(self):
        """Return a concise representation of the reaction."""
        return self.__repr__()
    def __repr__(self):
        """Return a concise representation of the reaction."""
        reaction_equation = ""
        counter = 0
        for compound in self.reactants :
            if int(compound["stoichiometric_coefficient"]) != 1:
                reaction_equation += ( " " + str(int(compound["stoichiometric_coefficient"])) + compound["compound"].unicode_formula) 
            else:
                reaction_equation += (" " + compound["compound"].unicode_formula)
            if compound["compound"].phase(self.T) != None:
                reaction_equation += ( "(" + compound["compound"].phase(self.T) + ")" )
            if counter < len(self.reactants) - 1:
                reaction_equation += " +"
            counter += 1    
        reaction_equation += " \u21cc"
        counter = 0
        for compound in self.products :
            if int(compound["stoichiometric_coefficient"]) != 1:
                reaction_equation += (" " + str(int(compound["stoichiometric_coefficient"])) + compound["compound"].unicode_formula)
            else :
                reaction_equation += (" " + compound["compound"].unicode_formula)
            if compound["compound"].phase(self.T) != None:
                reaction_equation += ( "(" + compound["compound"].phase(self.T) + ")" )
            if counter < len(self.products) - 1:
                reaction_equation += " +"   
            counter += 1  
        return reaction_equation
    def __add__(self , other):
        new_compounds_name = []
        new_reactants = []
        new_products = []
        for compound in (self.compounds + other.compounds ):
            compound_name = compound["compound"].formula
            if not compound_name in new_compounds_name :
                new_compounds_name.append(compound_name)
        for compound_name in new_compounds_name :
            stoichiometric_coefficient = 0
            for compound in (self.compounds + other.compounds ) :
                if compound["compound"].formula == compound_name :
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
    def __iter__(self):
        for compound in self.compounds:
            yield compound
class Enviroment():
    """
    Represents a chemical environment containing multiple reactions and compounds.

    The `Enviroment` class acts as a container for multiple `Reaction` objects,
    automatically managing compound lists, concentration aggregation, and access
    to kinetic or stoichiometric information for simulation or analysis.

    Attributes:
        reactions (list[Reaction]): List of `Reaction` objects within the environment.
        compounds_concentration (list[dict]): List of dictionaries, each with:
            - "compound" (Compound): Compound object.
            - "concentration" (float): Current concentration value.
        compounds (list[Compound]): Unique list of all compounds appearing in any reaction.
        T (float): System temperature in Kelvin.
    """
    def _check_if_reaction(self , reaction):
        """
        Validate whether the provided object is a `Reaction` instance.

        Args:
            reaction (Reaction): Object to validate.

        Returns:
            bool: True if the object is a valid Reaction.

        Raises:
            ValueError: If `reaction` is not an instance of `Reaction`.
        """
        if isinstance(reaction , Reaction):
            return True
        else:
            raise ValueError("Only Reaction objects can be added to Enviroment.")
        
    def __init__(self , *reactions , T=298):
        """
        Initialize the environment and add one or more reactions.

        Args:
            *reactions (Reaction): Variable number of Reaction objects.
            T (float, optional): Temperature of the environment (K). Default is 298 K.

        Raises:
            ValueError: If any argument is not a Reaction object.
        """
        self.reactions = []
        for reaction in reactions :
            if self._check_if_reaction(reaction):
                self.reactions.append(reaction)
        self.T = T
        self.compounds = []
        self.compounds_concentration = []
        for reaction in reactions:
            for compound in reaction.compounds:
                compounds = [i["compound"] for i in self.compounds_concentration]
                if compound["compound"] in compounds:
                    index_in_compounds_concentration = compounds.index(compound["compound"])
                    index_in_reaction = reaction.compounds.index(compound)
                    self.compounds_concentration[index_in_compounds_concentration]["concentration"] += reaction.compounds[index_in_reaction]["concentration"]
                else:
                    index_in_reaction = reaction.compounds.index(compound)
                    self.compounds_concentration.append({"compound" : compound["compound"] , "concentration" :reaction.compounds[index_in_reaction]["concentration"]})
                    self.compounds.append(compound["compound"])
    def __iadd__(self , reaction):
        """
        Add a reaction to the environment using the += operator.

        Args:
            reaction (Reaction): Reaction to add.

        Returns:
            Enviroment: The updated environment instance.

        Raises:
            ValueError: If `reaction` is not a valid Reaction object.
        """
        if self._check_if_reaction(reaction):
            self.reactions.append(reaction)
            self.compounds = []
            self.compounds_concentration = []
            for reaction in self.reactions:
                for compound in reaction.compounds:
                    compounds = [i["compound"] for i in self.compounds_concentration]
                    if compound["compound"] in compounds:
                        index_in_compounds_concentration = compounds.index(compound["compound"])
                        index_in_reaction = reaction.compounds.index(compound)
                        self.compounds_concentration[index_in_compounds_concentration]["concentration"] += reaction.compounds[index_in_reaction]["concentration"]
                    else:
                        index_in_reaction = reaction.compounds.index(compound)
                        self.compounds_concentration.append({"compound" : compound["compound"] , "concentration" :reaction.compounds[index_in_reaction]["concentration"]})
                        self.compounds.append(compound["compound"])
            return self
    def __iter__(self):
        """
        Iterate through all reactions in the environment.

        Yields:
            Reaction: Each reaction in the environment.
        """
        for reaction in self.reactions:
            yield reaction
    def add(self , reaction):
        """
        Add a new reaction to the environment manually.

        Args:
            reaction (Reaction): The reaction to add.

        Raises:
            ValueError: If `reaction` is not a valid Reaction object.
        """
        if self._check_if_reaction(reaction):
            self.reactions.append(reaction)
            self.compounds = []
            self.compounds_concentration = []
            for reaction in self.reactions:
                for compound in reaction.compounds:
                    compounds = [i["compound"] for i in self.compounds_concentration]
                    if compound["compound"] in compounds:
                        index_in_compounds_concentration = compounds.index(compound["compound"])
                        index_in_reaction = reaction.compounds.index(compound)
                        self.compounds_concentration[index_in_compounds_concentration]["concentration"] += reaction.compounds[index_in_reaction]["concentration"]
                    else:
                        index_in_reaction = reaction.compounds.index(compound)
                        self.compounds_concentration.append({"compound" : compound["compound"] , "concentration" :reaction.compounds[index_in_reaction]["concentration"]})
                        self.compounds.append(compound["compound"])
    @property
    def reaction_by_index(self):
        """
        Map each reaction’s reactants and products to their indices in the environment’s compound list.

        Returns:
            list[list[list[int]]]: A list of [reactants_index, products_index] for each reaction.
        """
        _reactions_by_index = []
        for rxn in self.reactions :
            reatants_index = []
            for reactant in rxn.reactants:
                index = self.compounds.index(reactant["compound"])
                reatants_index.append(index)
            products_index = []
            for product in rxn.products:
                index = self.compounds.index(product["compound"])
                products_index.append(index)
            _reactions_by_index.append([reatants_index , products_index])
        return _reactions_by_index
    @property
    def stoichiometric_coefficient_by_reaction(self):
        """
        Return stoichiometric coefficients for all reactions.

        Returns:
            list[list[list[float]]]: A list of [reactant_coefficients, product_coefficients] per reaction.
        """
        _stoichiometric_coefficient_by_reaction = []
        for rxn in self.reactions :
            reatants_index = []
            for reactant in rxn.reactants:
                reatants_index.append(reactant["stoichiometric_coefficient"])
            products_index = []
            for product in rxn.products:
                products_index.append(product["stoichiometric_coefficient"])
            _stoichiometric_coefficient_by_reaction.append([reatants_index , products_index])
        return _stoichiometric_coefficient_by_reaction
    @property
    def rate_constants(self):
        """
        Get all forward and backward rate constants for each reaction.

        Returns:
            list[list[float]]: Each entry is [kf, kb] for a reaction.
        """
        _rate_constants = []
        for rxn in self.reactions :
            _rate_constants.append([rxn.kf , rxn.kb])
        return _rate_constants
    @property
    def rate_dependency_by_reaction(self):
        """
        Return the kinetic order (rate dependency) for reactants and products in each reaction.

        Returns:
            list[list[list[float]]]: Each entry contains:
                [ [reactant_rate_dependencies], [product_rate_dependencies] ]
        """
        _rate_dependency_by_reaction = []
        for rxn in self.reactions :
            reatants_index = []
            for reactant in rxn.reactants:
                reatants_index.append(reactant["rate_dependency"])
            products_index = []
            for product in rxn.products:
                products_index.append(product["rate_dependency"])
            _rate_dependency_by_reaction.append([reatants_index , products_index])
        return _rate_dependency_by_reaction
    @property
    def compounds_unicode_formula(self):
        """
        Get the Unicode formulas of all compounds in the environment.

        Returns:
            list[str]: List of compound Unicode formula strings.
        """
        return [compound.unicode_formula for compound in self.compounds]

    @property
    def concentrations(self):
        """
        Retrieve the current concentrations of all compounds in the environment.

        Returns:
            list[float]: List of compound concentrations in the same order as `self.compounds`.
        """
        return [dict["concentration"] for dict in self.compounds_concentration]

    @concentrations.setter
    def concentrations(self , value):
        """
        Update the concentration values for all compounds.

        Args:
            value (list[float]): New concentration values corresponding to each compound.

        Raises:
            ValueError: If input is not a list or its length doesn’t match compound count.
        """
        if not(type(value) == list and len(value) == len(self.compounds_concentration)):
            raise ValueError("The concentrations property should be a list and half the same lenght as the number of compounds")
        for i in range(len(self.compounds_concentration)):
            self.compounds_concentration[i]["concentration"] = value[i]
            
    def change_temperature():
        pass

    def __len__(self):
        """
        Get the number of reactions currently in the environment.

        Returns:
            int: Count of reactions.
        """
        return len(self.reactions)