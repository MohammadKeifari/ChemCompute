import re
import math
import numpy as np
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
        return self.unicode_formula == value.unicode_formula 
        
class Reaction:
    """
    Represents a reversible chemical reaction with kinetic and equilibrium parameters.

    This class stores all information about a chemical reaction, including reactants,
    products, their stoichiometric coefficients, rate dependencies, and rate constants.
    It can be created manually or parsed from reaction strings written in either
    a simple or complex syntax.

    The class supports temperature-dependent calculations through thermodynamic properties.
    When the temperature is changed, rate constants and equilibrium constants are
    automatically updated using the Arrhenius and van't Hoff equations.

    Attributes:
        reactants (list[dict]): List of reactant dictionaries, each containing:
            - "stoichiometric_coefficient" (float)
            - "compound" (Compound)
            - "rate_dependency" (float)
        products (list[dict]): List of product dictionaries with the same structure.
        K (float): Equilibrium constant of the reaction.
        kf (float): Forward rate constant.
        kb (float): Backward rate constant.
        T (float): Reaction temperature in Kelvin. Setting this property automatically
            updates K, kf, and kb based on thermodynamic parameters.
        enthalpy (float): Enthalpy change of the reaction (J/mol). Used in van't Hoff
            equation for temperature-dependent equilibrium constant calculations.
        entropy (float): Entropy change of the reaction (J/(mol·K)).
        activation_energy_forward (float): Activation energy of the forward reaction
            (J/mol). Used in Arrhenius equation for temperature-dependent rate constant.
        activation_energy_backward (float): Activation energy of the backward reaction
            (J/mol). Used in Arrhenius equation for temperature-dependent rate constant.
        compounds (list[dict]): All involved species (reactants and products) with
            their concentration and type ("reactant" or "product").
    """
    def __init__(self,
                 reactants : list[dict] ,
                 products : list[dict] ,
                 reactants_concentration : list[float] ,
                 products_concentration : list[float] ,
                 K : float = 1,
                 enthalpy : float = 0,
                 entropy : float = 0,
                 kf : float = 1,
                 kb : float = 1,
                 activation_energy_forward : float = 0,
                 activation_energy_backward : float = 0,
                 T : float = 298):
                 
        """
        Initialize a Reaction instance.

        Args:
            reactants (list[dict]): List of reactant definitions.
            products (list[dict]): List of product definitions.
            reactants_concentration (list[float]): Initial concentrations of reactants.
            products_concentration (list[float]): Initial concentrations of products.
            K (float, optional): Equilibrium constant. Defaults to 1.
            kf (float, optional): Forward rate constant. Defaults to 1.
            kb (float, optional): Backward rate constant. Defaults to 1.
            T (float, optional): Temperature in Kelvin. Defaults to 298.
            enthalpy (float, optional): Enthalpy of the reaction. Defaults to 0.
            entropy (float, optional): Entropy of the reaction. Defaults to 0.
            activation_energy_forward (float, optional): Activation energy of the forward reaction. Defaults to 0.
            activation_energy_backward (float, optional): Activation energy of the backward reaction. Defaults to 0.
        """

        self.K = K
        self.kf = kf
        self.kb = kb
        self.reactants = reactants
        self.products = products
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.activation_energy_forward = activation_energy_forward
        self.activation_energy_backward = activation_energy_backward
        self._T = T
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
    def from_string_complex_syntax(cls, reaction_str: str,
                                   concentrations: list[float] = None,
                                   K: float = 1,
                                   enthalpy: float = 0,
                                   entropy: float = 0,
                                   kf: float = 1,
                                   kb: float = 1,
                                   activation_energy_forward: float = 0,
                                   activation_energy_backward: float = 0,
                                   T: float = 298):
                                   
        """
        Create a Reaction object from a string with complex syntax.

        The complex syntax allows compound names to include numbers and symbols
        such as + or -, and supports stoichiometric and rate order annotations.

        Format:
            "A & 2_B & ... > 3_C & 2_D_-1 & ..."
            - Prefix number = stoichiometric coefficient (default = 1)
            - Suffix number = rate dependency (default = 1)
            - Compounds may contain digits and signs (+, -, ( )).
            - Phases can be specified as .s, .l, .g, or .aq

        Example:
            "Fe(CN)6-3 & Ce+2 > Fe(CN)6-4 & Ce+3"

        Args:
            reaction_str (str): Reaction formula.
            concentrations (list[float], optional): Concentrations of reactants and products in order.
            K (float, optional): Equilibrium constant. Defaults to 1.
            kf (float, optional): Forward rate constant. Defaults to 1.
            kb (float, optional): Backward rate constant. Defaults to 1.
            T (float, optional): Temperature in Kelvin. Defaults to 298.
            enthalpy (float, optional): Enthalpy of the reaction. Defaults to 0.
            entropy (float, optional): Entropy of the reaction. Defaults to 0.
            activation_energy_forward (float, optional): Activation energy of the forward reaction. Defaults to 0.
            activation_energy_backward (float, optional): Activation energy of the backward reaction. Defaults to 0.

        Returns:
            Reaction: Parsed Reaction instance.

        Raises:
            ValueError: If the reaction string format is invalid.
        """

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
        return cls(inputed_reactants ,
                   inputed_products,
                   reactants_concentration ,
                   products_concentrations ,
                   K ,
                   enthalpy ,
                   entropy , 
                   kf , 
                   kb , 
                   activation_energy_forward , 
                   activation_energy_backward ,
                   T)
    @classmethod
    def from_string_simple_syntax(cls,
                                  reaction_str: str,
                                  concentrations: list[float] = None,
                                  K: float = 1,
                                  enthalpy: float = 0,
                                  entropy: float = 0,
                                  kf: float = 1,
                                  kb: float = 1,
                                  activation_energy_forward: float = 0,
                                  activation_energy_backward: float = 0,
                                  T: float = 298):
                                  
                                  
        """
        Create a Reaction object from a string using simple syntax.

        The simple syntax only allows alphabetic compound names (no +, -, or numbers inside names).
        It also supports optional stoichiometric and rate dependency annotations.

        Format:
            "A + 2B + ... > 3C + 2D-1 + ..."
            - Prefix number = stoichiometric coefficient (default = 1)
            - Suffix number = rate dependency (default = 1)
            - Phase can be added as .s, .l, .g, or .aq

        Example:
            "2A.g + B.g2 > C.l-1"

        Args:
            reaction_str (str): Reaction formula.
            concentrations (list[float], optional): Reactant/product concentrations.
            K (float, optional): Equilibrium constant. Defaults to 1.
            kf (float, optional): Forward rate constant. Defaults to 1.
            kb (float, optional): Backward rate constant. Defaults to 1.
            T (float, optional): Temperature in Kelvin. Defaults to 298.
            enthalpy (float, optional): Enthalpy of the reaction. Defaults to 0.
            entropy (float, optional): Entropy of the reaction. Defaults to 0.
            activation_energy_forward (float, optional): Activation energy of the forward reaction. Defaults to 0.
            activation_energy_backward (float, optional): Activation energy of the backward reaction. Defaults to 0.

        Returns:
            Reaction: Parsed Reaction instance.

        Raises:
            ValueError: If the input reaction string does not match valid format.
        """
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
        return cls(inputed_reactants ,
                   inputed_products,
                   reactants_concentration ,
                   products_concentrations ,
                   K ,
                   enthalpy ,
                   entropy , 
                   kf , 
                   kb , 
                   activation_energy_forward , 
                   activation_energy_backward ,
                   T)
    @property
    def T(self):
        """
        Get the reaction temperature.
        
        Returns:
            float: Temperature in Kelvin.
        """
        return self._T
    
    @T.setter
    def T(self , value):
        """
        Set the reaction temperature and automatically update rate constants and equilibrium constant.
        
        When the temperature is changed, the following calculations are performed:
        - Rate constants (kf, kb) are updated using the Arrhenius equation
        - Equilibrium constant (K) is updated using the van't Hoff equation
        
        The Arrhenius equation: k = k₀ * exp(-Ea/R * (1/T - 1/T₀))
        The van't Hoff equation: K = K₀ * exp(-ΔH/R * (1/T - 1/T₀))
        
        Where:
        - Ea is the activation energy (J/mol)
        - ΔH is the enthalpy change (J/mol)
        - R is the gas constant (8.3145 J/(mol·K))
        - T₀ is the previous temperature
        - T is the new temperature
        
        Args:
            value (float): New temperature in Kelvin.
            
        Note:
            This method requires that enthalpy and activation energies are set
            (non-zero values) for accurate temperature-dependent calculations.
            If these are zero, the rate constants and equilibrium constant
            will remain unchanged.
        """
        new_kf = self.kf * math.exp((-self.activation_energy_forward/8.3145) * (1/value - 1/self._T))
        new_kb = self.kb * math.exp((-self.activation_energy_backward/8.3145) * (1/value - 1/self._T))
        self.kf = new_kf
        self.kb = new_kb
        new_K = self.K * math.exp((-self.enthalpy/8.3145) * (1/value - 1/self._T))
        self.K = new_K
        self._T = value
    
    def __str__(self):
        """
        Return a human-readable chemical equation.

        Returns:
            str: Reaction equation string with phase labels.
        """
        return self.__repr__()
    def __repr__(self):
        """
        Return a formatted reversible reaction equation.

        Returns:
            str: Chemical equation formatted with Unicode ⇌ and phases.
        """
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
        """
        Combine two Reaction objects into a single net reaction.

        The resulting reaction merges reactants and products, cancelling species
        that appear on both sides.

        Args:
            other (Reaction): Another Reaction instance.

        Returns:
            Reaction: New Reaction object 
        """
        new_compounds_name = []
        new_reactants = []
        new_products = []
        concentrations = []
        for compound in (self.compounds + other.compounds ):
            compound_name = compound["compound"].formula
            if not compound_name in new_compounds_name :
                new_compounds_name.append(compound_name)
        for compound_name in new_compounds_name :
            stoichiometric_coefficient = 0
            concentration = 0
            for compound in (self.compounds + other.compounds ) :
                if compound["compound"].formula == compound_name :
                    concentration += compound["concentration"]
                    if compound["type"] == "reactant" :
                        stoichiometric_coefficient += compound["stoichiometric_coefficient"]
                       
                    elif compound["type"] == "product" :
                        stoichiometric_coefficient -= compound["stoichiometric_coefficient"]
                        
            if stoichiometric_coefficient == 0:
                continue
            elif stoichiometric_coefficient > 0 :
                new_reactants.append(str(stoichiometric_coefficient) + "_" + compound_name)
                concentrations.append(concentration)
            else :
                new_products.append(str(-stoichiometric_coefficient) +  "_" + compound_name)
                concentrations.append(concentration)
        new_reaction = ""
        counter = 0
        enthalpy = self.enthalpy + other.enthalpy
        entropy = self.entropy + other.entropy
        K = (self.K * math.exp((-self.enthalpy/8.3145) * (1/298 - 1/self.T)))* (other.K * math.exp((-other.enthalpy/8.3145) * (1/298 - 1/other.T)))
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
         
        return Reaction.from_string_complex_syntax(reaction_str =new_reaction,
                                                   concentrations = concentrations,
                                                   enthalpy = enthalpy,
                                                   entropy = entropy,
                                                   K = K,
                                                   T = 298)

    def __iadd__(self , other):
        """
        In-place addition operator for reactions.

        Equivalent to self + other.

        Args:
            other (Reaction): Another Reaction instance.

        Returns:
            Reaction: Combined reaction or None if invalid.
        """
        return self.__add__(other)
    def __iter__(self):
        """
        Iterate over all species in the reaction.

        Yields:
            dict: Compound dictionary with concentration, type, and coefficient.
        """
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
                reaction.T = T
                self.reactions.append(reaction)
        self._T = T
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
    @property
    def T(self):
        """
        Get the environment temperature.
        
        Returns:
            float: Temperature in Kelvin.
        """
        return self._T
    
    @T.setter
    def T(self , value):
        """
        Set the environment temperature and propagate to all reactions.
        
        When the environment temperature is changed, all reactions in the
        environment are updated to the new temperature. Each reaction will
        automatically recalculate its rate constants and equilibrium constant
        based on its thermodynamic parameters (enthalpy, activation energies).
        
        Args:
            value (float): New temperature in Kelvin.
        """
        self._T = value
        for reaction in self.reactions:
            reaction.T = value
            
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
            reaction.T = self.T
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
    def stoichiometric_coefficient_array(self):
        """
        Generate the stoichiometric coefficient matrix for all reactions in the environment.

        This property constructs a matrix that represents how each compound participates
        in each reaction. Each row corresponds to a reaction, and each column corresponds
        to a compound in `self.compounds`.

        - Reactants are assigned **positive** stoichiometric coefficients.
        - Products are assigned **negative** stoichiometric coefficients.

        This matrix is often used in rate law calculations, reaction network modeling,
        and dynamic simulations of multi-reaction systems.

        Returns:
            numpy.ndarray: A 2D array of shape `(n_reactions, n_compounds)` where each
            entry `[i, j]` represents the stoichiometric coefficient of compound `j`
            in reaction `i`. Positive values indicate reactants, and negative values
            indicate products.

        Example:
            Suppose an environment contains:
                Reaction 1: A + 2B ⇌ C  
                Reaction 2: C ⇌ D + E  

            And `self.compounds = [A, B, C, D, E]`.

            Then:
                >>> env.stoichiometric_coefficient_array
                array([
                    [ 1,  2, -1,  0,  0],
                    [ 0,  0,  1, -1, -1]
                ])
        """
        _stoichiometric_coefficient_array = []
        for rxn in self.reactions :
            reaction_stoichiometric_coefficients = [0] * len(self.compounds) 
            for reactant in rxn.reactants:
                index = self.compounds.index(reactant["compound"])
                reaction_stoichiometric_coefficients[index] += reactant["stoichiometric_coefficient"]
            for product in rxn.products:
                index = self.compounds.index(product["compound"])
                reaction_stoichiometric_coefficients[index] += product["stoichiometric_coefficient"] * -1
            _stoichiometric_coefficient_array.append(reaction_stoichiometric_coefficients)
        output_array = np.array(_stoichiometric_coefficient_array)
        return output_array
    
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
    def rate_constants_array(self):
        """
        Retrieve all forward and backward rate constants for reactions in the environment.

        This property aggregates the kinetic constants from each reaction object and
        returns them as a NumPy array, where each row corresponds to a reaction.

        Each row contains two values:
        - The **forward rate constant (kf)** — associated with the forward reaction direction.
        - The **backward rate constant (kb)** — associated with the reverse reaction direction.

        This structure is useful for numerical solvers and kinetic simulations where
        reaction rates are computed using vectorized operations.

        Returns:
            numpy.ndarray: A 2D array of shape `(n_reactions, 2)`, where each entry
            `[i, 0]` is the forward rate constant `kf` and `[i, 1]` is the backward
            rate constant `kb` for reaction `i`.

        Example:
            Suppose an environment contains 2 reactions:
                Reaction 1: A ⇌ B     with kf = 0.3, kb = 0.1  
                Reaction 2: B ⇌ C     with kf = 0.5, kb = 0.2  

            Then:
                >>> env.rate_constants_array
                array([
                    [0.3, 0.1],
                    [0.5, 0.2]
                ])
        """
        _rate_constants = []
        for rxn in self.reactions :
            _rate_constants.append([rxn.kf , rxn.kb])
        output_array = np.array(_rate_constants)
        return output_array
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
    def rate_dependency_array(self):
        """
        Retrieve the kinetic order (rate dependency) of each compound for all reactions.

        This property constructs a 3D NumPy array representing how the reaction rate
        depends on the concentration of each compound, for both the forward and reverse
        directions of every reaction in the environment.

        For each reaction, two vectors are generated:
        - **Reactant rate dependencies** — indicate the kinetic order of each compound
            in the forward reaction rate expression.
        - **Product rate dependencies** — indicate the kinetic order of each compound
            in the backward (reverse) reaction rate expression.

        Each vector’s length matches the total number of compounds in the environment.
        Compounds not participating in a given reaction have a dependency value of `0`.

        Returns:
            numpy.ndarray: A 3D array of shape `(n_reactions, 2, n_compounds)`, where:
                - `[:, 0, :]` corresponds to reactant rate dependencies.
                - `[:, 1, :]` corresponds to product rate dependencies.

        Example:
            Suppose the environment contains compounds [A, B, C] and one reaction:
                A + 2B2 > C
                rate_forward ∝ [A]^1 [B]^2
                rate_backward ∝ [C]^1

            Then:
                >>> env.rate_dependency_array
                array([
                    [
                        [1, 2, 0],   # Reactant dependencies (A, B, C)
                        [0, 0, 1]    # Product dependencies (A, B, C)
                    ]
                ])
        """
        _rate_dependency_array = []
        for rxn in self.reactions :
            reactants_rate_dependency =  [0] * len(self.compounds)
            for reactant in rxn.reactants:
                index = self.compounds.index(reactant["compound"])
                reactants_rate_dependency[index] = reactant["rate_dependency"]
            products_rate_dependency= [0] * len(self.compounds)
            for product in rxn.products:
                index = self.compounds.index(product["compound"])
                products_rate_dependency[index] = product["rate_dependency"]
            _rate_dependency_array.append([reactants_rate_dependency , products_rate_dependency])
        output_array = np.array(_rate_dependency_array)
        return output_array
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
    def concentrations_array(self):
        """
        Retrieve the current concentrations of all compounds in the environment.

        This property provides a NumPy array containing the concentration values
        of every compound tracked in the environment. The order of concentrations
        directly corresponds to the order of compounds in `self.compounds`.

        This representation is useful for numerical computations, matrix operations,
        and kinetic simulations where concentration vectors are required.

        Returns:
            numpy.ndarray: A 1D array of compound concentrations (in mol/L or the
            system’s chosen units), ordered consistently with `self.compounds`.

        Example:
            Suppose the environment contains:
                self.compounds = [A, B, C]
                self.compounds_concentration = [
                    {"compound": A, "concentration": 0.5},
                    {"compound": B, "concentration": 0.2},
                    {"compound": C, "concentration": 0.8}
                ]

            Then:
                >>> env.concentrations_array
                array([0.5, 0.2, 0.8])
        """
        return np.array([dict["concentration"] for dict in self.compounds_concentration])
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
            raise ValueError("The concentrations property should be a list and have the same length as the number of compounds")
        for i in range(len(self.compounds_concentration)):
            self.compounds_concentration[i]["concentration"] = value[i]
            

    def __len__(self):
        """
        Get the number of reactions currently in the environment.

        Returns:
            int: Count of reactions.
        """
        return len(self.reactions)