import re
import math
import warnings
import numpy as np

POURBAIX_RESERVED_SPECIES = frozenset({"H+", "OH-"})


class XS:
    """
    Excess-species marker for concentration dictionaries.

    ``XS(amount)`` sets concentration and marks the species as excess (fixed activity).
    Bare ``XS()`` keeps the current amount and only marks excess (environment overrides).
    """

    __slots__ = ("amount",)

    def __init__(self, amount=None):
        self.amount = None if amount is None else float(amount)

    def __repr__(self) -> str:
        if self.amount is None:
            return "XS()"
        return f"XS({self.amount})"


def _is_xs(value) -> bool:
    """Return True when a concentration override marks excess."""
    return isinstance(value, XS) or (isinstance(value, str) and value.upper() == "XS")


def _species_formula(entry) -> str:
    compound = entry["compound"]
    return compound.formula if hasattr(compound, "formula") else str(compound)


def _resolve_concentration_spec(value, *, current: float = 0.0) -> tuple[float, bool]:
    """Parse a concentration dict value into ``(concentration, excess)``."""
    if _is_xs(value):
        if isinstance(value, XS) and value.amount is not None:
            return value.amount, True
        return current, True
    return float(value), False


def _apply_concentration_map(entries: list[dict], concentration_map: dict) -> None:
    """Apply formula->concentration mapping; missing species default to 0."""
    for entry in entries:
        formula = _species_formula(entry)
        spec = concentration_map.get(formula, 0.0)
        concentration, excess = _resolve_concentration_spec(spec, current=0.0)
        entry["concentration"] = concentration
        entry["excess"] = excess


def _apply_concentration_list(entries: list[dict], values: list) -> None:
    """Apply ordered concentration values; each slot may be numeric or :class:`XS`."""
    for entry, spec in zip(entries, values):
        concentration, excess = _resolve_concentration_spec(spec, current=0.0)
        entry["concentration"] = concentration
        entry["excess"] = excess


def _filter_pourbaix_reserved_concentrations(concentrations, *, T=298):
    """Drop reserved pH species from user-supplied concentration maps."""
    from ._mixing import _resolve_compound_key

    filtered = {}
    for key, value in concentrations.items():
        if isinstance(key, str):
            formula = key
        else:
            formula = _resolve_compound_key(key, T).formula
        if formula in POURBAIX_RESERVED_SPECIES:
            warnings.warn(
                f"Ignoring assigned concentration for reserved species {formula!r}; "
                "Pourbaix and constant-pH workflows control H+ and OH-.",
                stacklevel=3,
            )
            continue
        filtered[key] = value
    return filtered
class Compound: 
    """
    Represents a chemical compound with formula, physical properties, and optional superscript/subscript formatting.

    Attributes:
        formula (str): The chemical formula of the compound.
        unicode_formula (str): Unicode representation of the formula (with sub/superscripts if enabled).
        phase_point_list (list[dict]): A list of phase data points, each as {"temperature": float, "phase": str}.
        mp (float | None): Melting point of the compound (°C or K, depending on convention).
        bp (float | None): Boiling point of the compound.
        spectrum: Optional UV-Vis molar absorptivity spec (:class:`SpectrumSpec`).
        token (str): ``@c{id}`` slot for interpolating this object into ``from_string``.
    """

    def __init__(self , formula  , phase_point_list=None , mp=None, bp=None ,scription=True, charge=0, spectrum=None):
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
            charge (int, optional): Ionic charge for activity-coefficient calculations. Default 0.
            spectrum: Optional :class:`SpectrumSpec` (or compatible) molar absorptivity curve.
        
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
        self.charge = int(charge)
        self.spectrum = spectrum

    def set_spectrum(self, spectrum_spec):
        """Attach a UV-Vis molar absorptivity specification to this compound."""
        self.spectrum = spectrum_spec

    @property
    def token(self) -> str:
        """Register this compound and return an ``@c{id}`` slot for ``from_string``."""
        from ._interpolation import register_compound

        return register_compound(self)

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

    def __hash__(self):
        return hash(self.formula)


def _ensure_species_rate_defaults(species_list: list[dict]) -> None:
    """Fill missing stoichiometric coefficient and rate dependency (order = stoich)."""
    for species in species_list:
        stoich = float(species.get("stoichiometric_coefficient", 1))
        species.setdefault("stoichiometric_coefficient", stoich)
        species.setdefault("rate_dependency", stoich)


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
            - "rate_dependency" (float; defaults to the stoichiometric coefficient)
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
            their concentration, excess flag, and type ("reactant" or "product").
    """
    def __init__(self,
                 reactants : list[dict] ,
                 products : list[dict] ,
                 reactants_concentration : list[float] = None,
                 products_concentration : list[float] = None,
                 concentrations : dict = None,
                 K : float = 1,
                 enthalpy : float = 0,
                 entropy : float = 0,
                 kf : float = 1,
                 kb : float = 1,
                 activation_energy_forward : float = 0,
                 activation_energy_backward : float = 0,
                 T : float = 298,
                 infinite_K : bool = False):
                 
        """
        Initialize a Reaction instance.

        Args:
            reactants (list[dict]): List of reactant definitions. Missing
                ``rate_dependency`` defaults to the stoichiometric coefficient.
            products (list[dict]): List of product definitions. Missing
                ``rate_dependency`` defaults to the stoichiometric coefficient.
            reactants_concentration (list, optional): Initial reactant concentrations
                (legacy list form, one value per reactant in order). Values may be
                numeric or :class:`XS`.
            products_concentration (list, optional): Initial product concentrations
                (legacy list form, one value per product in order). Values may be
                numeric or :class:`XS`.
            concentrations (dict, optional): ``{formula: amount}`` for all species in the
                reaction. Missing species default to ``0``. Values may be numeric or
                :class:`XS` to mark excess (``XS(amount)`` sets amount and excess;
                bare ``XS()`` is for environment-style keep-current overrides only).
            K (float, optional): Equilibrium constant. Defaults to 1.
            kf (float, optional): Forward rate constant. Defaults to 1.
            kb (float, optional): Backward rate constant. Defaults to 1.
            T (float, optional): Temperature in Kelvin. Defaults to 298.
            enthalpy (float, optional): Enthalpy of the reaction. Defaults to 0.
            entropy (float, optional): Entropy of the reaction. Defaults to 0.
            activation_energy_forward (float, optional): Activation energy of the forward reaction. Defaults to 0.
            activation_energy_backward (float, optional): Activation energy of the backward reaction. Defaults to 0.
            infinite_K (bool, optional): If True, treat as irreversible (K effectively infinite).
                The reaction is driven to completion during equilibrium calculations.
        """

        self.K = K
        self.kf = kf
        self.kb = kb
        self.infinite_K = infinite_K
        self.reactants = reactants
        self.products = products
        _ensure_species_rate_defaults(self.reactants)
        _ensure_species_rate_defaults(self.products)
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.activation_energy_forward = activation_energy_forward
        self.activation_energy_backward = activation_energy_backward
        self._T = T
        self._K_ref = K
        self._T_ref = T
        self._adjust_thermodynamics = True
        self.rate_law = "mass_action"
        self.bio_params = {}
        self.compounds = []
        self._assign_species_concentrations(
            reactants_concentration,
            products_concentration,
            concentrations,
        )

    def _assign_species_concentrations(
        self,
        reactants_concentration,
        products_concentration,
        concentrations,
    ) -> None:
        from ._half_reaction import _reject_electron_formula

        for compound in self.reactants:
            spec = compound.get("compound")
            if isinstance(spec, str):
                _reject_electron_formula(spec)
            elif hasattr(spec, "formula"):
                _reject_electron_formula(spec.formula)
        for compound in self.products:
            spec = compound.get("compound")
            if isinstance(spec, str):
                _reject_electron_formula(spec)
            elif hasattr(spec, "formula"):
                _reject_electron_formula(spec.formula)

        if concentrations is not None:
            _apply_concentration_map(self.reactants, concentrations)
            _apply_concentration_map(self.products, concentrations)
        else:
            if reactants_concentration is None or products_concentration is None:
                raise ValueError(
                    "Provide concentrations={formula: amount} or both "
                    "reactants_concentration and products_concentration lists."
                )
            _apply_concentration_list(self.reactants, reactants_concentration)
            _apply_concentration_list(self.products, products_concentration)

        self.compounds = []
        for compound in self.reactants:
            reactant = compound.copy()
            reactant.setdefault("excess", False)
            reactant.update({"type": "reactant"})
            self.compounds.append(reactant)
        for compound in self.products:
            product = compound.copy()
            product.setdefault("excess", False)
            product.update({"type": "product"})
            self.compounds.append(product)
    @classmethod
    def from_string(cls, reaction_str: str,
                                   concentrations=None,
                                   K: float = 1,
                                   enthalpy: float = 0,
                                   entropy: float = 0,
                                   kf: float = 1,
                                   kb: float = 1,
                                   activation_energy_forward: float = 0,
                                   activation_energy_backward: float = 0,
                                   T: float = 298,
                                   infinite_K: bool = False):
                                   
        """
        Create a Reaction from a string.

        Accepts the same thermodynamic and kinetic keyword arguments as
        :meth:`__init__`, including ``infinite_K``.

        ``concentrations`` may be a list (legacy, reactants then products in order)
        or a ``{formula: amount}`` dict. Missing dict entries default to ``0``.
        Use :class:`XS` values to mark excess species (e.g. ``{"H2O": XS(55.5)}``
        or ``[XS(55.5), 1e-7, 1e-7]``).

        Format:
            "A & 2_B & ... > 3_C & 2_D_-1 & ..."
            - ``&`` separates species on each side; ``>`` separates reactants/products
            - Prefix ``n_`` = stoichiometric coefficient (default 1)
            - Suffix ``_n`` = rate order (defaults to the stoichiometric coefficient)
            - Phases: ``.s``, ``.l``, ``.g``, ``.aq``
            - Ionic charge inferred from trailing ``+`` / ``-`` in species names
            - Live objects: ``f"{water.token} > H+ & OH-"`` (see :attr:`Compound.token`)

        Example:
            "Fe(CN)6-3 & Ce+2 > Fe(CN)6-4 & Ce+3"
        """
        from ._formula import split_phase_suffix
        from ._interpolation import compound_from_parsed_name

        reformed_reaction = reaction_str.replace(" ","").split(">")
        splited_to_component_reaction = [component.split("&") for component in reformed_reaction] 
        inputed_reactants = []
        inputed_products = []
        component_counter = 0
        _species = r"(?:@c\d+|[A-Za-z0-9+.\-()]+)"
        for component in splited_to_component_reaction:
            counter = 0
            acceptable_pattern_for_section = re.compile(
                rf"^(\d+(?:\.\d+)?_{_species}_-?\d+(?:\.\d+)?|"
                rf"\d+(?:\.\d+)?_{_species}|"
                rf"{_species}_-?\d+(?:\.\d+)?|"
                rf"{_species})(\.s|\.g|\.l|\.aq)?$"
            )
            for section in component:

                if not bool(acceptable_pattern_for_section.match(section)):
                     raise ValueError("You can't make a reaction from string with this expression")

                core, phase = split_phase_suffix(section)
                phase_suffix = f".{phase}" if phase else ""
                splitted_section = core.split("_")
                lenght = len(splitted_section)
                if lenght == 3:
                    compound_info = {
                    "stoichiometric_coefficient" : float(splitted_section[0]),
                    "compound" : splitted_section[1] + phase_suffix,
                    "rate_dependency" : float(splitted_section[2])
                    }

                elif lenght == 2:
                    if re.match(r'^\d+(?:\.\d+)?$', splitted_section[0]):
                        stoich = float(splitted_section[0])
                        compound_info = {
                            "stoichiometric_coefficient" : stoich,
                            "compound" : splitted_section[1] + phase_suffix,
                            "rate_dependency" : stoich
                        }

                    else:
                        compound_info = {
                            "stoichiometric_coefficient" : 1,
                            "compound" : splitted_section[0] + phase_suffix,
                            "rate_dependency" : float(splitted_section[1])
                        }

                elif lenght == 1:
                    compound_info = {
                    "stoichiometric_coefficient" : 1,
                    "compound" : splitted_section[0] + phase_suffix,
                    "rate_dependency" : 1
                    }
                if component_counter == 0:
                    inputed_reactants.append(compound_info)
                elif component_counter == 1:
                    inputed_products.append(compound_info)
            component_counter +=1    
        for index, section in enumerate(inputed_reactants):
            inputed_reactants[index]["compound"] = compound_from_parsed_name(
                section["compound"], T=T
            )

        for index, section in enumerate(inputed_products):
            inputed_products[index]["compound"] = compound_from_parsed_name(
                section["compound"], T=T
            )
        if concentrations is None:
            concentrations = {}
        elif isinstance(concentrations, (list, tuple)):
            reactants_concentration = list(concentrations[: len(inputed_reactants)])
            products_concentrations = list(concentrations[len(inputed_reactants) :])
            return cls(
                inputed_reactants,
                inputed_products,
                reactants_concentration,
                products_concentrations,
                K=K,
                enthalpy=enthalpy,
                entropy=entropy,
                kf=kf,
                kb=kb,
                activation_energy_forward=activation_energy_forward,
                activation_energy_backward=activation_energy_backward,
                T=T,
                infinite_K=infinite_K,
            )
        if isinstance(concentrations, dict):
            return cls(
                inputed_reactants,
                inputed_products,
                concentrations=concentrations,
                K=K,
                enthalpy=enthalpy,
                entropy=entropy,
                kf=kf,
                kb=kb,
                activation_energy_forward=activation_energy_forward,
                activation_energy_backward=activation_energy_backward,
                T=T,
                infinite_K=infinite_K,
            )
        raise TypeError("concentrations must be a dict, list, or tuple.")
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
        if not getattr(self, "_adjust_thermodynamics", True):
            self._T = value
            return
        new_kf = self.kf * math.exp((-self.activation_energy_forward/8.3145) * (1/value - 1/self._T))
        new_kb = self.kb * math.exp((-self.activation_energy_backward/8.3145) * (1/value - 1/self._T))
        self.kf = new_kf
        self.kb = new_kb
        delta_g_ref = self.enthalpy - self._T * self.entropy
        delta_g_new = self.enthalpy - value * self.entropy
        if self.entropy != 0:
            self.K = self._K_ref * math.exp(
                -(delta_g_new / (8.3145 * value) - delta_g_ref / (8.3145 * self._T_ref))
            )
        else:
            self.K = self.K * math.exp((-self.enthalpy/8.3145) * (1/value - 1/self._T))
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
         
        return Reaction.from_string(reaction_str=new_reaction,
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

    def equilibrium(
        self,
        *,
        method: str = "newton",
        loss: str = "log_quotient",
        max_iter=None,
        learning_rate=None,
        tol=None,
        backtrack_beta: float = 0.5,
        min_concentration: float = 1e-12,
        quotient_error_limit=None,
        huber_delta: float = 1.0,
        return_details: bool = False,
    ):
        """
        Calculate equilibrium for this single reaction without manually building an Enviroment.
        """
        env = Enviroment(
            self,
            T=self.T,
            volume=1.0,
            adjust_thermodynamics=getattr(self, "_adjust_thermodynamics", True),
        )
        result = env.equilibrium(
            method=method,
            loss=loss,
            max_iter=max_iter,
            learning_rate=learning_rate,
            tol=tol,
            backtrack_beta=backtrack_beta,
            min_concentration=min_concentration,
            quotient_error_limit=quotient_error_limit,
            huber_delta=huber_delta,
            return_details=True,
        )
        self._last_equilibrium_result = result
        if return_details:
            return result
        return result.concentrations

    def kinetics(
        self,
        time,
        checkpoint_time=None,
        plot=False,
        directory="./plot.png",
        colors=None,
        *,
        accuracy: float = 1e-3,
    ):
        """Integrate kinetics for this single reaction without manually building an Enviroment."""
        env = Enviroment(self, T=self.T, volume=1.0)
        return env.kinetics(
            time,
            checkpoint_time=checkpoint_time,
            plot=plot,
            directory=directory,
            colors=colors,
            accuracy=accuracy,
        )

    @property
    def last_equilibrium_result(self):
        return getattr(self, "_last_equilibrium_result", None)
    
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
            - "excess" (bool): If True, that species' amount is held fixed in equilibrium.
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

    def _add_half_reaction_item(self, half_reaction):
        from ._half_reaction import HalfReaction, register_half_reactions

        if not isinstance(half_reaction, HalfReaction):
            raise ValueError("Expected a HalfReaction instance.")
        half_reaction.attach_or_create(self)
        register_half_reactions(self, [half_reaction])
        return half_reaction
        
    def __init__(
        self,
        *items,
        T=298,
        adjust_thermodynamics=True,
        activity_model=None,
        concentrations=None,
        volume=1.0,
        buffer=None,
        half_reactions=None,
        electrode_Eh=None,
    ):
        """
        Initialize the environment and add reactions and/or half-reactions.

        Args:
            *items (Reaction | HalfReaction): Reactions and half-reactions in any order.
            T (float, optional): Temperature of the environment (K). Default is 298 K.
            adjust_thermodynamics (bool, optional): If True, update K/kf/kb when T changes.
                If False, K/kf/kb remain constant (default True).
            activity_model (str, ActivityModel, bool, or None): Ionic activity model for
                equilibrium Q. None disables activity corrections.
            concentrations (dict, optional): Override or set species concentrations by
                formula string or Compound key. Overrides values summed from reactions.
                Use :data:`XS` to mark a species as excess without changing its amount,
                or call :meth:`set_excess` after construction. Excess is stored on the
                environment next to concentration, not on :class:`Compound`.
            volume (float, optional): Solution volume in litres. Default 1.0.
            buffer (list, set, or dict, optional): Species held at fixed concentration during
                equilibrium and kinetics. Use ``buffer=["H+"]`` with ``concentrations={"H+": ...}``
                for constant pH. Dict values set explicit targets; list entries snap from current
                concentrations after build.
            half_reactions (list[HalfReaction], optional): Additional half-reactions to register.
            electrode_Eh (float, optional): Fixed electrode potential vs SHE (V). When set,
                redox equilibrium constants follow the Nernst equation at this potential.

        Raises:
            ValueError: If any positional item is not a Reaction or HalfReaction.
        """
        from ._activity import normalize_activity_model
        from ._half_reaction import HalfReaction, register_half_reactions

        if volume <= 0:
            raise ValueError("Environment volume must be positive.")

        self.reactions = []
        self.half_reactions = []
        self._electrode_Eh = electrode_Eh
        self._T = T
        self.adjust_thermodynamics = adjust_thermodynamics
        self.charge_map = {}
        self._activity_model = normalize_activity_model(activity_model)
        self.volume = float(volume)
        for item in items:
            if isinstance(item, Reaction):
                item._adjust_thermodynamics = adjust_thermodynamics
                item.T = T
                self.reactions.append(item)
            elif isinstance(item, HalfReaction):
                item.T = T
                self._add_half_reaction_item(item)
            else:
                raise ValueError(
                    "Enviroment items must be Reaction or HalfReaction instances."
                )
        if half_reactions:
            for hr in half_reactions:
                hr.T = T
                hr.attach_or_create(self)
            register_half_reactions(self, half_reactions)
        self.compounds = []
        self.compounds_concentration = []
        self._rebuild_compound_list()
        if concentrations:
            self._apply_concentration_overrides(
                _filter_pourbaix_reserved_concentrations(concentrations, T=T)
            )
        self._buffer_spec = self._normalize_buffer_spec(buffer)
        self._resolve_buffer_targets()
        self._last_equilibrium_result = None

    @classmethod
    def from_compounds(
        cls,
        concentrations,
        *,
        T=298,
        volume=1.0,
        adjust_thermodynamics=True,
        activity_model=None,
        buffer=None,
    ):
        """
        Create an environment from compounds and concentrations with no reactions.

        Args:
            concentrations (dict): Mapping of Compound or formula str to concentration (mol/L).
                Values may be numeric or :class:`XS` to mark excess.
            buffer (list, set, or dict, optional): Species held at fixed concentration. See
                :meth:`__init__` for constant-pH usage with ``buffer=["H+"]``.
        """
        from ._activity import normalize_activity_model
        from ._mixing import _resolve_compound_key

        if volume <= 0:
            raise ValueError("Environment volume must be positive.")

        env = cls.__new__(cls)
        env.reactions = []
        env.half_reactions = []
        env._electrode_Eh = None
        env.adjust_thermodynamics = adjust_thermodynamics
        env.charge_map = {}
        env._activity_model = normalize_activity_model(activity_model)
        env._T = T
        env.volume = float(volume)
        env.compounds = []
        env.compounds_concentration = []
        for key, concentration in concentrations.items():
            compound = _resolve_compound_key(key, T)
            amount, excess = _resolve_concentration_spec(concentration, current=0.0)
            env.compounds.append(compound)
            env.compounds_concentration.append(
                {"compound": compound, "concentration": amount, "excess": excess}
            )
        env._buffer_spec = env._normalize_buffer_spec(buffer)
        env._resolve_buffer_targets()
        env._last_equilibrium_result = None
        return env

    def _normalize_buffer_spec(self, buffer):
        """Convert buffer input to {formula: explicit_target_or_None}."""
        from ._mixing import _resolve_compound_key

        if buffer is None:
            return {}
        if isinstance(buffer, (list, tuple, set)):
            spec = {}
            for key in buffer:
                compound = _resolve_compound_key(key, self.T)
                spec[compound.formula] = None
            return spec
        if isinstance(buffer, dict):
            spec = {}
            for key, value in buffer.items():
                compound = _resolve_compound_key(key, self.T)
                spec[compound.formula] = float(value) if value is not None else None
            return spec
        raise TypeError(
            "buffer must be a list, set, tuple, or dict of species keys, or None."
        )

    def _resolve_buffer_targets(self):
        """Resolve fixed concentrations for buffered species from spec and current state."""
        from ._mixing import _resolve_compound_key

        targets = {}
        for formula, explicit in self._buffer_spec.items():
            if formula in self.compound_labels:
                index = self.compound_labels.index(formula)
                if explicit is not None:
                    targets[formula] = explicit
                    self.compounds_concentration[index]["concentration"] = explicit
                else:
                    targets[formula] = self.concentrations[index]
            elif explicit is not None:
                compound = _resolve_compound_key(formula, self.T)
                self.compounds.append(compound)
                self.compounds_concentration.append(
                    {"compound": compound, "concentration": explicit, "excess": False}
                )
                targets[formula] = explicit
            else:
                raise ValueError(
                    f"Buffered species {formula!r} is not in the environment; "
                    "provide an explicit concentration in buffer={{...}}."
                )
        self._buffer_targets = targets

    def set_buffer(self, buffer):
        """
        Replace buffered species and re-resolve fixed concentrations.

        For constant pH, set ``[H+]`` first (via ``concentrations=`` or direct edit),
        then call ``set_buffer(["H+"])`` to snapshot the current value.

        Args:
            buffer: Same forms as the ``buffer`` parameter on :meth:`__init__`.
        """
        self._buffer_spec = self._normalize_buffer_spec(buffer)
        self._resolve_buffer_targets()

    @property
    def buffer_targets(self) -> dict[str, float]:
        """Fixed concentrations for buffered species (formula -> mol/L)."""
        return dict(getattr(self, "_buffer_targets", {}))

    @property
    def buffer_indices(self):
        """Compound indices aligned with :attr:`compounds` for buffered species."""
        labels = self.compound_labels
        return [labels.index(formula) for formula in self._buffer_targets]

    def _set_concentration(
        self,
        key,
        value,
        *,
        excess: bool = False,
        allow_pourbaix_reserved: bool = False,
    ) -> None:
        """Apply one concentration override; ``XS()`` keeps amount and marks excess."""
        from ._mixing import _resolve_compound_key

        compound = _resolve_compound_key(key, self.T)
        formula = compound.formula
        if formula in POURBAIX_RESERVED_SPECIES and not allow_pourbaix_reserved:
            raise ValueError(
                f"{formula!r} is reserved for Pourbaix / constant-pH workflows. "
                "Set it with set_buffer(['H+']) and apply_pourbaix_state(...), "
                "or pass allow_pourbaix_reserved=True internally."
            )

        mark_excess = excess or _is_xs(value)
        current = 0.0
        if formula in self.compound_labels:
            current = self.compounds_concentration[self.compound_labels.index(formula)]["concentration"]
        numeric_value, parsed_excess = _resolve_concentration_spec(value, current=current)
        mark_excess = mark_excess or parsed_excess
        if not _is_xs(value) or (isinstance(value, XS) and value.amount is not None):
            resolved_value = numeric_value
        else:
            resolved_value = None

        if formula in self.compound_labels:
            index = self.compound_labels.index(formula)
            if resolved_value is not None:
                self.compounds_concentration[index]["concentration"] = resolved_value
            elif not mark_excess:
                raise ValueError(
                    f"Concentration override for {formula!r} must be numeric or XS."
                )
            if mark_excess:
                self.compounds_concentration[index]["excess"] = True
        else:
            if _is_xs(value) and (not isinstance(value, XS) or value.amount is None):
                raise ValueError(
                    f"Cannot use XS() for {formula!r}; species is not in the environment."
                )
            self.compounds.append(compound)
            self.compounds_concentration.append(
                {
                    "compound": compound,
                    "concentration": resolved_value if resolved_value is not None else 0.0,
                    "excess": mark_excess,
                }
            )

    def _apply_concentration_overrides(self, concentrations, *, allow_pourbaix_reserved=False):
        """Apply concentration dict; overrides reaction-derived values."""
        for key, value in concentrations.items():
            self._set_concentration(
                key,
                value,
                allow_pourbaix_reserved=allow_pourbaix_reserved,
            )

    def set_excess(self, concentrations) -> None:
        """
        Set concentrations and mark species as excess (fixed activity).

        Values may be numeric (mol/L) or :data:`XS` to keep the current amount and
        only mark the environment concentration entry as excess.

        Example::

            env.set_excess({"H2O": XS, "CaF2": 10.0})
        """
        for key, value in concentrations.items():
            self._set_concentration(key, value, excess=True)

    @classmethod
    def combine(cls, *terms):
        """Combine environments with optional (coefficient, env) terms."""
        from ._mixing import combine_environments

        return combine_environments(*terms)

    def add_compounds(self, concentrations, *, volume=1.0, coefficient=1.0):
        """Return a new environment with an added concentration slug mixed in."""
        from ._mixing import add_compounds_to_environment

        return add_compounds_to_environment(
            self,
            concentrations,
            volume=volume,
            coefficient=coefficient,
        )

    def __add__(self, other):
        from ._mixing import ScaledEnviroment

        if isinstance(other, ScaledEnviroment):
            return self.combine((1.0, self), other)
        if isinstance(other, Enviroment):
            return self.combine((1.0, self), (1.0, other))
        return NotImplemented

    def __mul__(self, coefficient):
        from ._mixing import ScaledEnviroment

        return ScaledEnviroment(coefficient, self)

    def __rmul__(self, coefficient):
        from ._mixing import ScaledEnviroment

        return ScaledEnviroment(coefficient, self)
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
            reaction._adjust_thermodynamics = self.adjust_thermodynamics
            reaction.T = value

    @property
    def activity_model(self):
        """Current ActivityModel instance, or None for ideal solution."""
        return self._activity_model

    @activity_model.setter
    def activity_model(self, value):
        from ._activity import normalize_activity_model

        self._activity_model = normalize_activity_model(value)

    def set_spectrum(self, formula: str, spectrum_spec):
        """Attach a SpectrumSpec to the compound(s) in this environment with ``formula``."""
        found = False
        for compound in self.compounds:
            if compound.formula == formula:
                compound.spectrum = spectrum_spec
                found = True
        for reaction in self.reactions:
            for entry in reaction.reactants + reaction.products:
                species = entry["compound"]
                if getattr(species, "formula", None) == formula:
                    species.spectrum = spectrum_spec
                    found = True
        if not found:
            raise ValueError(f"Compound {formula!r} is not in the environment.")

    def copy(self):
        """Return a deep copy of this environment for titration and other workflows."""
        import copy as copy_module

        from ._mixing import rewire_half_reactions, rewire_reaction_compounds

        formula_to_compound = {
            compound.formula: copy_module.deepcopy(compound) for compound in self.compounds
        }

        new_env = Enviroment.__new__(Enviroment)
        new_env.adjust_thermodynamics = self.adjust_thermodynamics
        new_env.charge_map = dict(self.charge_map)
        new_env._activity_model = self._activity_model
        new_env._T = self._T
        new_env.volume = self.volume
        new_env.reactions = copy_module.deepcopy(self.reactions)
        rewire_reaction_compounds(new_env.reactions, formula_to_compound)
        new_env.compounds_concentration = [
            {
                "compound": formula_to_compound[entry["compound"].formula],
                "concentration": float(entry["concentration"]),
                "excess": bool(entry.get("excess", False)),
            }
            for entry in self.compounds_concentration
        ]
        new_env.compounds = [entry["compound"] for entry in new_env.compounds_concentration]
        new_env._buffer_spec = dict(getattr(self, "_buffer_spec", {}))
        new_env._buffer_targets = dict(getattr(self, "_buffer_targets", {}))
        new_env._last_equilibrium_result = None
        new_env.half_reactions = copy_module.deepcopy(getattr(self, "half_reactions", []))
        rewire_half_reactions(new_env.half_reactions, formula_to_compound)
        new_env._electrode_Eh = getattr(self, "_electrode_Eh", None)
        for hr in new_env.half_reactions:
            if hr._reaction_index is not None and hr._reaction_index < len(new_env.reactions):
                pass
            else:
                hr._reaction_index = None
                hr.attach_or_create(new_env)
        for reaction in new_env.reactions:
            reaction._adjust_thermodynamics = new_env.adjust_thermodynamics
        return new_env

    @property
    def electrode_Eh(self):
        """User-imposed electrode potential vs SHE (V), or None for coupled solve."""
        return getattr(self, "_electrode_Eh", None)

    def set_electrode_potential(self, Eh: float):
        """Fix electrode potential during equilibrium (Pourbaix / potentiostat)."""
        self._electrode_Eh = float(Eh)

    def clear_electrode_potential(self):
        """Clear imposed electrode potential and restore linked reaction K values."""
        from ._half_reaction import apply_electrode_potential

        self._electrode_Eh = None
        for hr in getattr(self, "half_reactions", None) or []:
            idx = hr._reaction_index
            if idx is None or idx >= len(self.reactions):
                continue
            rxn = self.reactions[idx]
            if hasattr(rxn, "_K_ref"):
                rxn.K = rxn._K_ref

    def register_half_reactions(self, half_reactions):
        """Register additional half-reactions on this environment."""
        from ._half_reaction import register_half_reactions

        for hr in half_reactions:
            hr.attach_or_create(self)
        register_half_reactions(self, half_reactions)
        self._rebuild_compound_list()

    def electrode_potential(self, half_reaction=None):
        """Return Nernst E (V vs SHE) for one or all half-reactions."""
        hrs = getattr(self, "half_reactions", None) or []
        if not hrs:
            return None
        if half_reaction is not None:
            return half_reaction.E_at(self)
        values = [hr.E_at(self) for hr in hrs]
        return float(sum(values) / len(values))

    def buffer_diagnostics(self, equilibrium_concentrations=None):
        """Compute buffer capacity and Henderson-Hasselbalch diagnostics."""
        from ._buffer import buffer_diagnostics

        return buffer_diagnostics(self, equilibrium_concentrations)

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
            reaction._adjust_thermodynamics = self.adjust_thermodynamics
            reaction.T = self.T
            self.reactions.append(reaction)
            self._rebuild_compound_list()
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
            reaction._adjust_thermodynamics = self.adjust_thermodynamics
            self.reactions.append(reaction)
            self._rebuild_compound_list()

    def _rebuild_compound_list(self):
        """Rebuild aggregated compound list from all reactions."""
        self.compounds = []
        self.compounds_concentration = []
        for reaction in self.reactions:
            for compound in reaction.compounds:
                compounds = [i["compound"] for i in self.compounds_concentration]
                index_in_reaction = reaction.compounds.index(compound)
                entry_excess = compound.get("excess", False)
                if compound["compound"] in compounds:
                    index_in_compounds_concentration = compounds.index(compound["compound"])
                    self.compounds_concentration[index_in_compounds_concentration]["concentration"] += reaction.compounds[index_in_reaction]["concentration"]
                    if entry_excess:
                        self.compounds_concentration[index_in_compounds_concentration]["excess"] = True
                    kept = self.compounds_concentration[index_in_compounds_concentration]["compound"]
                    incoming = compound["compound"]
                    if getattr(kept, "spectrum", None) is None and getattr(incoming, "spectrum", None) is not None:
                        kept.spectrum = incoming.spectrum
                else:
                    self.compounds_concentration.append(
                        {
                            "compound": compound["compound"],
                            "concentration": reaction.compounds[index_in_reaction]["concentration"],
                            "excess": bool(entry_excess),
                        }
                    )
                    self.compounds.append(compound["compound"])

        for entry in self.compounds_concentration:
            compound = entry["compound"]
            if compound.charge != 0 and compound.formula not in self.charge_map:
                self.charge_map[compound.formula] = compound.charge

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

    @property
    def compound_labels(self):
        """
        Ordered compound formula labels aligned with ``concentrations``.

        Returns
        -------
        list[str]
            Chemical formulas in the same order as ``self.compounds``.
        """
        return [compound.formula for compound in self.compounds]

    @property
    def concentrations_dict(self) -> dict[str, float]:
        """
        Map compound formula labels to current concentrations.

        Returns
        -------
        dict[str, float]
            ``{formula: concentration}`` for every compound in the environment.
        """
        return dict(zip(self.compound_labels, self.concentrations))

    @property
    def excess_dict(self) -> dict[str, bool]:
        """Map formula labels to whether that species is marked excess in this environment."""
        return {
            label: bool(entry.get("excess", False))
            for label, entry in zip(self.compound_labels, self.compounds_concentration)
        }

    @property
    def excess_indices(self):
        """Compound indices whose environment concentration is marked excess."""
        return [
            index
            for index, entry in enumerate(self.compounds_concentration)
            if entry.get("excess", False)
        ]

    def equilibrium(
        self,
        *,
        method: str = "newton",
        loss: str = "log_quotient",
        max_iter=None,
        learning_rate=None,
        tol=None,
        backtrack_beta: float = 0.5,
        min_concentration: float = 1e-12,
        quotient_error_limit=None,
        huber_delta: float = 1.0,
        return_details: bool = False,
    ):
        """
        Calculate equilibrium concentrations for this environment.

        Parameters
        ----------
        method : str, optional
            Optimization method: ``"bgd"``, ``"sgd"``, or ``"newton"``. Default ``"newton"``.
        loss : str, optional
            Loss function: ``"log_quotient"``, ``"quotient_error"``, or ``"log_huber"``.
            Default ``"log_quotient"``.
        max_iter : int, optional
            Maximum iterations. Defaults depend on ``method``.
        learning_rate : float, optional
            Step size. Defaults depend on ``method``.
        tol : float, optional
            Residual convergence tolerance. Defaults depend on ``method``.
            Ignored when ``quotient_error_limit`` is set.
        backtrack_beta : float, optional
            Backtracking line search factor. Default ``0.5``.
        min_concentration : float, optional
            Floor for log computations. Default ``1e-12``.
        quotient_error_limit : float, optional
            Stop when every reaction satisfies ``|Q/K - 1| <= limit``.
            When set, ``tol`` is ignored. For example, ``0.01`` means within 1% of K.
        huber_delta : float, optional
            Delta parameter for the ``"log_huber"`` loss. Default ``1.0``.
        return_details : bool, optional
            If ``True``, return an :class:`EquilibriumResult` with per-reaction
            diagnostics. Default ``False``. Regardless of this flag, the full
            result is stored on ``last_equilibrium_result``. The result includes
            ``criterion_met``, which checks the final solution against
            ``quotient_error_limit`` (if set) or ``tol`` (otherwise).

        Returns
        -------
        list[float] or EquilibriumResult
            Equilibrium concentrations aligned with ``self.compounds``, or a full
            result object when ``return_details=True``.
        """
        from ._equilibrium import solve_equilibrium, EquilibriumResult

        if len(self.reactions) == 0:
            result = EquilibriumResult(
                concentrations=list(self.concentrations),
                compounds=list(self.compound_labels),
                reaction_extents=[],
                reaction_quotient_error=[],
                max_reaction_quotient_error=0.0,
                reaction_quotient_ratio=[],
                converged=True,
                stop_reason="no_reactions",
                iterations=0,
                criterion_met=True,
                criterion_type="no_reactions",
                criterion_value=0.0,
                criterion_limit=0.0,
                electrode_Eh=getattr(self, "electrode_Eh", None),
            )
            self._last_equilibrium_result = result
            if return_details:
                return result
            return result.concentrations

        result = solve_equilibrium(
            self,
            method=method,
            loss=loss,
            max_iter=max_iter,
            learning_rate=learning_rate,
            tol=tol,
            backtrack_beta=backtrack_beta,
            min_concentration=min_concentration,
            quotient_error_limit=quotient_error_limit,
            huber_delta=huber_delta,
            return_details=True,
        )
        self._last_equilibrium_result = result
        if return_details:
            return result
        return result.concentrations

    @property
    def last_equilibrium_result(self):
        """
        Most recent :class:`EquilibriumResult` from ``equilibrium()`` or
        ``apply_equilibrium()``.

        Returns
        -------
        EquilibriumResult or None
            None if no equilibrium calculation has been run yet.
        """
        return getattr(self, "_last_equilibrium_result", None)

    def apply_equilibrium(
        self,
        *,
        method: str = "newton",
        loss: str = "log_quotient",
        max_iter=None,
        learning_rate=None,
        tol=None,
        backtrack_beta: float = 0.5,
        min_concentration: float = 1e-12,
        quotient_error_limit=None,
        huber_delta: float = 1.0,
    ):
        """
        Calculate equilibrium and write concentrations back to this environment.

        Parameters
        ----------
        method, loss, max_iter, learning_rate, tol, backtrack_beta,
        min_concentration, quotient_error_limit, huber_delta
            Same as :meth:`equilibrium`.

        Returns
        -------
        EquilibriumResult
            Full result including ``concentrations``, ``q_over_k``,
            ``stop_reason``, and ``iterations``.
        """
        result = self.equilibrium(
            method=method,
            loss=loss,
            max_iter=max_iter,
            learning_rate=learning_rate,
            tol=tol,
            backtrack_beta=backtrack_beta,
            min_concentration=min_concentration,
            quotient_error_limit=quotient_error_limit,
            huber_delta=huber_delta,
            return_details=True,
        )
        self.concentrations = result.concentrations
        return result

    def kinetics(
        self,
        time,
        checkpoint_time=None,
        plot=False,
        directory="./plot.png",
        colors=None,
        *,
        accuracy: float = 1e-3,
    ):
        """
        Integrate reaction kinetics over time for this environment.

        Parameters
        ----------
        time : float
            Total simulation time.
        checkpoint_time : list[float], optional
            Times at which to record concentrations.
        plot : bool or str, optional
            Plotting mode: ``False``, ``"interactive"``, or ``"save"``.
        directory : str, optional
            File path when ``plot="save"``.
        colors : list, optional
            Plot colors, one per compound.
        accuracy : float, optional
            Integration time step. Default ``1e-3``.

        Returns
        -------
        list
            Checkpoint concentration snapshots.
        """
        from ._kinetics import integrate_kinetics
        from ._bio_kinetics import integrate_bio_kinetics, uses_bio_kinetics
        import warnings

        if getattr(self, "half_reactions", None):
            warnings.warn(
                "Half-reactions and electrode potential are not applied during kinetics; "
                "concentrations evolve with static reaction K values.",
                stacklevel=2,
            )

        if checkpoint_time is None:
            checkpoint_time = []

        integrator = integrate_bio_kinetics if uses_bio_kinetics(self) else integrate_kinetics
        return integrator(
            self,
            time=time,
            accuracy=accuracy,
            checkpoint_time=checkpoint_time,
            plot=plot,
            directory=directory,
            colors=colors,
        )

    def __len__(self):
        """
        Get the number of reactions currently in the environment.

        Returns:
            int: Count of reactions.
        """
        return len(self.reactions)