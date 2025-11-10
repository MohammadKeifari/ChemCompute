from ._general import Enviroment
import matplotlib
import random
from itertools import count
import numpy as np

class KineticalCalculator:
    """
    Simulates chemical reaction kinetics within an Enviroment instance.

    This class numerically integrates concentration changes of compounds over time
    using the rate constants and stoichiometric relationships defined in an `Enviroment` object.

    Attributes:
        accuracy (float): Time step for numerical integration (default: 1e-3).
        fitted (bool): Indicates whether the calculator has been linked to an `Enviroment` instance.
        enviroment (Enviroment): The fitted reaction environment (after calling `fit`).
        rate_constants (list[list[float]]): List of forward and backward rate constants for each reaction.
        reactions_by_index (list[list[list[int]]]): Index mapping of reactants and products per reaction.
        stoichiometric_coefficient_by_reaction (list[list[list[float]]]): Stoichiometric coefficients per reaction.
        rate_dependency_by_reaction (list[list[list[float]]]): Reaction rate dependencies on each compound.
        number_of_reactions (int): Number of reactions in the environment.
        concentrations (list[float]): Current concentration values for each compound in the environment.
    """
    def __init__(self , accuracy = 1e-3):
        """
        Initialize the kinetic calculator with a specified numerical accuracy.

        Args:
            accuracy (float, optional): Time step (Δt) for concentration updates.
                Smaller values yield higher accuracy but slower computation.
                Default is 1e-3.
        """
        self.accuracy = accuracy
        self.fitted = False
    def fit(self , enviroment):
        """
        Link the calculator to an existing `Enviroment` instance.

        Args:
            enviroment (Enviroment): The reaction environment containing all reactions and compounds.

        Raises:
            ValueError: If `enviroment` is not an instance of `Enviroment`.
        """
        if isinstance(enviroment, Enviroment):
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

    def calculate(self  , time , checkpoint_time = [] , plot = False , directory = "./plot.png"):
        """
        Numerically integrate the reaction kinetics over a specified time interval.

        This method simulates the time evolution of compound concentrations in the
        environment using a fixed time step defined by `self.accuracy`. It supports
        interactive plotting, saving plots to file, and recording concentration
        snapshots at specific checkpoint times.

        Args:
            time (float): Total simulation time in the same units as `self.accuracy`.
            checkpoint_time (list[float], optional): Times at which to record concentrations.
            plot (bool or str, optional): Plotting mode. Options:
                - False: Do not plot.
                - "interactive": Display real-time interactive plot.
                - "save": Save the plot to the specified `directory`.
            directory (str, optional): File path to save the plot if `plot="save"`.

        Returns:
            list[list]: List of checkpoints, where each entry is `[time, concentrations]`.
                The first entry includes the header `["time", compound_unicode_formulas]`.

        Raises:
            NameError: If the model has not been fitted to an environment (i.e., `fit` not called).
            ValueError: If an invalid plotting mode or directory is provided.
            ValueError: If `plot` is not one of [False, "save", "interactive"].

        Behavior:
            - Concentrations are clamped to zero if they become negative.
            - Supports recording concentrations at arbitrary checkpoint times.
            - Interactive plotting allows the user to type 'exit' to close the plot.

        Example:
            >>> kc = KineticalCalculator(accuracy=0.01)
            >>> kc.fit(env)
            >>> results = kc.calculate(time=10, checkpoint_time=[1,5,10], plot="interactive")
        """
        if not self.fitted :
            raise NameError("You should fit the model to an enviromt object before calculation")
        if not plot in [False , "save" , "interactive"]:
            raise ValueError("`plot` is not one of [False, 'save', 'interactive'].")
        
        if plot == "interactive" :
            matplotlib.use("TkAgg", force=True)
        elif plot == "save" :
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        
        colors = []
        checkpoints = []
        if plot != False:
            plt.figure()
            for i in self.enviroment.compounds :
                colors.append(((random.randint(0, 95)/100) , (random.randint(0, 95)/100) , (random.randint(0, 95)/100)))
                plt.xlabel("time")
                plt.ylabel("concentration")
        concentrations = self.enviroment.concentrations_array
        rate_dependencies = self.enviroment.rate_dependency_array
        stoichiometric_coefficient = self.enviroment.stoichiometric_coefficient_array
        rate_constants = self.enviroment.rate_constants_array
        time_interval = self.accuracy
        eps = 1e-300
        def calculate_rf():
            # rf = kf * prod_j concentrations[j]^order_fwd[j]
            # Use log form for numerical stability and speed
            log_c = np.log(concentrations + eps)
            log_prod = rate_dependencies[:, 0, :] @ log_c  # (R,)
            rf = np.exp(log_prod) * rate_constants[:, 0] * time_interval
            return rf
        def calculate_rb():
            log_c = np.log(concentrations + eps)
            log_prod = rate_dependencies[:, 1, :] @ log_c  # (R,)
            rb = np.exp(log_prod) * rate_constants[:, 1] * time_interval
            return rb
        def calculate_concentration_change():
            # rate vector for reactions
            rate = -(calculate_rf() - calculate_rb())  # (R,)
            # concentration change = stoich^T @ rate
            concentration_change = stoichiometric_coefficient.T @ rate  # (C,)
            return concentration_change.reshape(-1)
        t = 0
        for i in range(int(time/self.accuracy+1)):
            new_conentratinos = np.add(concentrations, calculate_concentration_change())
            new_conentratinos[new_conentratinos < 0] = 0
            if plot :
                for k in range(len(self.concentrations)):
                    plt.plot([t , t-self.accuracy],[new_conentratinos[k] , concentrations[k]] , color = colors[k])
            for checkpoint_t in checkpoint_time:
                
                if t <= checkpoint_t < t + self.accuracy:
                    checkpoints.append(new_conentratinos.copy())
            concentrations = new_conentratinos
            t += self.accuracy
        if plot == "interactive":
            for k in range(len(self.concentrations)):
                plt.plot([0 , 0],[0 , 0] , color = colors[k], label = self.enviroment.compounds[k].unicode_formula)
            plt.legend()
            plt.show(block = False)
            
            print("Type 'exit' to close the plot:")
            exited = False
            while not exited:
                cmd = input().strip().lower()
                if cmd == "exit":
                    plt.close()
                    break
                else:
                    print("Invalid input.")
        elif plot == "save" :
            for k in range(len(self.concentrations)):
                plt.plot([0 , 0],[0 , 0] , color = colors[k], label = self.enviroment.compounds[k].unicode_formula)
            plt.savefig(directory)
            plt.close('all')
            del plt
        checkpoints.append(concentrations)  
        return checkpoints

    def fit_calculate(self, enviroment, time, checkpoint_time=[], plot=False, directory="./plot.png"):
        """
        Fit the calculator to an environment and calculate reaction kinetics in one call.
        
        This is a convenience method that combines fit() and calculate() into a single
        operation. It first fits the calculator to the provided environment, then
        immediately performs the kinetic calculation.
        
        This method is equivalent to calling:
            calculator.fit(enviroment)
            results = calculator.calculate(time, checkpoint_time, plot, directory)
        
        Parameters
        ----------
        enviroment : Enviroment
            The reaction environment containing:
            - compounds: List of chemical compounds in the system
            - reactions: List of chemical reactions with rate constants
            - compounds_concentration: Initial concentrations of all compounds
            - rate_constants: Forward and backward rate constants for each reaction
            - reaction_by_index: Index mapping of reactants and products per reaction
            - stoichiometric_coefficient_by_reaction: Stoichiometric coefficients per reaction
            - rate_dependency_by_reaction: Reaction rate dependencies on each compound
        time : float
            Total simulation time in the same units as the calculator's accuracy.
        checkpoint_time : list[float], optional
            Times at which to record concentrations. Default is [].
        plot : bool or str, optional
            Plotting mode. Options:
            - False: Do not plot (default).
            - "interactive": Display real-time interactive plot.
            - "save": Save the plot to the specified `directory`.
        directory : str, optional
            File path to save the plot if `plot="save"`. Default is "./plot.png".
        
        Returns
        -------
        list[list]
            List of checkpoints, where each entry is `[time, concentrations]`.
            The first entry includes the header `["time", compound_unicode_formulas]`.
        
        Raises
        ------
        ValueError
            If `enviroment` is not an instance of `Enviroment`.
        ValueError
            If an invalid plotting mode or directory is provided.
        ValueError
            If `plot` is not one of [False, "save", "interactive"].
        
        Notes
        -----
        - This method automatically calls fit() before calculate(), so the
          environment does not need to be fitted beforehand
        - Concentrations are clamped to zero if they become negative during simulation
        - Supports recording concentrations at arbitrary checkpoint times
        - Interactive plotting allows the user to type 'exit' to close the plot
        - The method uses numerical integration with a fixed time step defined by
          the calculator's accuracy attribute
        
        Examples
        --------
        >>> from ChemCalc import Enviroment
        >>> kc = KineticalCalculator(accuracy=0.01)
        >>> env = Enviroment(...)  # Initialize environment with compounds and reactions
        >>> # Fit and calculate in one step
        >>> results = kc.fit_calculate(env, time=10, checkpoint_time=[1, 5, 10], plot="interactive")
        """
        self.fit(enviroment)
        return self.calculate(time, checkpoint_time, plot, directory)

    def calculate_responsively(self  , checkpoint_time = [] ,animation_update_interval = 0.1 ):
        """
        Simulate and visualize reaction kinetics dynamically using an interactive animation.

        This method continuously updates compound concentrations over time while
        optionally displaying a live animation using `matplotlib.animation.FuncAnimation`.
        It is especially useful for observing reaction progress in real-time.

        Args:
            checkpoint_time (list[float], optional): Specific times at which to record concentrations.
            animation_update_interval (float, optional): Time interval between animation updates (in seconds).
            plot (bool, optional): Whether to visualize the reaction dynamically.
                - If True, a live plot will be shown and user can interact with it:
                    * 'stop'  - pause the animation
                    * 'resume' - resume the animation
                    * 'exit'  - close the plot and stop the simulation
                - If False, no plot is generated.

        Returns:
            list[list]: A list of checkpoint data of the form `[time, concentrations]`.
                - Each entry is a snapshot of concentrations at the corresponding time.
                - The last entry corresponds to the final concentrations at the end of the simulation.

        Raises:
            NameError: If the model has not been fitted to an environment (i.e., `fit` not called).

        Notes:
            - Negative concentrations are clamped to zero during simulation.
            - The method relies on an internal `calculate_concentration_change` function
            to compute instantaneous changes in concentrations per reaction step.
            - This approach differs from `calculate()` by providing a live, responsive animation
            rather than a static plot or final checkpoint data.
            - Checkpoints specified in `checkpoint_time` are captured even during animation.

        Example:
            >>> kc = KineticalCalculator(accuracy=0.01)
            >>> kc.fit(env)
            >>> results = kc.calculate_responsively(checkpoint_time=[1,5,10], plot=True)
        """
        plot = True
        if not self.fitted :
            raise NameError("You must fit the model to an Enviroment before calculation.")
        if plot :
             matplotlib.use("TkAgg", force=True)
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation

        plt.figure()
        colors = []
        checkpoints = []
        if plot:
            for i in self.enviroment.compounds :
                colors.append(((random.randint(0, 95)/100) , (random.randint(0, 95)/100) , (random.randint(0, 95)/100)))
                plt.xlabel("time")
                plt.ylabel("concentration")
                
        
        concentrations = self.enviroment.concentrations_array
        rate_dependencies = self.enviroment.rate_dependency_array
        stoichiometric_coefficient = self.enviroment.stoichiometric_coefficient_array
        rate_constants = self.enviroment.rate_constants_array
        time_interval = self.accuracy
        eps = 1e-300
        def calculate_rf():
            log_c = np.log(concentrations + eps)
            log_prod = rate_dependencies[:, 0, :] @ log_c
            rf = np.exp(log_prod) * rate_constants[:, 0] * time_interval
            return rf
        def calculate_rb():
            log_c = np.log(concentrations + eps)
            log_prod = rate_dependencies[:, 1, :] @ log_c
            rb = np.exp(log_prod) * rate_constants[:, 1] * time_interval
            return rb
        def calculate_concentration_change():
            rate = -(calculate_rf() - calculate_rb())
            concentration_change = stoichiometric_coefficient.T @ rate
            return concentration_change.reshape(-1)
        
        time = count()
        def animate(i):
            """Animation update loop for real-time kinetics visualization."""

            nonlocal concentrations
            t = self.accuracy * next(time)
            
            new_conentratinos = np.add(concentrations, calculate_concentration_change())
            new_conentratinos[new_conentratinos < 0] = 0
        
            for k in range(len(self.concentrations)):
                plt.plot([t , t-self.accuracy],[new_conentratinos[k] , concentrations[k]] , color = colors[k])
                
            for checkpoint_t in checkpoint_time:
                
                if t <= checkpoint_t < t + self.accuracy:
                    checkpoints.append(new_conentratinos.copy())
            concentrations = new_conentratinos

        ani = FuncAnimation(plt.gcf() , animate , interval = animation_update_interval , cache_frame_data=False)        
        if plot :
            for k in range(len(self.concentrations)):
                plt.plot([0 , 0],[0 , 0] , color = colors[k], label = self.enviroment.compounds[k].unicode_formula)

            plt.legend()
            plt.show(block = False)
            print("Type 'exit' to close / 'stop' to pause / 'resume' to continue:")
            while True:
                cmd = input().strip().lower()
                if cmd == "exit":
                    ani.pause()
                    plt.close()
                    break
                elif cmd == "stop":
                    ani.pause()
                elif cmd == "resume":
                    ani.resume()
                else:
                    print("invalid_input")       
                    
        checkpoints.append(concentrations)  
        return checkpoints