import os
import sys
import warnings
import numpy as np
import pandas as pd

# append constants file
sys.path.append(
    os.path.dirname(
        os.path.dirname(os.path.realpath(__file__))))
from constants import *

# append additional modules
sys.path.append(THERMODYNAMICS_PATH)
from thermodynamics.equilibrium import EquilibriumCurve

# ignore warnings
warnings.filterwarnings('ignore')

class McCabeThiele:

    # equilibrium curve
    elv_module = EquilibriumCurve()

    def __init__(self) -> None:
        """
        class for the implementation of the McCabe-Thiele
        method
        """

        # load equilibrium curve
        self.load_vle_curve()

    def load_vle_curve(self) -> None:
        """
        loads the vapor liquid equilibrium curve
        """
        # loads the vapor liquid equilibrium curve
        self.vle = pd.read_csv(EQUILIBRIUM_DATA)

    def set_feed_temperature(self, t: float) -> None:
        """
        sets the temperature of the feed stream

        Parameters
        ----------
        t : float
            temperature of the feed stream [K]
        """
        # saves feed temperature
        self.T_f = t

    def set_column_pressure(self, p: float) -> None:
        """
        sets the pressure of the column

        Parameters
        ----------
        p : float
            pressure of the column [kPa]
        """

        # sets the pressure for the equilibrium module
        self.elv_module.set_pressure(p)

    def set_mixture_components(self, components: list) -> None:
        """
        sets the binary mixture components 

        Parameters
        ----------
        components : list
            name of the components
        """
        # set components to elv module
        self.elv_module.set_mixture_components(components)

        # define light key
        self.elv_module.define_light_key()

    def set_mixture_composition(self, zF: float) -> None:
        """
        sets the composition of the binary mixture

        Parameters
        ----------
        zF : float
            molar fraction of the light-key component
        """
        # set feed mixture molar composition
        self.z_LK = zF
        self.z_HK = 1 - zF

    def calculate_boiling_point(self) -> None:
        """
        calculates the boiling point temperature of the feed
        """

        # calculate the bubble point temperature of the 
        # feed mixture
        self.Tb = self.elv_module.find_bubble_point(self.z_LK)

    def calculate_feed_enthalpy(self) -> None:
        """
        calculates the enthalpy of the feed stream, in [kJ/mol]
        """
        # this version only considers the liquid phase
        if self.T_f <= self.Tb:
            self.HF = (Cp / 1000) * (self.T_f - Tref)

    def calculate_q_line(self) -> None:
        """
        calculates the q-line of the feed
        """

        # the enthalpy of the vapor phase is ignores the sensible
        # heat
        self.q = (Hvap - self.HF) / (Hvap)

        # calculates the slope and the intercept
        self.m_q = (self.q) / (self.q - 1)
        self.b_q = (self.z_LK) / (self.q - 1)

        # create dataframe to store the q-line coordinates
        self.q_data = pd.DataFrame(np.arange(0.2, self.z_LK + 0.01, 0.01), columns=['x_q'])
        self.q_data['y_q'] = self.q_data['x_q'].apply(lambda x: self.m_q * x - self.b_q)

    def set_reflux_ratio(self, r: float) -> None:
        """
        sets the reflux ratio of the column

        Parameters
        ----------
        r : float
            reflux ratio
        """
        self.r = r

    def set_LK_compositions(self, design: list) -> None:
        """
        sets the design composition of the light key component
        in the distillate and in the bottom product

        Parameters
        ----------
        design : list
            list with the molar composition of the light
            key in the distillate and in the bottom product
        """

        # set distillate and bottom product molar composition
        self.xD = design[0]
        self.xB = design[1]

    def calculate_operational_lines(self) -> None:
        """
        calculates the operational lines of enrichment and
        of stripping sections
        """

        # calculate slope and intercept of the enrichment
        # operation line
        self.m_e = self.r / (self.r + 1)
        self.b_e = self.xD / (self.r + 1)

        # create dataframe with the coordinates of the enrichment operation line
        self.el_data = pd.DataFrame(np.arange(self.z_LK, self.xD + 0.01, 0.01), columns=['x_el'])
        self.el_data['y_el'] = self.el_data['x_el'].apply(lambda x: self.m_e * x + self.b_e)

        # calculate slope and intercept of the stripping
        # operation line
        self.m_s = (self.m_e * self.z_LK + self.b_e - self.xB) / (self.z_LK - self.xB)
        self.b_s = self.xB - (self.m_s * self.xB)

        # create dataframe with the coordinates of the stripping operation line
        self.sl_data = pd.DataFrame(np.arange(self.xB, self.z_LK + 0.01, 0.01), columns=['x_sl'])
        self.sl_data['y_sl'] = self.sl_data['x_sl'].apply(lambda x: self.m_s * x + self.b_s)

    def find_interpolation_pair(self, x: float) -> list:
        """
        finds the interpolation pair necessary to 
        create the "touch point" between operation curve and 
        equilibrium curve

        Parameters
        ----------
        x : float
            point to be interpolated

        Returns
        -------
        list
            coordinates of the interpolation
            pair
        """

        # iterate between y values and check if the desired values
        # is between them
        for i in range(self.vle.shape[0] - 1):
            ymin = self.vle.loc[i, f'y_{self.elv_module.LK}']
            ymax = self.vle.loc[i+1, f'y_{self.elv_module.LK}']

            if ((ymin < x) & (ymax > x)):
                
                xmin = self.vle.loc[i, f'x_{self.elv_module.LK}']
                xmax = self.vle.loc[i+1, f'x_{self.elv_module.LK}']

                return [[xmin, ymin], [xmax, ymax]]

        return None

    @staticmethod
    def interpolate(interpolation_pair: list, y: float) -> float:
        """
        interpolates the equilibrium curve to find next
        composition of liquid phase

        Parameters
        ----------
        interpolation_pair : list
            interpolation pair of the equilibrium curve
        y : float
            vapor molar fraction to be interpolated

        Returns
        -------
        float
            molar fraction of liquid phase 
            of the next stage
        """
        
        # extract interpolation points
        xmin = interpolation_pair[0][0]
        ymin = interpolation_pair[0][1]
        xmax = interpolation_pair[1][0]
        ymax = interpolation_pair[1][1]

        # interpolate point
        x = xmin + (((y - ymin) * (xmax - xmin)) / (ymax - ymin))

        return x

    def calculate_equilibrium_stages(self) -> None:
        """
        applies the McCabe-Thiele method to calculate
        the equilibrium stages
        """

        # create dataframe to save equilibrium stages
        self.eq_stages = pd.DataFrame([[self.xD, self.xD]], columns=['x', 'y'], index=[0])

        xi = self.eq_stages.tail(1).x.values[0]
        while xi >= self.xB:

            # interpolate to find vapor composition in equilibrium
            pair = self.find_interpolation_pair(x = xi)
            y_eq = self.interpolate(pair, xi)

            # find liquid phase of next stage using the right
            # operational line
            if y_eq > self.z_LK:

                # use enrichment line
                y_next = self.m_e * y_eq + self.b_e

            else:
                # use stripping line
                y_next = self.m_s * y_eq + self.b_s

            # add the two new points to dataset
            self.eq_stages.loc[self.eq_stages.shape[0], :] = y_eq, xi
            self.eq_stages.loc[self.eq_stages.shape[0], :] = y_eq, y_next

            # update xi
            xi = y_next

    def set_initial_flow_rate(self, f: float) -> None:
        """
        sets the feed stream molar flow rate, in [kmol/h]

        Parameters
        ----------
        f : float
            feed stream molar flow rate, in [kmol/h]
        """
        self.F = f

    def calculate_streams(self) -> None:
        """
        calculates the molar flow rates of
        the distillate and of the bottom product
        by solving the mass balances
        """

        # set linear algebraic system
        A = np.array([[1, 1], [self.xD, self.xB]])
        B = np.array([self.F, self.z_LK * self.F])

        # solve linear system
        result = np.linalg.solve(A, B)

        # set parameters
        self.B = result[1]
        self.D = result[0]

    def make_column_design(self) -> None:
        """
        wraps the calculations up in one method
        """

        # calculate mixture boiling point
        self.calculate_boiling_point()

        # calculate enthalpy of the feed stream
        self.calculate_feed_enthalpy()

        # calculate q-line
        self.calculate_q_line()

        # calculate operation lines
        self.calculate_operational_lines()

        # calculate equilibrium stages
        self.calculate_equilibrium_stages()

        # calculate molar flow rates
        self.calculate_streams()


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    # test class creation
    mccabe = McCabeThiele()

    # set feed information and design data
    mccabe.set_feed_temperature(327.6)
    mccabe.set_column_pressure(101.3)
    mccabe.set_mixture_components(['benzene', 'toluene'])
    mccabe.set_mixture_composition(zF=0.45)
    mccabe.set_reflux_ratio(r=4)
    mccabe.set_LK_compositions([0.95, 0.10])
    mccabe.set_initial_flow_rate(100)

    # makes column design
    mccabe.make_column_design()

    # test figure creation
    plt.figure(figsize=(10,10))
    plt.plot(mccabe.vle.x_benzene, mccabe.vle.x_benzene, 'k-')
    plt.plot(mccabe.vle.x_benzene, mccabe.vle.y_benzene, 'ro--', label='Equilibrium Curve')
    plt.plot(mccabe.q_data.x_q, mccabe.q_data.y_q, 'g-', label='Q-Line')
    plt.plot(mccabe.el_data.x_el, mccabe.el_data.y_el, 'b-', label='Enrichment Line')
    plt.plot(mccabe.sl_data.x_sl, mccabe.sl_data.y_sl, 'b-', label='Stripping Line')
    plt.plot(mccabe.eq_stages.x, mccabe.eq_stages.y, color='orange', ls='-' , label='Equilibrium Stages')
    plt.xlabel('Molar Fraction Benzene - Liquid', size=14)
    plt.ylabel('Molar Fraction Benzene - Vapor', size=14)
    plt.legend(loc='best')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.axvline(x=mccabe.xD, color='gray', ls='--')
    plt.axvline(x=mccabe.xB, color='gray', ls='--')
    plt.grid(True, alpha=0.3)
    plt.show()
