import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# append constants file
sys.path.append(
    os.path.dirname(
        os.path.dirname(os.path.realpath(__file__))))
from constants import *

# ignore warnings
warnings.filterwarnings('ignore')

class EquilibriumCurve:

    def __init__(self) -> None:
        """
        class for calculation of the liquid-vapor equilibrium curve
        of a binary mixture
        """

        # load physical-chemical properties table
        self.load_data()

    def load_data(self) -> None:
        """
        loads the physical and chemical dataset
        """

        self.data = pd.read_excel(PHYSICAL_CHEMICAL_SOURCE, sheet_name='vapor_pressure')

    def set_mixture_components(self, comp_list: list) -> None:
        """
        sets the components of the binary components

        Parameters
        ----------
        comp_list : list
            list of the components
        """
        self.c1 = comp_list[0]
        self.c2 = comp_list[1]

    def set_pressure(self, p: float) -> None:
        """
        sets the system pressure, in kPa

        Parameters
        ----------
        p : float
            the system's pressure
        """
        self.p = p

    def define_light_key(self) -> None:
        """
        tests the vapor pressure of both components
        and defines the light key
        """

        # the test temperature will be 100 ºC
        T = 100 + 273.15    # convert to K

        # calculate the vapor pressure of both components
        Pv1 = self.calculate_Pvap(T, self.c1)
        Pv2 = self.calculate_Pvap(T, self.c2)

        # the light key has the highest vapor pressure at the same temperature
        if Pv1 > Pv2:
            self.LK = self.c1
            self.HK = self.c2
        else:
            self.LK = self.c2
            self.HK = self.c1

    def calculate_Pvap(self, T: float, component: str) -> float:
        """
        calculates the vapor pressure of the selected
        component

        Parameters
        ----------
        T : float
            temperature of the system, in K
        component : str
            name of the component

        Returns
        -------
        float
            vapor pressure of the component, in kPa
        """

        # extract the constants of the vapor pressure correlation
        c_list = self.data.loc[self.data.component==component,:].values

        # calculate the logarithm of vapor pressure in Pa
        ln_Pvap = c_list[0][1] + (c_list[0][2] / T) + (c_list[0][3] * np.log(T)) + (c_list[0][4] * (T ** c_list[0][5]))

        # return vapor pressure in kPa
        return np.exp(ln_Pvap) / 1000

    def calculate_P_delta(self, T: float, x_lk: float) -> float:
        """
        calculates the pressure of the system using the Raoult's law
        and compares it with the real pressure 

        Parameters
        ----------
        T : float
            temperature of the system (in K)
        x_lk : float
            molar fraction of the light key

        Returns
        -------
        float
            absolute difference of the calculate pressure to real
            pressure
        """

        # calculate both components vapor pressure
        Pvap_lk = self.calculate_Pvap(T, self.LK)
        Pvap_hk = self.calculate_Pvap(T, self.HK)

        # Raoult's law
        P_calc = (x_lk * Pvap_lk) + ((1 - x_lk) * Pvap_hk)

        # absolute diference
        return abs(self.p - P_calc)

    def find_bubble_point(self, x_lk: float) -> float:
        """
        calculates the temperature of the bubble point
        as function of the molar fraction of the light key

        Parameters
        ----------
        x_lk : float
            molar fraction of th light key

        Returns
        -------
        float
            temperature of the bubble point
        """

        # initial guess of temperature is 100 ºC
        T = 100 + 373.15

        # find bubble point temperature
        t = fsolve(self.calculate_P_delta, x0=T, args=x_lk)

        return t[0]

    def calculate_equilibrium_curve(self) -> None:
        """
        calculates the vapor-liquid equilibrium curve
        """

        # create array with the molar fraction of the light-key
        # in the liquid phase
        x_lk = np.arange(0, 1.05, 0.05)

        # create equilibrium dataframe
        elv_data = pd.DataFrame(x_lk, columns=[f'x_{self.LK}'])

        # calculate molar fraction of heavy key in the liquid phase
        elv_data[f'x_{self.HK}'] = 1 - elv_data[f'x_{self.LK}']

        # calculate the temperature  of bubble point for all concentrations
        elv_data['T'] = elv_data[f'x_{self.LK}'].apply(self.find_bubble_point)

        # calcule the equilibrium constant
        elv_data[f'K_{self.LK}'] = elv_data['T'].apply(self.calculate_K, component=self.LK)
        elv_data[f'K_{self.HK}'] = elv_data['T'].apply(self.calculate_K, component=self.HK)

        # calculate the molar fraction of the components in the vapor phase
        elv_data[f'y_{self.HK}'] = elv_data[f'K_{self.HK}'] * elv_data[f'x_{self.HK}']
        elv_data[f'y_{self.LK}'] = elv_data[f'K_{self.LK}'] * elv_data[f'x_{self.LK}']

        # save equilibrium data
        elv_data.to_csv(EQUILIBRIUM_DATA, index=False)

    def calculate_K(self, T: float, component: str) -> float:
        """
        calculates the equilibrium constant for phase
        equilibria, according to Raoult's law

        Parameters
        ----------
        T : float
            temperature of the system, in K
        component : str
            name of the component

        Returns
        -------
        float
            equilibrium constant
        """

        # calculate the vapor pressure
        Pvap = self.calculate_Pvap(T=T, component=component)

        # return equilibrium constant
        return Pvap / self.p

    

if __name__ == '__main__':


    # test class creation
    elv = EquilibriumCurve()

    # sets the components of the mixture
    elv.set_mixture_components(['benzene', 'toluene'])

    # sets the system's pressure
    elv.set_pressure(p=101.325)

    # test calculation of vapor pressure
    elv.calculate_Pvap(373.15, 'benzene')

    # test definition of light key
    elv.define_light_key()

    # test Raoult's law
    elv.calculate_P_delta(T=323.15, x_lk=0.5)

    # test bubble point calculation
    elv.find_bubble_point(1)

    # test equilibrium curve calculation
    elv.calculate_equilibrium_curve()

    print(0)
