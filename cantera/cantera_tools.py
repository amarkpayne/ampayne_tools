from cantera import ck2cti
import cantera as ct
import numpy as np
from copy import deepcopy
from collections import OrderedDict
import bokeh.plotting as bp
import bokeh.palettes as bpp
import bokeh.models as bm


def run_ck2cti(input_file, thermo_file, output_file):
    parser = ck2cti.Parser()
    parser.convertMech(input_file, thermoFile=thermo_file, outName=output_file, quiet=True, permissive=True)
    print 'Done'


class CanteraSimulation:

    def __init__(self, mechfile, intial_mol_fracs, temp_0, pressure_0, energy='off', atol=1e-23, rtol=1e-8):
        self.mechfile = mechfile
        self.initial_mol_fracs = intial_mol_fracs
        self.temp_0 = temp_0
        self.pressure_0 = pressure_0
        self.solution = None
        self.reactor_net = None
        self.energy = energy
        self.atol = atol
        self.rtol = rtol

        self.reintialize_simulation()
        self.identifier_lookup = {}
        for i, spc in enumerate(self.solution.species()):
            self.identifier_lookup[spc.name] = i

        self.rxn_identifier_lookup = {}
        for i, rxn in enumerate(self.solution.reactions()):
            self.rxn_identifier_lookup[rxn.equation] = i

    def reintialize_simulation(self, setup_sensitivities=False, thermo_sensitivity=False):
        if not self.solution:
            self.solution = ct.Solution(infile=self.mechfile)

        self.solution.TPX = self.temp_0, self.pressure_0, self.initial_mol_fracs
        cantera_reactor = ct.IdealGasConstPressureReactor(contents=self.solution, energy=self.energy)
        self.reactor_net = ct.ReactorNet([cantera_reactor])
        self.reactor_net.atol = self.atol
        self.reactor_net.rtol = self.rtol

        if setup_sensitivities:
            for m in range(len(self.solution.reactions())):
                cantera_reactor.add_sensitivity_reaction(m)
                
            if thermo_sensitivity:
                for m in range(len(self.solution.species())):
                    cantera_reactor.add_sensitivity_species_enthalpy(m)
            # self.reactor_net.atol_sensitivity = 1e-23
            # self.reactor_net.rtol_sensitivity = 1e-8

    def find_equilibrium(self, constrainted_state_vars, intitial_conditions=None):
        """

        :param constrainted_state_vars: One of the following: 'TP','TV','HP','SP','SV','UV'
        :param intitial_conditions: tuple of T0, P0, mole fraction dictionary (defaults __init__ values)
        :return: None
        """
        if intitial_conditions:
            self.solution.TPX = intitial_conditions[0], intitial_conditions[1], intitial_conditions[2]
        else:
            self.solution.TPX = self.temp_0, self.pressure_0, self.initial_mol_fracs

        self.solution.equilibrate(constrainted_state_vars)

    def find_t50(self, specie, criteria='mass_frac'):
        """

        :param specie: Cantera species name (string)
        :param criteria: 'mol_frac' (50% reduction in X) or 'mass_frac' (50% reduction in Y)
        :return: Time [sec] of 50% conversion from initial amount
        """

        self.reintialize_simulation()
        spc_index = self.identifier_lookup[specie]
        if criteria == 'mol_frac':
            x0 = self.solution.X[spc_index]
        elif criteria == 'mass_frac':
            x0 = self.solution.mass_fraction_dict()[specie]
        else:
            raise ValueError('Invalid criteria: {0}'.format(criteria))

        while True:
            self.reactor_net.step()

            if criteria == 'mol_frac':
                x1 = self.solution.X[spc_index]
            else:
                x1 = self.solution.mass_fraction_dict()[specie]

            if x1 < x0*0.5:
                return self.reactor_net.time

    def return_highest_concentrated_species(self, n):
        highest_spc_indices = deepcopy(self.solution.X).argsort()[-n:][::-1]
        concentration_dict = OrderedDict()

        for index in highest_spc_indices:
            concentration_dict[self.solution.species()[index].name] = self.solution.X[index]

        return concentration_dict

    def get_concentration_profiles(self, final_time, spc_list, n_time_points=100):
        self.reintialize_simulation()
        time_points = np.linspace(0, final_time, n_time_points)
        spc_indices = np.array(map(lambda x: self.identifier_lookup[x], spc_list))
        spc_concentrations = np.zeros((spc_indices.shape[0], time_points.shape[0]))

        for i, t in enumerate(time_points):
            self.reactor_net.advance(t)
            spc_concentrations[:, i] = np.log10(self.solution.X.take(spc_indices))

        return time_points, spc_concentrations

    @staticmethod
    def plot_concentration_profiles(times, profiles, spc_list):
        fig = bp.figure(x_axis_label='Time [sec]', y_axis_label='Log10(X)')
        color_palette = bpp.Accent6
        legend_items = []

        for i in range(profiles.shape[0]):
            legend_items += [(spc_list[i], [fig.line(times, profiles[i, :], color=color_palette[i], line_width=5)])]

        legend = bm.Legend(items=legend_items, location=(0, -30))

        fig.add_layout(legend, 'right')

        bp.show(fig)

    def make_temperature_time_plot(self, tf, n_time_points=100):
        self.reintialize_simulation()
        temps = np.zeros(n_time_points)
        times = np.linspace(0, tf, n_time_points)

        for i, t in enumerate(times):
            self.reactor_net.advance(t)
            temps[i] = self.solution.T

        fig = bp.figure(x_axis_label='Time [sec]', y_axis_label='Temperature [K]')
        color_palette = bpp.Accent6
        fig.line(times, temps, color=color_palette[0], line_width=5)

        bp.show(fig)

    def get_rops(self, time_list):
        self.reintialize_simulation()
        stoichiometry_matrix = self.solution.product_stoich_coeffs() - self.solution.reactant_stoich_coeffs()
        rop_matrix = np.zeros((len(self.solution.species()), len(self.solution.reactions()), len(time_list)))
        for i, t in enumerate(time_list):
            self.reactor_net.advance(t)
            rop_matrix[:, :, i] = np.dot(stoichiometry_matrix, np.diag(self.solution.net_rates_of_progress))

        return rop_matrix

    def plot_highest_rops(self, specie, tf, n=5, n_time_points=100):
        time_list = np.linspace(0, tf, n_time_points)
        time_list = time_list[1:]
        rop_matrix = self.get_rops(time_list)
        rxn_set = set()
        spc_index = self.identifier_lookup[specie]
        color_palette = bpp.Paired[12]

        for i, _ in enumerate(time_list):
            rxn_dict = self.get_highest_rop_rxns(n, rop_matrix, specie, i)

            for label in rxn_dict.keys():
                rxn_set.add(label)

        fig = bp.figure(x_axis_label='Time [sec]', y_axis_label='Net ROP for {0} [kmol/m3-s]'.format(specie), plot_width=800,
                        plot_height=500)
        legend_items = []
        j = 0
        for r_label in rxn_set:
            legend_items += [(r_label,
                              [fig.line(time_list,
                                        rop_matrix[spc_index, self.rxn_identifier_lookup[r_label], :],
                                        line_width=5, color=color_palette[j])
                               ]
                              )
                             ]
            j += 1
        legend = bm.Legend(items=legend_items, location=(0, -30))
        fig.add_layout(legend, 'right')

        bp.show(fig)

    def get_highest_rop_rxns(self, n, rop_matrix, specie_label, time_index):
        spc_index = self.identifier_lookup[specie_label]
        rop_for_specie = abs(deepcopy(rop_matrix[spc_index, :, time_index]))
        highest_rxn_indices = rop_for_specie.argsort()[-n:][::-1]

        rxn_dict = OrderedDict()
        for rxn_index in highest_rxn_indices:
            rxn_dict[self.solution.reactions()[rxn_index].equation] = rop_matrix[spc_index, rxn_index]

        return rxn_dict

    def run_sensitivity_analysis(self, specie, time):
        self.reintialize_simulation(setup_sensitivities=True)
        self.reactor_net.advance(time)
        print 'Done simulating reactor'
        s_list = map(lambda x: self.reactor_net.sensitivity(specie, x), range(len(self.solution.reactions())))
        sensitivities_mult = np.array(s_list)
        sensitivities = sensitivities_mult

        return sensitivities
    
    def get_full_sensitivity_data(self, tf=None, times='all'):
        sensitivity_list = []
        
        if tf is None:
            if times == 'all':
                raise ValueError('Must Specify a value for tf if all time points are desired')
                
            else:
                tf = times[-1]
        
        if times == 'all':
            time_list = [0]
            while self.reactor_net.time < tf:
                t_step = self.reactor_net.time - time_list[-1]
                time_list.append(self.reactor_net.time)

                if tf > self.reactor_net.time + 5*t_step:
                    self.reactor_net.step()
                else:
                    self.reactor_net.advance(tf)    
                sensitivity_list.append(self.reactor_net.sensitivities())

            time_list.append(self.reactor_net.time)
            time_list = time_list[2:]  # Remove padded zeros at the start of the list
            
        else:
            time_list = []
            for t in times:
                self.reactor_net.advance(t)
                sensitivity_list.append(self.reactor_net.sensitivities())
                time_list.append(self.reactor_net.time)
            
        n_time_points = len(sensitivity_list)
        sensitivity_tensor = np.zeros((self.solution.n_species+2, self.solution.n_reactions + self.solution.n_species, n_time_points))

        for t, sens_matrix in enumerate(sensitivity_list):
            sensitivity_tensor[:, :, t] = sens_matrix

        n_spcs = self.solution.n_species
        n_rxns = self.solution.n_reactions
        component_labels = [self.reactor_net.component_name(i) for i in range(n_spcs + 2)]
        parameter_labels = [self.reactor_net.sensitivity_parameter_name(i) for i in range(n_spcs + n_rxns)]
        sensitivity_data = SensitivityData(sensitivity_tensor, component_labels, parameter_labels, time_list)

        return sensitivity_data

    def get_highest_sensitive_to_params(self, sensitivity, n, specie_name):
        sensitivities_copy = abs(deepcopy(sensitivity))
        sens_indices = sensitivities_copy.argsort()[-n:][::-1]
        ordered_rxn_dict = OrderedDict()
        cantera_rxn_list = self.solution.reactions()

        for index in reversed(sens_indices):
            ordered_rxn_dict[cantera_rxn_list[index].equation] = sensitivity[index]

        fig = self.sensitivity_bar_chart(ordered_rxn_dict, specie_name)
        bp.show(fig)

        return ordered_rxn_dict

    @staticmethod
    def sensitivity_bar_chart(ordered_rxns_dictionary, sensitive_specie):
        rxn_list = ordered_rxns_dictionary.keys()
        sensitivity_values = ordered_rxns_dictionary.values()

        figure = bp.figure(y_range=rxn_list,
                           x_axis_label='Normalized Sensitivity of {0} to ki'.format(sensitive_specie))
        figure.hbar(y=np.arange(0.5, len(rxn_list)), height=0.5, left=0, right=sensitivity_values)
        return figure

    
class SensitivityData:
    
    def __init__(self, data_tensor, component_labels, parameter_labels, time_list):
        self.data = data_tensor
        self.components = component_labels
        self.parameters = parameter_labels
        self.times = time_list


def main():
    pass


if __name__ == '__main__':
    main()
