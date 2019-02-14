'''
Uses the stored dataframe variables from leiva_2step_MM.py then creates updated variables based on a common tangent construction of the Gibbs free energy.

All 2phase variables have the 2p suffix appended.
'''

from __future__ import division
from sys import argv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp,log
from scipy.misc import factorial
from time import sleep
import argparse
import matplotlib as mpl
import string
from plotting_2l import Plotting
from leiva_2step_MM import Two_layer
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import numpy.polynomial.polynomial as poly

class Two_phase():
    def __init__(self,arg_dict):
       self.arg_dict = arg_dict
       self.tl_inst = Two_layer(arg_dict) 
       self.tl_inst.solution()
       self.df = self.tl_inst.df # Import the big dataframe from two layer.
       self.args = list(self.df)
       self.no_rows = self.df.shape[0] # Total row count.
       self.scaling = 10 # Factor to multiply by during interpolation step.
       self.args_2p = [(arg , arg.split('_')[0] + '_tp_' + arg.split('_')[-1]) for arg in self.args] # Tuple containing all the original args, together with the 
       self.val_list = []
       self.get_variables()
       self.convex_G = {} # Dictionary that is filled with list of index values where two phase coexistence occurs. Key is the interaction parameter.
       self.dict_of_intervals = {} # Contains all convex interval indices.
       self.dict_of_polys = {} # Contains all fitted polynomials. Return none if there 
       
    def get_variables(self):
       for entry in self.args:
           if len(entry.split('_')) > 1:
               value = entry.split('_')[-1]             
               if value not in self.val_list:
                   self.val_list.append(value)

    def get_convex_points(self):
        # Put in the stuff to allow common tangent to be calculated.
        for key in list(self.df):
            if key.startswith('dxdmu'):
                reduced_key = key.split('_')[-1] # Interaction parameter.
                self.convex_G[reduced_key] = np.array(self.df.index[self.df[key] >= 0]) # List of indices where G is convex, i.e. non-two phase regions.
            
    def get_secondderiv_nonzero(self):
        # self.list_of_tangents is a list of tuples, denoting how the common tangents are constructed
        self.get_convex_points()

        for key, indices in self.convex_G.iteritems():
            print 'indices=', indices
            self.list_of_intervals = [] # Will contain all of indices used for common tangent.
            differences = np.diff(indices)
#                print 'indices=', indices
                # Only perform operation where two phase coexistence is found.
            cross_over_points = np.argwhere([differences[i] > 1 for i in range(len(differences))])
            self.no_cross_overs = cross_over_points.size
            if len(cross_over_points) > 0:
                last_index_plus = 0
                for index in cross_over_points:
                    self.list_of_intervals.append(indices[last_index_plus : index[0]]) # Slice up.
                    last_index_plus = index[0] + 1
                self.list_of_intervals.append(indices[last_index_plus:]) # Gets final interval
            else:
                self.list_of_intervals.append('All') # In this case every point is convex.

            self.dict_of_intervals[key] = self.list_of_intervals

    def fit_polynomial(self, poly_order = 6):
        # Fits cubic polynomials to each of the intervals defined within list_of_intervals.
#                print self.list_of_intervals
        for key, intervals in self.dict_of_intervals.iteritems():
            print intervals
            if intervals == ['All']:
                self.dict_of_polys[key] = None
                self.df['spline_'+key] = np.linspace(self.df['G_' + key].iloc[0], self.df['G_' + key].iloc[-1], len(self.df.index))

                self.df['G_sub_' + key] = self.df['G_' + key] - self.df['spline_'+key] # Visulaisation of fit
                plt.plot(self.df['x'],self.df['G_sub_' +key],label ='G(x)')

                plt.legend()
            else:
#                self.fine_x = np.linspace(0,1,self.no_rows * self.scaling) # Generate a fitted function with more datapoints.
#                self.interpolated_G = interp1d(self.df['x'], self.df['G_' + k])(self.df['x']) # Generate finer mesh.
                self.df['spline_'+key] = np.linspace(self.df['G_' + key].iloc[0], self.df['G_' + key].iloc[-1], len(self.df.index))

                self.df['G_sub_' + key] = self.df['G_' + key] - self.df['spline_'+key] # Visulaisation of fit
                print self.df['x']
                self.dict_of_polys[key] = [] # Empty list to contain all fitted polynomials.
                for index,interval in enumerate(intervals):
                    xarr = self.df['x'].iloc[interval].values.astype(np.float64)
                    Garr = self.df['G_sub_' + key].iloc[interval].values.astype(np.float64)
                    pnomial = poly.polyfit(x = xarr,y = Garr, deg = poly_order)   # coefficients         
                    self.dict_of_polys[key].append(pnomial) # Empty list to contain all fitted polynomials.
                    self.df['f_'+str(index)] = poly.polyval(x = self.df['x'],c = pnomial,tensor=False)
                    plt.plot(self.df['x'].iloc[interval],self.df['f_'+str(index)].iloc[interval],label=index)
                plt.plot(self.df['x'],self.df['G_sub_' +key],label ='G(x)')

                plt.legend()
            
            plt.show()         
                
    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def straight_line(self,x):
        # Get a straight line tangent using a reference point on the curve as input.
        self.idx = self.find_nearest(self.fine_x,x)# Finish later!
        self.linear_function = np.array([self.interpolated_G[self.idx] + self.interpolated_dG[self.idx]*(x2 - x) for x2 in self.fine_x])
        return(self.linear_function)
        
    def min_intersection_1(x1,x2):
        # Determine the min rms distance between the fitting functions.
        diff = self.interpolated_G[0:self.min_indices[0]] - self.straight_line(x1)[0:self.min_indices[0]]
        return(diff)
    
#    def get_common_tangent(self):
#        self.get_secondderiv_zero()
#
#        print self.dict_of_intervals
# 
#        for k, index_pairs in self.dict_of_intervals.iteritems():
#            self.min_indices = [] # Indices indicating the minima.
#            self.df['G_tp_' + k] = self.df['G_' + k] # Copy of G
#            self.df['G_tp_' + k].iloc[replace_list] = tangent
#            self.df['G_raw_' + k] = self.df['G_' + k] - self.df['G_spline_' + k]

#            self.fine_x = np.linspace(0,1,self.no_rows * self.scaling) # Generate a fitted function with more datapoints.
#            self.interpolated_G = interp1d(self.df['x'], self.df['G_raw_' + k])(self.fine_x) # Generate finer mesh.
#            self.interpolated_dG = np.gradient(self.interpolated_G, self.fine_x, edge_order = 2) # Take gradient on the finer mesh.
            
#            for low_i, high_i in index_pairs:
 #               self.min_indices.append(self.df['dG_raw_' + k].iloc[low_i:high_i].min) # Find the minim#um of G (approximately).

 #           solution = least_squares(self.min_intersection_1(x1),self.min_indices[0]).x
 #           plt.plot(self.fine_x,solution)
 #           plt.show()
##           self.straight_line(x)
            
#            for point1, point2 in zip(index_pairs[0],index_pairs[1]):
#                input_point = np.array(self.min_indices[0,1])
#                rescaled_input = input_point * self.scaling # Equate point to the one on the interpolated function.
                
#                self.df['G_spline_' + k] = np.linspace(self.df['G_tp_'+k].iloc[0], self.df['G_tp_'+k].iloc[-1], len(self.df.index))

#                    self.df['G_raw_tp_' + k] = self.df['G_tp_'+k] - self.df['G_spline_' + k]                                 
#            plt.plot(self.df['x'],self.df['G_raw_' + k],label = k + ' raw')    
#            plt.plot(self.df['x'],self.df['G_raw_tp_'+k],label = k +' tangent')
#        else:
#            self.df[self.args_2p[1]] = self.df[self.args_2p[0]] # Put synonym names for two phase varaibles, where 1st order phase transition does not happen.
#        plt.legend()
#        plt.show()
                         
        
if __name__ == '__main__':
    # This is important if you wish to use the code well. Refer to the details on the "argparse" module for more information. On one line an entire plot, with one or several variables altered at time, can be generated. You need to change the argument "nargs" to alter the number of command line arguments.
    # There have been some strange numerical issues if the command line variables with multiple arguments are not entered from lowest to highest. Until I get to the bottom of this better to enter them from lowest to highest, i.e. "--g -0.5 -0.45 -0.4" for example.
    # You need to change "label" to the name of the varible that will be plotted, in this case "--label g".
    # Example: say you wanted a single plot, where g was varied as described above. You could enter the following on the command line.
    # "python leiva_2step_MM.py --M 300 --T 300 --nargs 3 --g -0.5 -0.45 -0.4 --label g"
    
    parser = argparse.ArgumentParser(description='Put in some energy parameters')
    parser.add_argument('--E0', type=float, help='Adds a point term, in kT', default = [-4.505], nargs='+')
    parser.add_argument('--g', type=float, help='In plane interaction parameter.', default = [-0.45], nargs='+')
    parser.add_argument('--delta1', type=float, help='Next layer interaction parameter', default = [1.12], nargs='+')
    parser.add_argument('--delE', type=float, help='Point term separation', default = [0.0], nargs='+')
    parser.add_argument('--loc', type=int, help='Legend_location', default = 0)
    parser.add_argument('--M',type=int, help = 'removable particles in each sublattice', default = 5)
    parser.add_argument('--T',type=float, help = 'temperature in Kelvin', default = [293.0], nargs='+')
    parser.add_argument('--nargs',type= int, help = 'number of arguments', default = 3)
    #parser.add_argument
    parser.add_argument('--label', type=str,help= 'variable to be labelled in plots.', default = 'g')
    parser.add_argument('--alpha4', type=float, help='amplitude of exponential decay, coupled to E0', default = [0.0], nargs = '+')
    parser.add_argument('--beta4', type=float, help='decay constant, coupled to E)', default = [0.0], nargs = '+')
    parser.add_argument('--alpha3', type=float, help='amplitude of exponential decay, coupled to g', default = [0.0], nargs = '+')
    parser.add_argument('--beta3', type=float, help='decay constant, coupled to g', default = [0.0], nargs = '+')
    parser.add_argument('--alpha1', type=float, help='amplitude of exponential decay, coupled to g', default = [0.0], nargs = '+')
    parser.add_argument('--beta1', type=float, help='decay constant, coupled to g', default = [0.0], nargs = '+')
#    parser.add_argument('--sigma', type =float, help='quadruplet parameter from Bazant.',default=[0.0], nargs = '+')
    
    args = parser.parse_args(argv[1:])
    
    E0 = args.E0 # Not important, just for ease of syntax.
    delE = args.delE
    g = args.g
    delta1 = args.delta1
    T = args.T
    M = args.M
    alpha4 = args.alpha4
    beta4 = args.beta4
    beta3 = args.beta3
    alpha3 = args.alpha3
    alpha1 = args.alpha1
    beta1 = args.beta1
#    sigma= args.sigma
    
    input_pars = ['E0','delE','g','delta1','alpha4','beta4','alpha3','beta3','alpha1','beta1','T','M']
    arg_dict = dict(vars(args))
    del arg_dict['loc']
    del arg_dict['nargs']

    for key,value in arg_dict.iteritems():
        if key in input_pars and key != 'M':
            if len(value) == args.nargs:
                var = key

    for key,value in arg_dict.iteritems():
        if key != 'label' and key != 'M':
            while len(value) < args.nargs:
                arg_dict[key].append(value[-1]) # This allows the variables that aren't changes to be specified once and then duplicated between plots/
    
    print 'Full list: ', arg_dict
            
    print 'E0 = ', E0
    
    counter = 0
    print list(arg_dict)
    two_phase_obj = Two_phase(arg_dict) # Instantiates claa.
    two_phase_obj.get_secondderiv_nonzero()
    two_phase_obj.fit_polynomial()
