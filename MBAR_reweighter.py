#!/usr/bin/python

from pickle import load, dump
from argparse import ArgumentParser
from sys import stdout
from pymbar import MBAR
from pymbar import timeseries
from numpy import *
import matplotlib.pyplot as plt
import pylab
import matplotlib.cm as cm

#######################################################################
#                     Notes and in-file settings                      #
#######################################################################

#For dihedrals CV with n16N PLUM, use -m 58*n_peptides

#in-file settings
pickle_datafiles  = True
pickle_free_ener  = True
pickle_f_i_tables = True
VERBOSE           = False
kB                = 0.001987 # Boltzmann constant in kcal/mol/K
temp_units        = 'K'

#######################################################################
#                              Functions                              #
#######################################################################

def logSum(log_terms):
   """Compute the log of a sum of terms whose logarithms are provided.

   REQUIRED ARGUMENTS  
      log_terms is the array (possibly multidimensional) containing the logs of the terms to be summed.

   RETURN VALUES
      log_sum is the log of the sum of the terms.

   """

   # compute the maximum argument
   max_log_term = log_terms.max()

   # compute the reduced terms
   terms = exp(log_terms - max_log_term)

   # compute the log sum
   log_sum = log( terms.sum() ) + max_log_term

   # return the log sum
   return log_sum

def histogram_wham(beta_k, U_kn, N_k, N_bins = 100, bin_width = None, maximum_iterations = 5000, relative_tolerance = 1.0e-8, initial_f_k = None):
   """Construct an initial guess of the f_k by histogram reweighting (specifically, WHAM [2]).

   ARGUMENTS
     beta_k (numpy K array) - inverse temperatures (in units of 1/energy)
     U_kn (numpy K x N_max array) - potential energies (in units of energy)
     N_k (numpy K array of numpy.int32) - number of samples per states, N_max = N_k.max()

   OPTIONAL ARGUMENTS
     N_bins (int) - if specified, sets the number of bins to use (default: 100)
     bin_width (float) - if specified, sets the bin width (overrides N_bins) (defulat: None)
     maximum_iterations (int) - maximum number of iterations to use
     relative_tolerance (floeat) - relative convergence tolerance (default: 1.0e-8)

   RETURNS
     f_k (numpy K array) - guess at initial state dimensionless free energies

   REFERENCE
     [2] Kumar S, Bouzida D, Swensen RH, Kollman PA, and Rosenberg JM. The weighted histogram analysis method
     for free-energy calculations on biomolecules. I. The Method. J. Comput Chem. 13:1011, 1992.
   
   """

   # Get sizes
   K = N_k.size
   N_max = N_k.max()

   # Create a list of indices of all configurations in kn-indexing.
   mask_kn = zeros([K,N_max], dtype=bool)
   for k in range(0,K):
      mask_kn[k,0:N_k[k]] = True
   # Create a list from this mask.
   sample_indices = where(mask_kn)

   # Construct histogram bins
   M = N_bins # number of energy bins
   SMALL = 1.0e-6
   U_min = U_kn[sample_indices].min()
   U_max = U_kn[sample_indices].max()
   U_max += (U_max - U_min) * SMALL # increment by a bit
   delta_U = (U_max - U_min) / float(M)
   if (bin_width != None):
      delta_U = bin_width
      M = int(ceil((U_max - U_min) / bin_width))
      #print "Using %d bins to achieve bin width of %f" % (M, delta_U)
   else:
      #print "Bin width is %f energy units to achieve N_bins = %d" % (delta_U, M)
      pass
   U_m = U_min + delta_U * (0.5 + arange(0,M,dtype = float64))
   H_m = zeros([M], float64)
   # assign snapshots to energy bins
   bin_kn = zeros([K,N_max], int32) # bin_kn[k,n] is bin index of U_kn[k,n]
   bin_kn[sample_indices] = array(((U_kn[sample_indices] - U_min) / delta_U), int32)
   H_m = bincount(bin_kn[sample_indices])
   
   # compute logs of various quantities
   LOG_ZERO = -1000.0 # replacement for log(0)
   log_H_m = ones([M], float64) * LOG_ZERO
   for m in range(M):
      if (H_m[m] > 0):
         log_H_m[m] = log(H_m[m])
   log_N_k = ones([K], float64) * LOG_ZERO
   for k in range(K):
      if (N_k[k] > 0):
         log_N_k[k] = log(N_k[k])

   # initialize free energies
   f_k = zeros([K], float64)
   if (initial_f_k != None):
      f_k = initial_f_k.copy()

   # iterate
   f_k_new = zeros([K], float64)
   max_delta = 0.0
   for iteration in range(maximum_iterations):
      # print "iteration %d" % iteration
      
      # Form auxiliary matrices, used in both self-consistent iteration and Newtom-Raphson.
      # log space
      log_w_mi = zeros([M,K], float64)
      for m in range(M):
         # denominator = \sum_k N_k exp[f_k - \beta_k U_m]
         exp_arg_k = log_N_k[:] + f_k[:] - beta_k[:]*U_m[m]
         log_w_mi[m,:] = exp_arg_k[:] - logSum(exp_arg_k[:])
      # real space
      w_mi = zeros([M,K], float64)
      for m in range(M):
         exp_arg_k = f_k[:] - beta_k[:]*U_m[m]
         max_arg = exp_arg_k.max()
         numerators_i = N_k[:] * exp(exp_arg_k[:] - max_arg)
         w_mi[m,:] = numerators_i[:] / numerators_i.sum()
         
      # Compute new estimates of log weights using self-consistent iteration.
      for i in range(K):
         # Compute new estimate of log weights.
         f_k_new[i] = f_k[i] + log(N_k[i]) - logSum(log_H_m[:] + log_w_mi[:,i])
      
      # shift to ensure f_k_new[0] = 0.0
      f_k_new -= f_k_new[0]
      # store difference
      Delta_f_k = f_k_new - f_k
      # update f_k
      f_k = f_k_new.copy()

      # If all f_k are zero, terminate
      if all(f_k == 0.0):
        break
      
      # Terminate when max((f - fold) / f) < relative_tolerance for all nonzero f.
      max_delta = (abs(Delta_f_k) / (abs(f_k_new)).max()).max()
      #print "iteration %8d relative max_delta = %8e" % (iteration, max_delta)
      if isnan(max_delta) or (max_delta < relative_tolerance):
         break

   # if unconverged
   if (iteration == maximum_iterations):
      print "Did not converge in %d iterations (final relative tolerance %e)" % (maximum_iterations, max_delta)
      print "f_k = "
      print f_k
      return f_k
   
   # summary
   #print "Converged to relative tolerance of %e (convergence tolerance %e) in %d iterations" % (max_delta, relative_tolerance, iteration)
   #print "f_k = "
   #print f_k

   # return the estimate of the dimensionless free energies
   return f_k

#######################################################################
#                         Parse input options                         #
#######################################################################

#Set up argument parser
parser = ArgumentParser()
parser.add_argument("T_file", help="A file with one line containing a space-separated list of temperatures of replicas, and a second line containing a list of temperatures at which to output results. If not present, the first line will be used for this again.")
parser.add_argument("PE_file", help="1 line per snapshot, with a space-separated list of the potential energy in each temperature.")
parser.add_argument("CV_file", nargs='+', help="1 file per CV containing 1 line per snapshot, with a space-separated list of the collective variable in each temperature.")
parser.add_argument('-i', nargs='?', help="Value to use as statistical inefficiency. Set to 1 for no subsampling.\
                                            Leave unset for this program to estimate i for you.")
parser.add_argument('-g', '--graph', help="Save graphs of the raw CV data histogram and of the reweighted CV probability distribution.", default=False, action='store_true')
parser.add_argument('-m', '--max-CV', nargs='?', help="Give the maximum value of the collective variable, if you want to normalise it.")
parser.add_argument('-b', '--n-bins', nargs='?', help="Number of bins per dimension in histogram and probability plots.", default=100)
parser.add_argument('-n', '--name-CV', nargs='*', help="Give the name of the collective variable which will appear in the graphs, in the output, and in the graphs' filenames. Spaces are supported if you surround with quotation marks.", default='Collective variables')

#Interpret input
args  = parser.parse_args()

N_bins = int(args.n_bins)
T_fn  = args.T_file
PE_fn = args.PE_file
CV_fns = args.CV_file
N_CVs = len(CV_fns)

try:
   stat_inefficiency = float(args.i)
except TypeError:
   stat_inefficiency = None

try:
   CV_hist_max = float(args.max_CV)
except TypeError:
   CV_hist_max = 1.0

CV_name = args.name_CV
if len(CV_name) != N_CVs:
   CV_name = [ 'Collective variable ' + str(d) for d in range(N_CVs) ]

graph = args.graph

args = None

#Count samples in PE_file
with open(PE_fn) as PE_file:
  for i, l in enumerate(PE_file):
      pass
N_samples = i + 1

#######################################################################
#          Produce arrays of temperature_k/beta_k of interest           #
#######################################################################

#Calculate temperature_k and beta_k
with open(T_fn, 'r') as T_file:
   temperature_in_strs = T_file.readline().split(' ')
   temperature_out_strs = T_file.readline().split(' ')
   try:
      assert ( temperature_out_strs[0] != '' )
   except AssertionError: #AssertionError will be thrown if second line is empty.
      temperature_out_strs = temperature_in_strs

N_replicas     = len(temperature_in_strs)
N_output_temps = len(temperature_out_strs)
temperature_k  = zeros([N_replicas+N_output_temps], float32)

#Set temperature_ks that data was taken at
for i,temp_str in enumerate(temperature_in_strs+temperature_out_strs):
   temperature_k[i] = float(temp_str)

beta_k = (kB * temperature_k)**(-1)

print ""

#######################################################################
#                Obtain {U,A}_kn_correlated from files                #
#######################################################################

data_pickle_fn = PE_fn+"+"+"_".join(CV_fns)+".pickle"
print "Reading in data:"
try: #Try to read the pickle; 1) assert that the user wants this and 2) try to open the file
   assert (pickle_datafiles == True)
   pickle_file = open(data_pickle_fn, 'rb')
   print "Loading PE data and", CV_name, "data from pickle file", data_pickle_fn
   (U_kn_correlated, A_ikn_correlated) = load( pickle_file )
   pickle_file.close()
except (IOError, AssertionError): #Obtain the PE and CV data from text files.
   print "(1 of 2) PE data"
   U_kn_correlated = zeros([N_replicas,N_samples], float32)
   with open(PE_fn, 'r') as PE_file:
      for i in range(N_samples):
         line = PE_file.readline()
         numbers = line.split()
         for k in range(N_replicas):
            U_kn_correlated[k][i] = float32( numbers[k] )

   print "(2 of 2)", CV_name
   A_ikn_correlated = zeros([N_CVs, N_replicas, N_samples], float32)
   for n, CV_fn in enumerate(CV_fns):
      with open(CV_fn, 'r') as CV_file:
         for i in range(N_samples):
            line = CV_file.readline()
            numbers = line.split()
            for k in range(N_replicas):
               A_ikn_correlated[n][k][i] = float32(numbers[k])/CV_hist_max

   if pickle_datafiles == True:
      print "Pickling resultant data with file name:", data_pickle_fn
      pickle_file = open(data_pickle_fn, 'wb')
      dump( (U_kn_correlated, A_ikn_correlated), pickle_file )
      pickle_file.close()

print ""

#######################################################################
#           Subsample {U,A}_kn_correlated to be uncorrelated          #
#######################################################################

print "Subsampling to achieve uncorrelated data"
if stat_inefficiency == None:
   print "(1 of 2) Calculating statistical inefficiency (i = ",
   stdout.flush()
   for d in range(N_CVs):
      statnew = timeseries.statisticalInefficiencyMultiple(A_ikn_correlated[d])
      stat_inefficiency = max([stat_inefficiency, statnew])
   print stat_inefficiency, ")"
else:
   print "(1 of 2) Using given statistical inefficiency (i =", str(stat_inefficiency) + ")"

indices = timeseries.subsampleCorrelatedData(U_kn_correlated[0,:], g = stat_inefficiency)
N_uncorrelated_samples = len(indices)

print "(2 of 2) Subsampling to achieve", N_uncorrelated_samples, "samples per replica"

U_kn  = zeros([      N_replicas+N_output_temps,N_uncorrelated_samples], float32)
A_ikn = zeros([N_CVs,N_replicas+N_output_temps,N_uncorrelated_samples], float32)

for k in range(N_replicas):
   U_kn[k] = U_kn_correlated[k][indices]
   for d in range(N_CVs):
      A_ikn[d][k] = A_ikn_correlated[d][k][indices]

print ""

#######################################################################
#                     Create histogram of CV data                     #
#######################################################################

if graph:
   if N_CVs == 1:

      every=5
      print "Creating", CV_name[0], "histogram graph"
      plt.figure(figsize=(14,5))
      plt.hist( [dataset for dataset in A_ikn[0][:N_replicas:every] ], N_bins, fill=False, histtype='step',
            label=[ str(temperature)+'K' for temperature in temperature_k[:N_replicas:every] ],
                 color= cm.Dark2(linspace(0, 1, N_replicas/every))                                  )
      plt.legend(loc=9)
      plt.grid(True)
      axis_bounds = plt.axis()
      plt.axis( [ axis_bounds[1], axis_bounds[0], axis_bounds[2], axis_bounds[3] ] )
      plt.title('Histogram of ' + str(CV_name[0]) + ' dataset')
      plt.xlabel(CV_name[0])
      plt.ylabel('Population of bin')
      graph_file_name = str(CV_name[0]) + ' Histogram'
      print "Saving as \'" + graph_file_name + "\'"
      plt.savefig(graph_file_name)
      print ""
      exit()

   #if N_CVs == 2:
      #print "Creating", CV_name, "2D histogram graph of lowest temperature."
      #heatmap, xedges, yedges = histogram2d(A_ikn[0][0], A_ikn[1][0], bins=N_bins)
      #fig = plt.figure()
      #ax = fig.add_subplot(1,1,1)
      #ax.set_ylabel(CV_name[0])
      #ax.set_xlabel(CV_name[1])
      #imageplot = ax.imshow(heatmap)
      #imageplot.set_interpolation('nearest')

      #plt.show()

#######################################################################
#                           Calculate u_kln                           #
#######################################################################

print "Calculating u_kln"

# u_kln[k,l,n] is reduced potential energy of sample n in temperature replica k evaluated at temperature l
u_kln = zeros([N_replicas+N_output_temps,N_replicas+N_output_temps,N_uncorrelated_samples], float32)
N_k = zeros([N_replicas+N_output_temps], int32)

for k in range(N_replicas):
   N_k[k] = N_uncorrelated_samples
   for l in range(N_replicas+N_output_temps):
      u_kln[k,l, 0:N_k[k]] = beta_k[l] * U_kn[k, 0:N_k[k]]

print ""

#######################################################################
#              Instantiate MBAR and compute expectations              #
#######################################################################

print "Obtaining free energies from MBAR"
free_E_pickle_fn = "free_E@"+str(stat_inefficiency)+".pickle"
try: #Try to read the pickle; 1) assert that the user wants this and 2) try to open the file
   assert (pickle_free_ener == True)
   pickle_file = open(free_E_pickle_fn, 'rb')
   print "(1 of 2) Loading free energy estimates from pickle file:", free_E_pickle_fn
   f_k = load( pickle_file )
   mbar = MBAR(u_kln, N_k, use_optimized = True, verbose=VERBOSE, initial_f_k = f_k, maximum_iterations = 0)
   pickle_file.close()
   print "MBAR f_k"
   print mbar.f_k

except (AssertionError, IOError):
   print "(1 of 2) Calculating free energy differences"
   f_k = histogram_wham ( beta_k[:N_replicas], U_kn[:N_replicas], N_k[:N_replicas] )
   print "WHAM f_k"
   print f_k
   f_k_estimate = f_k.resize(N_replicas + N_output_temps)
   mbar = MBAR(u_kln, N_k, use_optimized = True, verbose=VERBOSE, initial_f_k = f_k)
   print "MBAR f_k"
   print mbar.f_k

   if pickle_free_ener == True:
      print "Pickling free energy estimates with file name:", free_E_pickle_fn
      pickle_file = open(free_E_pickle_fn, 'wb')
      dump( mbar.f_k , pickle_file )
      pickle_file.close()

print "(2 of 2) Calculating ensemble average values and uncertainties."
#A_k_estimated = ndarray(N_CVs, dtype=ndarray)
#for i, A_kn in enumerate(A_ikn):
   #A_k_estimated[i] = mbar.computeExpectations(A_kn)[0]
print ""

#######################################################################
#                       Compute N-D PMF for CV                        #
#######################################################################

print "Computing PMF for CVs"
print "(1 of 2) Calculating histograms"

u_kn = ones( (N_replicas, N_uncorrelated_samples), int32 )

CV_hist_min = ndarray(N_CVs, dtype=ndarray)
CV_hist_max = ndarray(N_CVs, dtype=ndarray)
bin_kn      = zeros( u_kn.shape , int32 )

for i, A_kn in enumerate(A_ikn):
   CV_hist_min[i] = min(A_kn[:N_replicas].flatten())
   CV_hist_max[i] = max(A_kn[:N_replicas].flatten())+0.001

#Iterate through bins and assign each bin to all relevant A_kns at once.
#If a bin is empty, expand it until it's not.
bin_index = 0
bin_count = []
bin_number = zeros( [N_bins]*N_CVs, dtype=int )

bin_bounds = array([ linspace ( CV_hist_min[d], CV_hist_max[d], num=N_bins+1 ) for d in range(N_CVs) ])

coords_low = [ 0 for d in range(N_CVs) ]

#Binning loop
for i in range(N_bins**N_CVs):
   coords_high = [ (i/N_bins**(N_CVs-d-1))-N_bins*(i/(N_bins**(N_CVs-d)))+1 for d in range (N_CVs) ]
   coords_low = [ [ (i/N_bins**(N_CVs-d-1))-N_bins*(i/(N_bins**(N_CVs-d))) ] for d in range (N_CVs) ]
   bin_high    = [ bin_bounds[d][coords_high[d]] for d in range (N_CVs) ]
   bin_low     = [ bin_bounds[d][coords_low[d] ] for d in range (N_CVs) ]

   in_this_bin = ones([N_replicas, N_uncorrelated_samples], dtype=bool)
   for d in range(N_CVs):
      in_this_bin &= ( (A_ikn[d][:N_replicas] >= bin_low[d]) & (A_ikn[d][:N_replicas] < bin_high[d]) )

   indices = where(in_this_bin == True)
   counts = len(indices[0])

   if counts > 0:
      bin_number[coords_low] = bin_index

      bin_count.append( counts )
      bin_kn[indices] = bin_index
      bin_index += 1
   else:
      bin_number[coords_low] = -1

f_i_pickle_fn = "f_i-of-"+"_".join(CV_fns)+str(stat_inefficiency)+"@"+str(N_bins)+".pickle"
try: #Try to read the pickle; 1) assert that the user wants this and 2) try to open the file
   assert (pickle_f_i_tables == True)
   pickle_file = open(f_i_pickle_fn, 'rb')
   print "(2 of 2) Loading f_i estimates from pickle file:", f_i_pickle_fn
   (f_i_tables, f_i_uncerts_tables) = load( pickle_file )
   pickle_file.close()

except (AssertionError, IOError):
   print "(2 of 2) Reweighting"
   f_is = []
   f_i_uncerts = []

   for i, target_beta in enumerate(beta_k[N_replicas:]):
      u_kn = target_beta * U_kn
      (f_i, d2f_ij) = mbar.computePMF(u_kn, bin_kn, bin_index)
      f_i = append(f_i, inf)
      d2f_ij = append(d2f_ij, 0.0)

      f_is.append(f_i)
      f_i_uncerts.append(d2f_ij)

   f_i_tables         = zeros([N_output_temps]+[N_bins]*N_CVs, dtype = float)
   f_i_uncerts_tables = zeros([N_output_temps]+[N_bins]*N_CVs, dtype = float)

   for i in range(N_output_temps):
      f_i_tables[i][:] = f_is[i][ bin_number[:] ]
      f_i_uncerts_tables[i][:] = f_i_uncerts[i][ bin_number[:] ]

   if pickle_f_i_tables == True:
      print "Pickling free energy estimates with file name:", f_i_pickle_fn
      pickle_file = open(f_i_pickle_fn, 'wb')
      dump( (f_i_tables, f_i_uncerts_tables), pickle_file )
      pickle_file.close()

p_i_tables = exp(-f_i_tables)
for i in range(N_output_temps): #normalise p_i
   p_i_tables[i] = p_i_tables[i]/p_i_tables[i].sum()

p_i_uncerts_tables = f_i_uncerts_tables*p_i_tables

if graph:
   if N_CVs == 1:
      print ""
      print "Creating", CV_name[0], "probability distribution graph"
      plt.clf()
      #print p_i_uncerts_tables
      for i, p_i in enumerate(p_i_tables):
         plt.errorbar(bin_bounds[0][0:-1], p_i, label=str(temperature_k[N_replicas+i])+'K',
                 linestyle='steps',
                 #yerr=p_i_uncerts[i],
                 elinewidth=1.0,
                 color= cm.Dark2(linspace(0, 1, N_output_temps))[i],
                 linewidth=1.0,
                 barsabove=True)

      #Graph settings
      plt.legend(loc=1)
      axis_bounds = plt.axis()
      plt.axis( [ axis_bounds[1], axis_bounds[0], axis_bounds[2], axis_bounds[3] ] )
      plt.grid(True)
      plt.title('Probability distribution of ' + str(CV_name[0]))
      ylims = plt.ylim()
      plt.ylim(0.0,ylims[1])
      plt.xlabel(CV_name[0])
      plt.ylabel('Probability of state')

      graph_file_name = str(CV_name[0]).replace(" ","_") + '_p_i' + str(stat_inefficiency) + ".png"
      print "Saving as \'" + graph_file_name + "\'"
      plt.savefig(graph_file_name)

   if N_CVs == 2:
      print ""
      #plt.clf()

      X = zeros( [N_bins+1]*N_CVs, dtype=float )
      Y = zeros( [N_bins+1]*N_CVs, dtype=float )

      for i in range(N_bins+1):
         for j in range(N_bins+1):
            X[i][j] = bin_bounds[0][j]
            Y[i][j] = bin_bounds[1][i]

      for i in range(N_output_temps):
         print "Creating", CV_name, "2D probability distribution graph for", str(temperature_k[N_replicas+i]) + temp_units
         prob_max = abs(p_i_tables[i]).max()
         free_max = abs(f_i_tables[i][where(f_i_tables[i] != inf)]).max()
         free_max = 5

         fig = plt.figure()

         ax = fig.add_subplot(111)
         plt.axis([X.min(), X.max(), Y.min(), Y.max()])

         mygnuplot = cm.gnuplot
         mygnuplot.set_over(color=(1,1,1))

         graph = ax.pcolor(X,Y,f_i_tables[i], vmin=0, vmax=free_max, cmap=mygnuplot)
         fig.colorbar(graph)

         ax.grid(True)
         ax.set_ylabel(CV_name[0])
         ax.set_xlabel(CV_name[1])
         ax.set_title('Free energy Ramachandran plot of Gly-Ala-Gly peptide')

         plt.savefig('free_e_' + str(stat_inefficiency) + ".png")

print ""

#######################################################################
#                  Write the requested data to file                   #
#######################################################################

#averages_fn = CV_name.replace(" ","_") + '_ensemble_average'
#probability_distribution_fn = CV_name.replace(" ","_") + '_probability_distribution'

#with open( averages_fn, 'w' ) as outfile:
   #print "Saving datafile of ensemble average values as", averages_fn
   #for i in range(N_replicas, N_replicas+N_output_temps):
      ##outfile.write ('{0:<{len}} '.format(temperature_k[i], len=6))
      #outfile.write ( '{0:{len}} {1:{len}}'.format(A_k_estimated[i], dA_k_estimated[i], len=10) )
      #outfile.write ("\n")

##with open( probability_distribution_fn, 'w' ) as outfile:
   ##print "Saving datafile of probability distribution as", probability_distribution_fn

   ##for i in range(N_replicas, N_replicas+N_output_temps):
      ##outfile.write ('{0:<{len}} '.format(temperature_k[i], len=6))
      ##outfile.write ( '{0:{len}} {1:{len}}'.format(A_k_estimated[i], dA_k_estimated[i], len=10) )
      ##outfile.write ("\n")
   ##outfile.close()
