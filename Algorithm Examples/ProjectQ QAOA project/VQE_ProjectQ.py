## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Variational Quantum Eigensolver
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from projectq import MainEngine
from projectq.ops import QubitOperator, Measure, All, Rx, Ry

import numpy as np
from scipy.optimize import minimize

from collections import Counter

class VQE(object):
    '''
    The Variational Quantum Eigensolver algorithm

    VQE is an object that encapsulates the VQE algorithm (functional
    minimization). The main components of the VQE algorithm are a minimizer
    function for performing the functional minimization, a function that takes a
    set of parameters to initialize an ansatz state, and a Hamiltonian of which
    to calculate the expectation value.

    Using this object:

        1) initialize with `inst = VQE(minimizer)` where `minimizer` is a
        function that performs a gradient free minimuzation--i.e.
        scipy.optimize.minimize(., ., method='Nelder-Mead')

        2) call `inst.run(ansatz, hamiltonian, initial_parameters)`. Returns
        the optimal parameters and minimum expectation
    '''
    

    def __init__(self, minimizer = minimize, minimizer_args = [], minimizer_kwargs = {'method': 'Nelder-Mead'}):
        self.minimizer = minimizer
        self.minimizer_args = minimizer_args
        self.minimizer_kwargs = minimizer_kwargs

        '''
        Set up minimiser to be used in main VQE algorithm.
        Default is scipy's optimize.minimize. 
            
        Parameters:
            
            :minimizer: function that minimizes objective f(obj, param). For
            example the function scipy.optimize.minimize() needs at least two
            parameters, the objective and an initial point for the optimization.
            The args for minimizer are the cost function (provided by this class),
            initial parameters (passed to run() method). kwargs can be passed in below.
            
            :minimizer_args: (list) arguments to be passed to `minimizer` in
            form of *args (see `minimizer` above). Default = None
            
            :min_kwargs: (dict) arguments to be passed to `minimizer `in
            form of **kwargs (see `minimizer` above). Default = {'method': 'Nelder-Mead'}
            
            '''

        
        
    def run(self, ansatz, hamiltonian, initial_params, shots = 10000, eng = None):
        '''
        Main function. Applies the VQE algorithm with the given parameters and arguments
        using the minimizing function set up in self.__init__(). Finds the parameters
        that create the best ansatz state with lowest possible expectation value for the
        Hamiltonian.

        Parameters:

            :ansatz: (function: params, engine -> qureg) function that prepares a state to
            use to calculate expectation value of `hamiltonian`.

            :hamiltonian: (QubitOperator) the hamiltonian of which to calculate
            the expectation value.

            :initial_params: (list) initial parameters to provide to the minimizing function

            :shots: (int) number of times to sample in calculating expectation

            :eng: (MainEngine) the engine object we will be making calls to

        '''

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Check types
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        types = [QubitOperator, list, int]
        args = [hamiltonian, initial_params, shots]
        argstr = ['hamiltonian', 'initial_params', 'shots']

        for Type, arg in zip(types, args):
            if not isinstance(arg, Type):
                raise TypeError('Argument `{}` provided to expectation must be of type {}'.format(arg, Type))
            
        if shots <= 0:
            raise ValueError("`samples` variable must be a positive integer")

        
        if not callable(ansatz):
            raise TypeError('Argument `ansatz` provided to expectation must be a function of two arguments - the state preparation parameters and a QC engine.') 
    



        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Prepare Objective Function
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        def objective_function(params):
            exp = self.expectation(lambda eng: ansatz(params, eng), hamiltonian, shots, eng)
            
            self._current_expectation = exp # store for printing
            
            return exp
    

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Using minimzer to finally find best expectation value and params
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        final_result = self.minimizer(objective_function, initial_params, *self.minimizer_args, **self.minimizer_kwargs)

        return final_result


    

    def expectation(self, state_prep, hamiltonian, shots, eng):
        
        """
        Compute the expecation value of a given hamiltonian over the distribution
        generated by some given state_prep function.

        Parameters:
        
        
            :state_prep: The state preparation function used to calculate the
            expectation value.
        
            :hamiltonian: QubitOperator representing the hamiltonian of which
            to calculate the expectation value.
        
            :shots: The number of samples used to calculate the expectation value.
        
            :eng: The MainEngine object we will be making calls to.
        
            :return: A float representing the expectation value of hamiltonian given
            the state prepared by state_prep.
            
        
        """
        



        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Get expectation, term by term
        #
        # Hamiltonians are made up of Pauli gates: Is, Xs, Ys, and Zs
        #
        # Look through terms of Hamiltonian and add some correcting rotations to 
        # throw specific qubits back into the Z-basis (the native measuring basis)
        # when term is an X or a Y
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        expectation_total = 0
        for term, coefficient in hamiltonian.terms.items():

            rotations = []
            marked_qubits = []

            if term == ():
                expectation_total += coefficient

            for index, gate in term:
                marked_qubits.append(index)
                if gate == 'X':
                    rotations.append((Ry(-np.pi/2), index))
                elif gate == 'Y':
                    rotations.append((Rx(np.pi/2), index))



        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # run through circuit `shots` number of times
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            results = []
            for _ in range(shots):

            # Prepare state by instantiating state_prep
                q = state_prep(eng)
                
            # Apply correction rotations for 'X' and 'Y' gates in Hamiltonian
                for rotation, qubit_index in reversed(rotations):
                    rotation | q[qubit_index]

            # Measure and flush engine
                All(Measure) | q
                eng.flush()

            # Append each bitstring resulting from measurement to results
                results.append([int(qubit) for qubit in q])




        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Convert results into a dictionary with bitstring tuples as keys,
        # number of occurences as values
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            bitstring_tuples = list(map(tuple, results))
            freq = Counter(bitstring_tuples)

        

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Add or subtract from expectation value depending on number of 1's in bitstring
        # If number of 1s is odd, subtract from total expectation value. If even, add.
        #
        # Reasoning is that reading out a |0 > corresponds to an eigenvalue of +1, and
        # reading out a |1 > corresponds to an eigenvalue of -1. We're multiplying the
        # eigenvalues of each individual term, so if we had an odd number of terms we'd
        # have an overall -1 multiplying the entire state
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            expectation = 0
            for bitstring, count in freq.items():
                number_of_1s = 0
                for b, bit in enumerate(bitstring):
                    if bit == 1 and b in marked_qubits:
                        number_of_1s +=1
                if number_of_1s % 2 == 0:
                    expectation += coefficient*float(count)/shots
                else:
                    expectation -= coefficient*float(count)/shots

            expectation_total += expectation
        

        return expectation_total.real

        
    



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##             Test
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == "__main__":


    from projectq.ops import H, X
    import time
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # First define a state_prep function. This should take our initialized state from
    # |0 0 > ----> |+ 1 >
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def state_prep(eng):
        qureg = eng.allocate_qureg(2)
        H | qureg[0]
        X | qureg[1]
        return qureg

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Now let's create some example Hamiltonians.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ham_example1 = QubitOperator('X0 Z1', 3.0)
    ham_example2 = QubitOperator('X0 Z1', 3.0) + \
                           QubitOperator('Y1') + \
                           QubitOperator('Z1', 2)
    ham_example3 = QubitOperator('Z0')


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Let's just create an engine instance
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eng_example = MainEngine()


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # We should find expectation value to be ~ -3.0 if we use ham_example1,
    # and ~ -5.0 if we use ham_example2
    # Let's check it out:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\nTesting expectation method. Should get -5.0...')
    start = time.time()
    exp = VQE().expectation(state_prep, ham_example2, 1000, eng_example)
    end = time.time()
    print('\nExpectation (calculated in', end - start, 'seconds):\n', exp, '\n')




    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Let's test our run method
    #
    # First let's define some state_prep/ansatz functions 
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def ansatz_1(params, eng):
        qureg = eng.allocate_qureg(1)
        Ry(params[0]) | qureg[0]
        return qureg

    def ansatz_2(params, eng):
        qureg = eng.allocate_qureg(1)
        Rx(params[0]) | qureg[0]
        return qureg



    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Run a couple of examples, instantiate a VQE object
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vqe = VQE()
    print('\nTesting run method using Ry...')
    start = time.time()
    result = vqe.run(ansatz_1, ham_example3, [3], shots = 1000, eng = eng_example)
    end = time.time()
    print('Optimal eigenvalue:', result['fun'], '\nOptimal parameter:', result['x'], '\n(Calculated in', end - start, 'seconds)')

    print('\nTesting run method using Rx...')
    start = time.time()
    result = vqe.run(ansatz_2, ham_example3, [3], shots = 1000, eng = eng_example)
    end = time.time()
    print('Optimal eigenvalue:', result['fun'], '\nOptimal parameter:', result['x'], '\n(Calculated in', end - start, 'seconds)')
    
    

        

    
    






