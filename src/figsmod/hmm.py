# Standard protein HMM
import numpy as np
from pomegranate import *
from time import time


def initialize_hmm(k):    
    number_of_states = k

    # Define the states
    states = []
    i_dist = DiscreteDistribution({'A': 0.05, 'C': 0.05, 'D': 0.05, 'E': 0.05,
                                   'F': 0.05, 'G': 0.05, 'H': 0.05, 'I': 0.05,
                                   'K': 0.05, 'L': 0.05, 'M': 0.05, 'N': 0.05,
                                   'P': 0.05, 'Q': 0.05, 'R': 0.05, 'S': 0.05,
                                   'T': 0.05, 'V': 0.05, 'W': 0.05, 'Y': 0.05})
    i_dist.freeze()
    for j in range(1, number_of_states + 1):
        m = State(DiscreteDistribution({'A': 0.05, 'C': 0.05, 'D': 0.05, 'E': 0.05,
                                        'F': 0.05, 'G': 0.05, 'H': 0.05, 'I': 0.05,
                                        'K': 0.05, 'L': 0.05, 'M': 0.05, 'N': 0.05,
                                        'P': 0.05, 'Q': 0.05, 'R': 0.05, 'S': 0.05,
                                        'T': 0.05, 'V': 0.05, 'W': 0.05, 'Y': 0.05}), name=f'm{j}')
        i = State(i_dist, name=f'i{j}')
        d = State(None, name=f'd{j}')
        states.extend([m, i, d])

    model = HiddenMarkovModel('protein HMM')
    model.add_states(states)

    # Define the transitions
    model.add_transition(model.start, states[0], 1)

    for j in range(number_of_states - 1):
        m = states[j * 3]
        i = states[j * 3 + 1]
        d = states[j * 3 + 2]
        m_next = states[(j + 1) * 3]
        i_next = states[(j + 1) * 3 + 1]
        d_next = states[(j + 1) * 3 + 2]
        model.add_transition(m, m_next, 0.9)
        model.add_transition(m, i, 0.05)
        model.add_transition(m, d_next, 0.05)
        model.add_transition(i, m_next, 0.8)
        model.add_transition(i, i, 0.2)
        model.add_transition(i, d_next, 0.2)
        model.add_transition(d, m_next, 0.8)
        model.add_transition(d, d_next, 0.2)

    model.add_transition(states[-3], model.end, 1)
    model.add_transition(states[-1], model.end, 1)
    model.bake()

    return model


