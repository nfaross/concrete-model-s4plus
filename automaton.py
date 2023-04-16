import ast
import sage.all as sage

# read the output of 'quotient.g' from the file 'groebner_basis.txt'
with open("basis4.txt") as file:
    data = file.read()

basis = ast.literal_eval(data)

n = 4  # size of the magic unitary
variables = range(1, n * n + 1)
leading_terms = [tuple(monomials[0]) for monomials, _ in basis]

# construction of the non-deterministic finite automaton
initial_state = tuple()
transitions = []

# add self-loops to the initial state
for variable in variables:
    transitions.append((initial_state, initial_state, variable))

for word in leading_terms:
    # add a path corresponding to the leading term
    for i in range(len(word)):
        transitions.append((word[:i], word[:i+1], word[:i+1][-1]))

    # add self-loops to the final state
    for variable in variables:
        transitions.append((word, word, variable))

automaton = sage.Automaton(
    transitions,
    initial_states=[initial_state],
    final_states=leading_terms,
    input_alphabet=variables
)

# transform the automaton
automaton = automaton.minimization().complement()

# remove the 'dead state'
for state in automaton.states():
    is_dead_end = all(trans.to_state == state for trans in state.transitions)
    if not state.is_final and is_dead_end:
        break

automaton.delete_state(state)

# print the final automaton, adjacency matrix and number of paths
automaton = automaton.relabeled()
M = automaton.adjacency_matrix(entry=lambda x: 1)

print(automaton)
print(M)
print(automaton.number_of_words())
