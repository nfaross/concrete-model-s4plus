try:
    import queue
except ImportError:
    import Queue as queue
import resource
import sys

from algebra import I, mult, M

# description of the finite automaton generating the quotient basis of I_4
# it has the form {state: [(next_state, evaluated_label), ...], ...}
transitions = {
    0: [(1, M[0][0]), (2, M[0][1]), (3, M[0][2]), (4, M[1][0]), (5, M[1][1]),
        (6, M[1][2]), (7, M[2][0]), (8, M[2][1]), (9, M[2][2])],

    1: [(9, M[2][2]), (8, M[2][1]), (6, M[1][2]), (5, M[1][1])],
    2: [(4, M[1][0]), (6, M[1][2]), (7, M[2][0]), (9, M[2][2])],
    4: [(2, M[0][1]), (3, M[0][2]), (8, M[2][1]), (9, M[2][2])],
    5: [(1, M[0][0]), (3, M[0][2]), (7, M[2][0]), (9, M[2][2])],

    8: [(15, M[0][0]), (14, M[0][2])],
    7: [(13, M[0][1]), (14, M[0][2])],
    6: [(12, M[0][0]), (11, M[2][0])],
    3: [(10, M[1][0]), (11, M[2][0])],

    15: [(8, M[2][1]), (9, M[2][2])],
    12: [(6, M[1][2]), (9, M[2][2])],
    13: [(7, M[2][0]), (9, M[2][2])],
    10: [(3, M[0][2]), (9, M[2][2])],

    14: [(11, M[2][0])],
    11: [(14, M[0][2])],

    9: [(16, M[0][0])],
    16: [(9, M[2][2])],
}

def log(msg):
    """
    Prints the current matrix size and resource usage. Uses the global
    variables 'rows', 'columns' and 'non_zero_entries'.
    """
    res = resource.getrusage(resource.RUSAGE_SELF)
    form = "{} - rows: {} - columns: {} - entries: {} - memory: {} - time: {}"
    print(form.format(
        msg, len(rows), len(columns), non_zero_entries,
        res.ru_maxrss, res.ru_utime
    ))

if len(sys.argv) != 2:
    print("Usage: %s <maximum degree m>" % sys.argv[0])
    sys.exit(1)

max_degree = int(sys.argv[1]) # maximum degree 'm' to check
last_degree = 0               # used for logging

# sparse matrix \Psi_m
rows = {}            # maps a row index to non-zero column indices
columns = {}         # maps a column index to non-zero row indices
row_id = {}          # maps a tensor to a row index
non_zero_entries = 0 # used for logging

# Algorithm 1: breadth-first search
queue = queue.Queue()
queue.put((0, {(I, I, I): 1}, 0))

while not queue.empty():
    state, monomial, degree = queue.get()

    # log intermediate results when starting to process the next degree
    if degree > last_degree:
        log("degree %d" % last_degree)
    last_degree = degree

    column = len(columns)    # current column index
    columns[column] = set()

    for tensor, coeff in monomial.items():
        # add a new row if the tensor has not be seen before
        if tensor not in row_id:
            row_id[tensor] = len(rows)
            rows[len(rows)] = set()

        row = row_id[tensor]  # current row index

        # insert (row, column) into the sparse matrix
        columns[column].add(row)
        rows[row].add(column)
        non_zero_entries += 1

    if degree < max_degree:
        for next_state, label in transitions[state]:
            queue.put((next_state, mult(monomial, label), degree + 1))

log("degree %d" % max_degree)


# Algorithm 2: matrix elimination
stack = []
rank = 0    # lower bound on the rank of \Psi_m

for i in rows:
    if len(rows[i]) == 1:
        stack.append(i)

while stack:
    i = stack.pop()

    # row already completely eliminated
    if len(rows[i]) == 0:
        continue

    j, = rows[i]

    for k in columns[j]:
        if k != i:
            rows[k].remove(j)
            if len(rows[k]) == 1:
                stack.append(k)

    columns[j] = {i}
    rank += 1

log("finished")

# print result
print("The rank is at least %d of %d." % (rank, len(columns)))
if rank == len(columns):
    print("No such polynomial up to degree %d exists." % max_degree)
else:
    print("Such a polynomial might exist.")
