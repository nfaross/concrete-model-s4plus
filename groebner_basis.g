LoadPackage("GBNP");

N := 4;  # size of the magic unitary

A := FreeAssociativeAlgebraWithOne(Rationals, N * N, "R");
Generators := GeneratorsOfAlgebra(A);
A_One := Generators[1];  # unit in A

# construct a matrix M containing the generators
M := [];
for i in [1..N] do
    Row := [];
    for j in [1..N] do
        Add(Row, Generators[N * (i - 1) + j + 1]);
    od;
    Add(M, Row);
od;

# define the magic unitary relations
Relations := [];

for i in [1..N] do
    # sum of elements in the same column
    Rel := -A_One;
    for j in [1..N] do
        Rel := Rel + M[i][j];
    od;
    AddSet(Relations, GP2NP(Rel));
    # sum of elements in the same row
    Rel := -A_One;
    for j in [1..N] do
        Rel := Rel + M[j][i];
    od;
    AddSet(Relations, GP2NP(Rel));
od;

for i1 in [1..N] do for j1 in [1..N] do
    for i2 in [1..N] do for j2 in [1..N] do
        if i1 = i2 and j1 = j2 then
            # square of an element
            Rel := M[i1][j1] * M[i1][j1] - M[i1][j1];
            AddSet(Relations, GP2NP(Rel));
        elif i1 = i2 or j1 = j2 then
            # product in the same row / column
            Rel := M[i1][j1] * M[i2][j2];
            AddSet(Relations, GP2NP(Rel));
        fi;
    od; od;
od; od;

# compute and display the Groebner basis
GB := SGrobner(Relations);
Display(GB);
