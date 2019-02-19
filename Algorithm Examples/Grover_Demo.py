## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#  Grover's Algorithm using ProjectQ
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from projectq import MainEngine
from projectq.ops import X, Z, H, Measure, All, ControlledGate
from projectq.meta import Loop
import math



# (x OR y) AND (~x OR y) AND (x OR ~y)
# xy = 00, 01, 10, 11

eng=MainEngine()
r=eng.allocate_qureg(2)
ancilla=eng.allocate_qureg(3)
q=eng.allocate_qubit()
num_it=int(math.pi/4*math.sqrt(2**len(r)))

# Step 1: Superposition
All(H) | r
X | q
H | q

# Step 2: Implement function as a circuit
with Loop(eng, num_it):
    X | r[0]
    X | r[1]
    X | ancilla[0]
    ControlledGate(X, 2) | (r[0], r[1], ancilla[0])
    X | r[0]
    X | r[1]


    X | r[1]
    X | ancilla[1]
    ControlledGate(X, 2) | (r[0], r[1], ancilla[1])
    X | r[1]


    X | r[0]
    X | ancilla[2]
    ControlledGate(X, 2) | (r[0], r[1], ancilla[2])
    X | r[0]

    # Step 3: Tag the correct state (00, 01, 10, 11)
    ControlledGate(X, 3) | (ancilla[0], ancilla[1], ancilla[2], q)

    # Uncompute (Step 2 in reverse)
    X | r[0]
    X | ancilla[2]
    ControlledGate(X, 2) | (r[0], r[1], ancilla[2])
    X | r[0]

    X | r[1]
    X | ancilla[1]
    ControlledGate(X, 2) | (r[0], r[1], ancilla[1])
    X | r[1]

    X | r[0]
    X | r[1]
    X | ancilla[0]
    ControlledGate(X, 2) | (r[0], r[1], ancilla[0])
    X | r[0]
    X | r[1]

    # Step 4: Amplitude amplification
    All(H) | r
    All(X) | r
    ControlledGate(Z, 1) | (r[0], r[1])
    All(X) | r
    All(H) | r

# Step 5: Measure
Measure | r
Measure | ancilla
Measure | q

eng.flush()
print([int(qubit) for qubit in r])






























    
    


























































































































