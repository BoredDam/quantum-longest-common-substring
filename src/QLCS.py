#!pip install qiskit qiskit-ibm-runtime qiskit-aer numpy
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit.library import PermutationGate
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import SamplerV2 as Sampler
import numpy as np


def match_operator(n: int, char_size: int) -> QuantumCircuit:
    qx = QuantumRegister(n * char_size, "x")
    qy = QuantumRegister(n * char_size, "y")
    qout = QuantumRegister(n, "o")
    qc = QuantumCircuit(qx, qy, qout)

    for i in range(n * char_size):
        qc.cx(qx[i], qy[i])
        qc.x(qy[i])

    for i in range(0, n):
        ctrls = list(range(i * char_size, i * char_size + char_size))
        qc.mcx([*qy[ctrls]], qout[i])
    return qc


def extension_operator(n: int, i: int) -> QuantumCircuit:
    """
    builds the quantum circuit for the extension operator for an array of size `n`
    for `i`-th lambda/L vector
    """

    shift = 2 ** (i - 1)

    qin = QuantumRegister(n, "in")
    qout = QuantumRegister(n, "out")
    qc = QuantumCircuit(qin, qout)

    for j in range(n - shift):
        qc.ccx(qin[j], qin[j + shift], qout[j])

    return qc


def quantum_rotation(n: int, s: int):
    """
    builds the rotation operator on `n` qubits by `s` positions
    """
    s %= n
    qc = QuantumCircuit(n)

    perm = [((k - s) % n) for k in range(n)]

    qc.append(PermutationGate(perm), range(n))
    return qc.decompose().to_gate(label=f"R_{s}")


def quantum_parametric_rotation(n: int) -> QuantumCircuit:
    """
    builds the quantum circuit for the parametric rotation gate on `n` qubits.
    The amount of control qubits it's equal to ceil(log2(n))
    """

    l = int(np.ceil(np.log2(n)))

    jr = QuantumRegister(l, "j")
    qr = QuantumRegister(n, "q")
    qc = QuantumCircuit(jr, qr)

    for i in range(l):
        qc = qc.compose(quantum_rotation(n, 2**i).control(1), [jr[i], *qr])
    return qc


def quantum_character_rotation(n: int, char_size: int) -> QuantumCircuit:
    """
    builds the quantum circuit for the parametric rotation on characters of
    size `char_size` qubits of `n` positions
    """

    l = int(np.ceil(np.log2(n)))

    jr = QuantumRegister(l, "j")
    qr = QuantumRegister(n * char_size, "q")
    qc = QuantumCircuit(jr, qr)

    for i in range(l):
        controlled_rot = quantum_rotation(n, 2**i).control(1)

        for bit in range(char_size):
            lines = [qr[k * char_size + bit] for k in range(n)]
            qc = qc.compose(controlled_rot, [jr[i], *lines])

    return qc


def copy_operator(n: int):
    """
    builds a copy gate on `n` qubits
    """

    qin = QuantumRegister(n, "in")
    qout = QuantumRegister(n, "out")
    qc = QuantumCircuit(qin, qout)

    for i in range(n):
        qc.cx(qin[i], qout[i])

    return qc.to_gate(label="COPY")


def or_operator(n: int):
    """
    builds a logical or gate on `n` qubits
    """

    qx = QuantumRegister(n, "x")
    qy = QuantumRegister(1, "y")
    qc = QuantumCircuit(qx, qy)

    qc.x(qx)
    qc.mcx(qx, qy)
    qc.x(qx)
    qc.x(qy)

    return qc.to_gate(label="OR")


def contr_bitwise_and_operator(n: int):
    """
    builds a controlled and operator on `n` qubits
    """

    qc_ctrl = QuantumRegister(1, "c")
    qL = QuantumRegister(n, "L")
    qin = QuantumRegister(n + 1, "qin")
    qout = QuantumRegister(n + 1, "qout")
    qc = QuantumCircuit(qc_ctrl, qL, qin, qout)

    for k in range(n):
        qc.mcx([qc_ctrl[0], qL[k], qin[k]], qout[k])

    return qc.to_gate(label="C-AND")


def FSM(n, char_size) -> QuantumCircuit:
    """
    builds the Fixed Substring Matching circuit for inputs of size `n` and
    characters of `char_size` qubits
    """

    d_len = int(np.ceil(np.log2(n)))
    L_num = d_len
    qx = QuantumRegister(n * char_size, "x")
    qy = QuantumRegister(n * char_size, "y")
    qd = QuantumRegister(d_len, "d")
    qD = QuantumRegister(d_len + 1, "D")
    qL = []

    for i in range(L_num):
        qL.append(QuantumRegister(n, "L[" + str(i) + "]"))

    qD = []
    for i in range(d_len + 1):
        qD.append(QuantumRegister(n + 1, "D[" + str(i - 1) + "]"))

    qr = QuantumRegister(1, "r")
    qc = QuantumCircuit(qx, qy, qd, *qL, *qD, qr)

    qc = qc.compose(
        match_operator(n, char_size).to_gate(label="MATCH"), [*qx, *qy, *qL[0]]
    )

    for i in range(L_num - 1):
        qc = qc.compose(
            extension_operator(n, i + 1).to_gate(label=f"EXT_{i}"),
            [*qL[i][:n], *qL[i + 1]],
        )

    for i in range(d_len):

        qc = qc.compose(
            contr_bitwise_and_operator(n), [qd[i], *qL[i], *qD[i], *qD[i + 1]]
        )

        controlled_rot = quantum_rotation(n + 1, 2**i).control(1)
        qc = qc.compose(controlled_rot, [qd[i], *qD[i + 1]])
        qc.x(qd[i])
        controlled_copy = copy_operator(n + 1).control(1)
        qc = qc.compose(controlled_copy, [qd[i], *qD[i], *qD[i + 1]])
        qc.x(qd[i])

    qc = qc.compose(or_operator(n + 1), [*qD[d_len], qr])

    return qc


def SFC(n, char_size) -> QuantumCircuit:
    """
    builds the Shared Fix Substring Check circuit for inputs of size `n` and
    characters of `char_size` qubits
    """

    d_len = int(np.ceil(np.log2(n)))
    L_num = d_len

    qx = QuantumRegister(n * char_size, "x")
    qy = QuantumRegister(n * char_size, "y")
    qd = QuantumRegister(d_len, "d")
    qD = QuantumRegister(d_len + 1, "D")
    qL = []

    for i in range(L_num):
        qL.append(QuantumRegister(n, "L[" + str(i) + "]"))

    qD = []
    for i in range(d_len + 1):
        qD.append(QuantumRegister(n + 1, "D[" + str(i - 1) + "]"))

    qr = QuantumRegister(1, "r")
    qc = QuantumCircuit(qx, qy, qd, *qL, *qD, qr)

    qc.x(qD[0])

    qc = qc.compose(FSM(n, char_size))
    return qc


def FPM(n, char_size) -> QuantumCircuit:
    """
    builds the Fixed Prefix Matching operator for inputs of size `n`
    """

    d_len = int(np.ceil(np.log2(n)))
    L_num = d_len
    qx = QuantumRegister(n * char_size, "x")
    qy = QuantumRegister(n * char_size, "y")
    qd = QuantumRegister(d_len, "d")
    qD = QuantumRegister(d_len + 1, "D")
    qL = []

    for i in range(L_num):
        qL.append(QuantumRegister(n, "L[" + str(i) + "]"))

    qD = []
    for i in range(d_len + 1):
        qD.append(QuantumRegister(n + 1, "D[" + str(i - 1) + "]"))

    qr = QuantumRegister(1, "r")
    qc = QuantumCircuit(qx, qy, qd, *qL, *qD, qr)

    # init
    qc.x(qD[0][0])

    qc = qc.compose(FSM(n, char_size))
    return qc


def grover_operator(n: int):
    """
    builds the diffusion operator from the grover algorithm
    """

    qx = QuantumRegister(n, "x")
    qc = QuantumCircuit(qx)

    qc.h(qx)
    qc.x(qx)

    if n == 1:
        qc.z(qx[0])
    else:
        qc.h(qx[n - 1])
        qc.mcx(qx[list(range(n - 1))], qx[n - 1])
        qc.h(qx[n - 1])

    qc.x(qx)
    qc.h(qx)

    return qc.to_gate(label="DIFFUSION")


def int_to_bin(num: int, code_len: int) -> str:
    """
    converts `num` in a binary string, with the specified `code_len`
    """

    num = bin(num)[2:]
    num = "0" * (code_len - len(num)) + num
    return num


def quantum_number_encode(num: int, code_len: int = 0):
    """
    encodes the given `num` into a gate. if specified, `code_len` lets
    you use more qubits than needed for the number specified by `num`
    """

    if num == 0:
        return QuantumCircuit(code_len).to_gate(label=" 0 ")

    if code_len < int(np.ceil(np.log2(num))):
        code_len = len(bin(num)[2:])

    bin_num = int_to_bin(num, code_len)
    qc = QuantumCircuit(code_len)

    for c, i in enumerate(bin_num):
        if i == "1":
            qc.x(code_len - 1 - c)

    return qc.to_gate(label=" " + str(num) + " ")


def encode_boolean_string(x: str) -> QuantumCircuit:
    """
    encodes the given boolean string `x` in a quantum circuit
    """
    n = len(x)
    qx = QuantumRegister(n, "x")
    qc = QuantumCircuit(qx)

    for c, i in enumerate(x):
        if i == "1":
            qc.x(qx[c])

    return qc


def quantum_step_circuit(
    x: str, y: str, n: int, d: int, char_size: int, grover_epochs: int
) -> QuantumCircuit:
    """
    builds the quantum circuit for the QLCS_test on the strings `x` and `y` of length `n`,
    checking for the common substring of length `d`, with the specified `char_size`
    """

    qi_amount = int(np.ceil(np.log2(n)))

    qi = QuantumRegister(qi_amount, "i")
    qj = QuantumRegister(qi_amount, "j")
    qx = QuantumRegister(n * char_size, "x")
    qy = QuantumRegister(n * char_size, "y")

    ancilla_amount = n * int(np.ceil(np.log2(n))) + (n + 1) * (
        int(np.ceil(np.log2(n + 1)))
    )

    qa = QuantumRegister(ancilla_amount, "a")
    qd = QuantumRegister(qi_amount, "d")

    qexp = QuantumRegister(2, "dom_exp")
    qr = QuantumRegister(1, "r")
    qo = QuantumRegister(1, "out")
    c = ClassicalRegister(1, "c")

    qc = QuantumCircuit(qexp, qi, qj, qx, qy, qa, qd, qr, qo, c)

    qc = qc.compose(quantum_number_encode(d, qi_amount), qd)
    qc = qc.compose(encode_boolean_string(x), qx)
    qc = qc.compose(encode_boolean_string(y), qy)

    qc.h(qexp[0])
    qc.barrier(label="init")

    # search phase
    qc.h(qj)
    qc.x(qo)
    qc.h(qo)

    for _ in range(grover_epochs):
        qc = qc.compose(
            quantum_character_rotation(n, char_size).to_gate(label="ROT"), [*qj, *qx]
        )
        qc = qc.compose(
            SFC(n, char_size).to_gate(label="SFC"), [*qx, *qy, *qd, *qa, qr]
        )

        qc.ccx(qexp[0], qr, qo)
        qc = qc.compose(
            SFC(n, char_size).inverse().to_gate(label="- SFC"), [*qx, *qy, *qd, *qa, qr]
        )
        qc = qc.compose(
            quantum_character_rotation(n, char_size).inverse().to_gate(label="- ROT"),
            [*qj, *qx],
        )

        qc = qc.compose(grover_operator(qi_amount + 1), [qexp[0], *qj])

    for j in qj:
        qc.measure(j, c)

    qc.barrier(label="SEARCH PHASE")
    qc.h(qexp[1])
    qc.h(qi)

    # verification phase
    for _ in range(grover_epochs):
        qc = qc.compose(
            quantum_character_rotation(n, char_size).to_gate(label="ROT"), [*qj, *qx]
        )
        qc = qc.compose(
            quantum_character_rotation(n, char_size).to_gate(label="ROT"), [*qi, *qx]
        )
        qc = qc.compose(
            quantum_character_rotation(n, char_size).to_gate(label="ROT"), [*qi, *qy]
        )
        qc = qc.compose(
            FPM(n, char_size).to_gate(label="FPM"), [*qx, *qy, *qd, *qa, qr]
        )

        qc.ccx(qexp[1], qr, qo)

        qc = qc.compose(
            FPM(n, char_size).inverse().to_gate(label="- FPM"), [*qx, *qy, *qd, *qa, qr]
        )
        qc = qc.compose(
            quantum_character_rotation(n, char_size).inverse().to_gate(label="- ROT"),
            [*qi, *qy],
        )
        qc = qc.compose(
            quantum_character_rotation(n, char_size).inverse().to_gate(label="- ROT"),
            [*qi, *qx],
        )
        qc = qc.compose(
            quantum_character_rotation(n, char_size).inverse().to_gate(label="- ROT"),
            [*qj, *qx],
        )

        qc = qc.compose(grover_operator(qi_amount + 1), [qexp[1], *qi])

    for i in qi:
        qc.measure(i, c)
    qc.barrier(label="VERIFICATION PHASE")

    # final check
    qc.h(qo)
    qc.x(qo)

    qc = qc.compose(
        quantum_character_rotation(n, char_size).to_gate(label="ROT"), [*qj, *qx]
    )
    qc = qc.compose(
        quantum_character_rotation(n, char_size).to_gate(label="ROT"), [*qi, *qx]
    )
    qc = qc.compose(
        quantum_character_rotation(n, char_size).to_gate(label="ROT"), [*qi, *qy]
    )
    qc = qc.compose(FPM(n, char_size).to_gate(label="FPM"), [*qx, *qy, *qd, *qa, qr])
    qc.cx(qr, qo)

    qc = qc.compose(
        FPM(n, char_size).inverse().to_gate(label="- FPM"), [*qx, *qy, *qd, *qa, qr]
    )
    qc = qc.compose(
        quantum_character_rotation(n, char_size).inverse().to_gate(label="- ROT"),
        [*qi, *qy],
    )
    qc = qc.compose(
        quantum_character_rotation(n, char_size).inverse().to_gate(label="- ROT"),
        [*qi, *qx],
    )
    qc = qc.compose(
        quantum_character_rotation(n, char_size).inverse().to_gate(label="- ROT"),
        [*qj, *qx],
    )

    qc.measure(qo, c)
    qc.barrier(label="FINAL CHECK")

    return qc


def run(qc: QuantumCircuit, shots=1024):
    """
    returns the result of the execution of the `qc` quantum circuit
    """
    sim_backend = AerSimulator(
        method="matrix_product_state"
    )  # the only viable method for circuits of this size with no entanglement
    sim_backend.set_max_qubits(qc.num_qubits)  # as big as needed

    pm = generate_preset_pass_manager(backend=sim_backend, optimization_level=2)
    isa_qc = pm.run(qc)
    sampler = Sampler(mode=sim_backend)

    job = sampler.run([isa_qc], shots=shots)
    results = job.result()
    data = results[0].data
    counts = data.c.get_counts()

    return counts


RUNS_PER_CIRCUIT = 2


def grover_iter_search(x: str, y: str, n: int, d: int, char_size: int):
    """
    implements a pseudo-BBHT approach to the grover search
    """
    m = 1
    lam = 1.34
    while True:
        j = int(m)
        qc = quantum_step_circuit(x, y, n, d, char_size, grover_epochs=j)
        c = run(qc, RUNS_PER_CIRCUIT)
        if "1" in list(c.keys()):
            return True

        m = min(lam * m, np.sqrt(n))
        if m >= np.sqrt(n):
            break

    return False


def quantum_LCS_test(x: str, y: str, n: int, d: int, char_size: int) -> bool:
    """
    returns the result of the QLCS-test on string `x`, `y` of size `n` for substring of length `d`
    on characters of size `char_size`
    """
    # trivial cases
    if d == 0:
        return True
    if d == n:
        return x == y

    return grover_iter_search(x, y, n, d, char_size)


def quantum_LCS(x: str, y: str, n: int, char_size: int, verbose: bool = False) -> int:
    """
    resolves the longest common substring problem on strings `x`, `y`,
    of length `n`, with characters of `char_size` bits.
    When `verbose=True`, every test result is printed.
    """
    l = 0
    r = n
    while l < r:
        d = (l + r + 1) // 2
        if verbose:
            print(f"testing for common substring of length {d}:", end=" ")

        if quantum_LCS_test(x, y, n, d, char_size):
            l = d
            if verbose:
                print("found!")
        else:
            r = d - 1
            if verbose:
                print("not found!")
    return l


"""
an alphabet of size n uses characters made by
n binary digits. 

n = 2 => alphabet = {'00', '01', '10', '11'}
n = 3 => alphabet = {'000', '001', '010', '011',  
                     '100', '101', '110', '111'}

X and Y have to be of a size which is a power of 2, 
and always end with disjointed closing characters.

!!! IMPORTANT !!!
you have to keep any 2 characters from your alphabet as suffix
of each input. not following will lead to deterministically incorrect results.

on our examples, we use '00' and '11' or '000' and '111' as closing characters
"""


def classical_lcs_substring_len(x: str, y: str, char_size: int = 2) -> int:
    """
    for testing purpose
    """
    xs = [x[i : i + char_size] for i in range(0, len(x), char_size)]
    ys = [y[i : i + char_size] for i in range(0, len(y), char_size)]

    n = len(xs)
    m = len(ys)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    best = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if xs[i - 1] == ys[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
                best = max(best, dp[i][j])
    return best
