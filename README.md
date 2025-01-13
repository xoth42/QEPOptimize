## Quantum Entanglement Purification Optimizer


This repository is intended to organize and reform the code from QuantumHardware and qevo_optimizer,and qevo, while being readable, extendable, and abstracted as much as possible for ease of later advancements. 

# Building Plan:

The optimizer is run using two parts: 
1.  QEPO: The main runner for the genetic optimizer
2.  Configurable: this manages all parameters needed to run the optimizer, and has methods to set up various types of optimizers on various types of hardware



Then, a UI will be created that aims to meet these goals:

1. Simple interface, with no unnecessary user-chosen parameters
2. Clear visuals and data representation
3. Ability to Choose an IBM quantum computer, and its hardware specifications are loaded automatically
4. If desired by the user, let user see and edit all parameters before running

## From YipiaoWu/QuantumHardware:

 Base functionality for optimizing circuits using genetic algorithms
and an efficient representation of purification circuits
adapted to IBM quantum computer 

SUMMARY OF CONSTRAINTS AND STRATEGIES

1. **Qubit Reordering and Long-Range Entanglement**:
   - In quantum systems where qubits are physically distant, long-range entanglement must be established through intermediate steps.
   - The `NoisyBellSwap` gate is used to **reorder qubits** and move quantum information between distant qubits by SWAP operations, ensuring qubits can interact when needed for entanglement purification.
   - SWAP gates introduce noise, so their use must be minimized or strategically placed to reduce fidelity loss.

2. **Gate Connectivity Constraints**:
   - In IBM superconducting quantum computers, qubit connectivity is limited
    to nearest-neighbor interactions.
   - Only **adjacent qubits** can directly interact using two-qubit gates (e.g., CNOT gates). For non-adjacent qubits, operations like **SWAP gates** must be used to bring them into proximity.
   - The optimizer ensures that two-qubit operations (like `PauliNoiseBellGate` or `CNOTPerm`) only occur between qubit pairs that are physically connected.

3. **Noise Models and Calibration Data**:
   - Realistic noise models are incorporated based on **IBM calibration data** for the specific hardware used. These include:
     - **Two-qubit gate errors**: Each gate (especially CNOT and SWAP gates) introduces a level of noise depending on hardware calibration.
     - **Readout errors**: Measurement operations are noisy, with varying error rates based on qubit performance.
   - These noise parameters (e.g., T1, T2, and gate time) are pulled directly from the calibration data, ensuring that circuit performance predictions closely match real hardware behavior.

4. **Circuit Depth and Latency Constraints**:
   - Due to **cumulative noise and decoherence** over time the depth of a quantum circuit (i.e., the number of sequential operations) affects its reliability .
   - Higher-depth circuits are more prone to errors, especially in systems with short **coherence times (T1/T2)**.
   - Certain operations must consider their **gate durations** to mitigate these effects, such as reducing the number of sequential gates or parallelizing operations where possible.
   - The optimizer strives to minimize circuit depth while maintaining performance, striking a balance between gate fidelity and noise accumulation.

5. **Randomized Measurement Operations**:
   - Measurement operations on qubits introduce uncertainty due to **readout errors** and the state they leave the qubits in post-measurement.
   - Some measurement operations result in a **maximally mixed state**, especially when no reset mechanism is applied after measurement (`NoisyBellMeasureNoisyReset`).
   - The optimizer accounts for this uncertainty by considering **measurement noise** and ensuring that critical measurements are handled in a way that reduces the impact of these errors on overall circuit fidelity.

STRATEGY:
- This genetic algorithm framework evolves quantum circuits by mutating, selecting, and optimizing individuals, respecting the above constraints.
- It aims to find circuits that maximize the **fidelity of entanglement purification** while minimizing noise, ensuring the circuits are executable on noisy intermediate-scale quantum (NISQ) devices.

