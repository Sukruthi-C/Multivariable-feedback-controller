
## Project Overview
This repository contains the final project for EECS 565, focusing on the control design for the Reactive Ion Etching (RIE) process. The project aims to improve the performance of the RIE process by employing a multivariable feedback controller. The design approach relies on using a Linear Quadratic Regulator (LQR) and Loop Transfer Recovery (LTR).

### Multivariable Feedback Controller Design

The project addresses the challenges associated with controlling a plant with strong interactions using a decentralized controller. A multivariable controller was designed to achieve the following specifications:
- Integral control in both loops to achieve zero steady state error in response to step commands.
- Minimized response to noise while maintaining fast response to step commands.
- Prevented throttle overshoot beyond 5%.
- Kept the response to step commands in \( V_{bias} \) small.

#### State Feedback Controller Design
- Utilized LQR with augmented integrators to choose feedback gains.
- Penalized both plant states and augmented integrator states to optimize performance.
- Achieved better step responses compared to decentralized control.

#### Observer Design Using LTR
- Designed an observer to approximate the state feedback design.
- Adjusted observer parameters to balance state feedback margins and sensor noise response.

### Multivariable Stability Margins
- Assessed stability margins for the state feedback design and compared them with the observer-based design.
- Evaluated singular values of loop transfer functions, sensitivity functions, and complementary sensitivity functions.

### Reverse Engineering the Multivariable Controller
- Analyzed the structure and transfer functions of the multivariable controller.
- Designed a controller using equivalent compensator structures for better understanding and performance.

### Additional Actuator Integration
- Introduced Oxygen (\( O_2 \)) as an additional actuator to regulate all three plasma variables: \( [F] \), \( V_{bias} \), and Pressure.
- Identified a linear model for the plant with the new actuator and designed a MIMO controller using LQG techniques.

## Project Structure

- `src/`: Contains the source code for the controller design and simulation.
- `data/`: Includes datasets and identified models used in the project.
- `plots/`: Generated plots for step responses, Bode plots, and stability margins.

## Getting Started

### Prerequisites
- MATLAB or Octave for running the simulations and control design.
- Control Systems Toolbox.

### Running the Simulations
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/rie-process-control.git
   cd Multivariable-feedback-controller
2. Run the main script:
   ```
   Part1and2.m

