# Set-Based Reinforcement Learning 

This folder contains all classes for set-based reinforcement learning [1].
Additionally, state-of-the-art adversarial methods [2] have been implemented to compare the performances. 

The folder contains the abstract agent class `RLagent`, of which a set/point-based DDPG `DDPGagent` [3] or TD3 `TD3agent` [4] can be instantiated. 
The agent consists of an actor `actor`, one or two critic(s) `critc`, and a replay buffer `buffer`. 
The control environment class `ctrlEnvironment` models the state transition.

Examples for set-based reinforcement learning are provided in `./cora/examples/nn/rl` for the Quadrocopter 1D benchmark `example_neuralNetwork_rl_RLagent_Quad1D.m'`and the Inverted Pendulum benchamrk. 
These examples are only trained for a single random seed and do not show the exact results as in [1], due to shorter training times. 

References:
- [1] Wendl, M. et al. 'Training Verifiably Robust Agents Using Set-Based Reinforcement Learning', 2024
- [2] Pattanaik, A. et al. 'Robust Deep Reinforcement Learning with Adversarial Attacks', Int. Conf. on Autonomous Agents and Multiagent Systems (AAMAS) 2018 
- [3] Lillicrap, T. et al. 'Continuous control with deep reinforcement learning', Int. Conf. on Learning Representations (ICLR), 2016 
- [4] Fujimoto, S. et al. 'Addressing Function Approximation Error in Actor-Critic Methods', Int. Conf. on Machine Learning (ICML), 2018

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>
